
#' General imputation framework.
#'
#' @inheritParams vruv4
#' @param impute_func A function that takes as input a matrix names
#'     \code{Y} that has missing values and returns a matrix called
#'     \code{Yhat} of the same dimension of \code{Y} with the missing
#'     values filled in. If \code{do_variance = TRUE}, then
#'     \code{impute_func} should also return \code{sig_diag} --- a
#'     vector of column-specific variance estimates. The default is a
#'     wrapper for \code{\link[softImpute]{softImpute}}.  I provide a
#'     few functions in this package. \code{\link{knn_wrapper}}
#'     performs k-nearest-neighbors
#'     imputation. \code{\link{missforest_wrapper}} performs random
#'     forest imputation. \code{\link{softimpute_wrapper}} performs
#'     nuclear norm minimization imputation.
#' @param impute_args A list of additional parameters to pass to
#'     \code{impute_func}.
#' @param do_variance A logical. Does \code{impute_func} also return
#'     estimates of the column-specific variances?
#' @param k The rank of the underlying matrix. Used by
#'     \code{\link{hard_impute}} if that is the value of
#'     \code{impute_func}. If not provided, will be estimated
#'     by \code{\link[sva]{num.sv}}.
#' @return \code{beta2hat} The estimates of the coefficients of the
#'     covariates of interest that do not correspond to control genes.
#'
#'     \code{betahat_long} The estimates of the coefficients. Those
#'     corresponding to control genes are set to 0.
#'
#'     \code{sebetahat} If \code{do_variance = TRUE}, then these are
#'     the "standard errors" of \code{beta2hat} (but not really).
#'
#'     \code{tstats} If \code{do_variance = TRUE}, then these are
#'     the "t-statistics" of \code{beta2hat} (but not really).
#'
#'     \code{pvalues} If \code{do_variance = TRUE}, then these are
#'     the "p-values" of \code{tstats} (but not really).
#'
#'
#' @author David Gerard
#'
#' @export
ruvimpute <- function(Y, X, ctl, k = NULL, impute_func = em_miss,
                      impute_args = list(), cov_of_interest = ncol(X),
                      include_intercept = TRUE, do_variance = FALSE) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.numeric(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::assert_that(is.numeric(X))
    assertthat::assert_that(is.vector(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(all(cov_of_interest >= 1 & cov_of_interest <= ncol(X)))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.function(impute_func))
    assertthat::assert_that(is.list(impute_args))
    assertthat::assert_that(is.null(impute_args$Y21))
    assertthat::assert_that(is.null(impute_args$Y31))
    assertthat::assert_that(is.null(impute_args$Y32))
    assertthat::assert_that(is.logical(do_variance))


    rotate_out <- rotate_model(Y = Y, X = X, k = k, cov_of_interest =
                               cov_of_interest, include_intercept =
                               include_intercept, limmashrink = FALSE,
                               do_factor = FALSE)

    k <- rotate_out$k

    Y21 <- rotate_out$Y2[, ctl, drop = FALSE]
    Y22 <- rotate_out$Y2[, !ctl, drop = FALSE]
    Y31 <- rotate_out$Y3[, ctl, drop = FALSE]
    Y32 <- rotate_out$Y3[, !ctl, drop = FALSE]
    R22 <- rotate_out$R22

    if (identical(impute_func, em_miss)) {
        impute_args$Y21 <- Y21
        impute_args$Y31 <- Y31
        impute_args$Y32 <- Y32
        impute_args$k   <- k
        impute_args$gls <- TRUE
        emout <- do.call(what = em_miss, args = impute_args)
        Y22hat <- emout$Y22hat
    } else if (identical(impute_func, hard_impute)) {
        impute_args$Y21 <- Y21
        impute_args$Y31 <- Y31
        impute_args$Y32 <- Y32
        impute_args$k   <- k
        Y22hat <- do.call(what = hard_impute, args = impute_args)
    } else if (identical(impute_func, impute_ruv_reproduce) ) {
        impute_args$Y21 <- Y21
        impute_args$Y31 <- Y31
        impute_args$Y32 <- Y32
        impute_args$k   <- k
        Y22hat <- do.call(what = impute_ruv_reproduce, args = impute_args)
    } else {
        impout <- impute_block(Y21 = Y21, Y31 = Y31, Y32 = Y32, impute_func = impute_func,
                               impute_args = impute_args, do_variance = do_variance)
        Y22hat   <- impout$Y22hat
        sig_diag <- impout$sig_diag
    }

    R22inv <- backsolve(R22, diag(nrow(R22)))
    beta2hat <- R22inv %*% (Y22 - Y22hat)

    return_list <- list()
    return_list$beta2hat <- beta2hat

    betahat_long <- matrix(0, nrow = nrow(beta2hat), ncol = ncol(Y))
    betahat_long[, !ctl] <- beta2hat
    return_list$betahat_long <- betahat_long


    if (do_variance) {
        mult_matrix <- solve(t(X) %*% X)[cov_of_interest, cov_of_interest, drop = FALSE]
        sebetahat <- sqrt(outer(diag(mult_matrix), sig_diag))
        tstats <- beta2hat / sebetahat
        pvalues <- 2 * (stats::pt(q = -abs(tstats), df = nrow(X) - ncol(X)))

        return_list$sebetahat <- sebetahat
        return_list$tstats    <- tstats
        return_list$pvalues   <- pvalues
    }

    return(return_list)
}


#' Same as ruvimpute but only does ruvem and comes up with estimates of standard errors.
#'
#' This is more friendly than \code{\link{ruvimpute}}.
#'
#' @inheritParams ruvimpute
#' @param gls A logical. Should we use GLS (\code{TRUE}) or OLS
#'     (\code{FALSE})?
#' @author David Gerard
#' @export
#' @seealso \code{\link{ruvimpute}}, \code{\link{em_miss}}.
ruvem <- function(Y, X, ctl, k = NULL, impute_args = list(),
                  cov_of_interest = ncol(X), include_intercept = TRUE,
                  gls = TRUE) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.numeric(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::assert_that(is.numeric(X))
    assertthat::assert_that(is.vector(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(all(cov_of_interest >= 1 & cov_of_interest <= ncol(X)))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.list(impute_args))
    assertthat::assert_that(is.null(impute_args$Y21))
    assertthat::assert_that(is.null(impute_args$Y31))
    assertthat::assert_that(is.null(impute_args$Y32))


    rotate_out <- rotate_model(Y = Y, X = X, k = k, cov_of_interest =
                               cov_of_interest, include_intercept =
                               include_intercept, limmashrink = FALSE,
                               do_factor = FALSE)

    k <- rotate_out$k
    ncontrols <- sum(ctl)

    Y21 <- rotate_out$Y2[, ctl, drop = FALSE]
    Y22 <- rotate_out$Y2[, !ctl, drop = FALSE]
    Y31 <- rotate_out$Y3[, ctl, drop = FALSE]
    Y32 <- rotate_out$Y3[, !ctl, drop = FALSE]
    R22 <- rotate_out$R22

    impute_args$Y21 <- Y21
    impute_args$Y31 <- Y31
    impute_args$Y32 <- Y32
    impute_args$k   <- k
    impute_args$gls <- TRUE
    emout <- do.call(what = em_miss, args = impute_args)
    Y22hat <- emout$Y22hat

    R22inv <- backsolve(R22, diag(nrow(R22)))
    beta2hat <- R22inv %*% (Y22 - Y22hat)

    return_list <- list()
    return_list$beta2hat <- beta2hat
    betahat_long <- matrix(0, nrow = nrow(beta2hat), ncol = ncol(Y))
    betahat_long[, !ctl] <- beta2hat
    return_list$betahat_long <- betahat_long

    alphahat_final <- matrix(NA, nrow = nrow(emout$alpha), ncol = ncol(emout$alpha))
    alphahat_final[, ctl] <- emout$alpha[, 1:ncontrols]
    alphahat_final[, !ctl] <- emout$alpha[, (ncontrols + 1):ncol(emout$alpha)]
    sig_diag_final <- rep(NA, length = length(emout$sig_diag))
    sig_diag_final[ctl] <- emout$sig_diag[1:ncontrols]
    sig_diag_final[!ctl] <- emout$sig_diag[(ncontrols + 1):length(emout$sig_diag)]

    ## betahat <- R22inv %*%
    ##     (rotate_out$Y2 - emout$Z[1:length(cov_of_interest), , drop = FALSE] %*% alphahat_final)

    ## Get Z1 hat ------------------------------------------------
    Y1  <- rotate_out$Y1
    if (!is.null(Y1)) {
        R12 <- rotate_out$R12
        R11 <- rotate_out$R11
        Q   <- rotate_out$Q
        beta1_ols <- solve(R11) %*% (Y1 - R12 %*% rotate_out$betahat_ols)
        resid_top <- Y1 - R12 %*% betahat_long - R11 %*% beta1_ols
        if (gls) {
            Z1  <- solve(alphahat_final %*% diag(1 / sig_diag_final) %*% t(alphahat_final)) %*%
                alphahat_final %*% diag(1 / sig_diag_final) %*% t(resid_top)
        } else {
            Z1  <- solve(alphahat_final %*%  t(alphahat_final)) %*% alphahat_final %*% t(resid_top)
        }
        Zhat <- Q %*% rbind(t(Z1), emout$Z)
    } else {
        Q   <- rotate_out$Q
        Zhat <- Q %*% emout$Z
    }

    degrees_freedom   <- nrow(X) - ncol(X) - k
    lmout             <- limma::squeezeVar(sig_diag_final[!ctl], df = degrees_freedom)
    df_limma          <- degrees_freedom + lmout$df.prior
    sigma_limma       <- rep(NA, length = length(sig_diag_final))
    sigma_limma[!ctl] <- lmout$var.post

    XZ                      <- cbind(X, Zhat)
    mult_matrix             <- solve(t(XZ) %*% XZ)[cov_of_interest, cov_of_interest, drop = FALSE]
    sebetahat               <- t(sqrt(outer(sig_diag_final, diag(mult_matrix), FUN = "*")))
    sebetahat[, ctl]        <- NA
    sebetahat_limma         <- t(sqrt(outer(sigma_limma, diag(mult_matrix), FUN = "*")))
    sebetahat_limma[, ctl]  <- NA
    tstats          <- betahat_long / sebetahat
    tstats_limma    <- betahat_long / sebetahat_limma
    pvalues         <- 2 * (stats::pt(q = -abs(tstats), df = degrees_freedom))
    pvalues_limma   <- 2 * (stats::pt(q = -abs(tstats_limma), df = df_limma))

    return_list$alphahat        <- alphahat_final
    return_list$Zhat            <- Zhat
    return_list$sebetahat       <- sebetahat
    return_list$tstats          <- tstats
    return_list$pvalues         <- pvalues
    return_list$sebetahat_limma <- sebetahat_limma
    return_list$tstats_limma    <- tstats_limma
    return_list$pvalues_limma   <- pvalues_limma
    return_list$XZ              <- XZ
    return_list$mult_matrix     <- mult_matrix
    return_list$df              <- degrees_freedom
    return_list$df_limma        <- df_limma

    return(return_list)
}


#' Constructs an overall matrix, then applies a given imputation function.
#'
#' @param Y21 A matrix of size k by m
#' @param Y31 A matrix of size (n - k) by m
#' @param Y32 A matrix of size (n - k) by (p - m)
#' @inheritParams ruvimpute
#'
#' @return \code{Y22} A matrix of size k by (p - m). This is the
#'     imputed submatrix.
#'
#'     \code{sig_diag} If \code{do_variance = FALSE}, then this is
#'     \code{NULL}. Else, it is a vector of column-specific variance
#'     estimates.
#'
#' @author David Gerard
impute_block <- function(Y21, Y31, Y32, impute_func,
                         impute_args = list(), do_variance = FALSE) {
    assertthat::assert_that(is.matrix(Y21))
    assertthat::assert_that(is.matrix(Y31))
    assertthat::assert_that(is.matrix(Y32))
    assertthat::are_equal(nrow(Y31), nrow(Y32))
    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::assert_that(is.function(impute_func))
    assertthat::assert_that(is.list(impute_args))
    assertthat::assert_that(is.null(impute_args$Y))
    assertthat::assert_that(is.logical(do_variance))

    p <- ncol(Y31) + ncol(Y32)
    n <- nrow(Y21) + nrow(Y31)
    m <- ncol(Y21)
    k <- nrow(Y21)

    Y                       <- matrix(NA, nrow = n, ncol = p)
    Y[1:k, 1:m]             <- Y21
    Y[(k + 1):n, 1:m]       <- Y31
    Y[(k + 1):n, (m + 1):p] <- Y32

    impute_args$Y <- Y
    impout <- do.call(what = impute_func, args = impute_args)

    if (do_variance) {
        Yhat <- impout$Yhat
        sig_diag <- impout$sig_diag
    } else {
        Yhat <- impout
        sig_diag <- NULL
    }

    assertthat::assert_that(is.matrix(Yhat))
    assertthat::are_equal(dim(Yhat), dim(Y))

    Y22hat <- Yhat[1:k, (m + 1):p, drop = FALSE]
    return(list(Y22hat = Y22hat, sig_diag = sig_diag))
}






#' A wrapper for using the softImpute function from the softImpute package.
#'
#' @param Y A matrix with missing values.
#' @param max_rank The same as the \code{max.rank} option in
#'     \code{\link[softImpute]{softImpute}}.
#'
#' @return A matrix with the missing values imputed.
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso \code{\link[softImpute]{softImpute}}
softimpute_wrapper <- function(Y, max_rank) {
    ## lout <- softImpute::lambda0(x = Y)
    softout <- softImpute::softImpute(x = Y, rank.max = max_rank,
                                      lambda = 0, maxit = 1000,
                                      type = "svd")
    dbout <- softImpute::deBias(x = Y, svdObject = softout)
    cout <- softImpute::complete(x = Y, object = dbout)
    return(cout)
}

#' Wrapper for missForest package.
#'
#' @param Y A matirx with missing values.
#'
#' @return A matrix with the missing values imputed.
#'
#' @seealso \code{\link[missForest]{missForest}}
#'
#' @author David Gerard
#'
#' @export
missforest_wrapper <- function(Y) {
    trash <- utils::capture.output(impout <- missForest::missForest(xmis = Y))
    return(impout$ximp)
}

#' Wrapper for impute.knn
#'
#' @param Y A matirx with missing values.
#'
#' @return A matrix with the missing values imputed.
#'
#' @seealso \code{\link[impute]{impute.knn}}
#'
#' @author David Gerard
#'
#' @export
knn_wrapper <- function(Y) {
    impout <- impute::impute.knn(data = Y, colmax = 1, rowmax = 1)
    return(impout$data)
}



#' My version of hard imputation that begins at the ruv estimates
#'
#' @author David Gerard
#'
#' @param Y21 Top left of matrix.
#' @param Y31 Bottom left of matrix.
#' @param Y32 Top right of matrix.
#' @param tol The tolerance for stopping.
#' @param maxit The maximum number of iterations to run.
#' @param init_type Which version of RUV should we start?
#'     \code{"ruv3"} or \code{"ruv4"}?
#' @param k The rank of the mean matrix.
#'
#' @return The top right of the matrix.
#'
#' @export
hard_impute <- function(Y21, Y31, Y32, k, tol = 10 ^ -5, maxit = 1000,
                           init_type = c("ruv4", "ruv3")) {

    init_type <- match.arg(init_type)
    m <- ncol(Y21)
    num_cov <- nrow(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)
    degrees_freedom <- nrow(Y31) - k


    ## Factor analysis on Y31 ------------------------------------------------

    if (init_type == "ruv3") {
        pcout <- pca_naive(Y = Y31, r = k)
        alpha1 <- t(pcout$alpha)
        Z3 <- pcout$Z
        sig_diag1 <- pcout$sig_diag

        ## Regression to get alpha2 ----------------------------------------------
        alpha2 <- solve(t(Z3) %*% Z3) %*% t(Z3) %*% Y32
        sig_diag2 <- colSums((Y32 - Z3 %*% alpha2) ^ 2) / degrees_freedom
    } else if (init_type == "ruv4") {
        pcout <- pca_naive(cbind(Y31, Y32), r = k)
        alpha1 <- t(pcout$alpha[1:m, , drop = FALSE])
        alpha2 <- t(pcout$alpha[(m + 1):p, , drop = FALSE])
        Z3 <- pcout$Z
        sig_diag1 <- pcout$sig_diag[1:m]
    }

    ## Regression to get Z2 --------------------------------------------------
    limma_out1 <- limma::squeezeVar(var = sig_diag1,
                                    df = degrees_freedom)
    sig_diag1_temp <- limma_out1$var.post
    Z2 <- Y21 %*% diag(1 / sig_diag1_temp) %*% t(alpha1) %*%
        solve(alpha1 %*% diag(1 / sig_diag1_temp) %*% t(alpha1))

    ## get inital values for Y22
    Y22init <- Z2 %*% alpha2
    Y22new <- Y22init

    Y2 <- cbind(Y21, Y22init)
    Y3 <- cbind(Y31, Y32)
    Ycomp <- rbind(Y2, Y3)

    ismiss <- matrix(FALSE, nrow = n, ncol = p)
    ismiss[1:num_cov, (m + 1):p] <- TRUE

    iter_index <- 1
    err <- tol + 1
    while(err > tol & iter_index < maxit) {
        Y22old <- Y22new
        svout <- irlba::irlba(Ycomp, nv = k)
        ## svout$u %*% diag(x = svout$d, nrow = length(svout$d), ncol = length(svout$d)) %*%
        ##     t(svout$v)
        ## lowrank_est <- tcrossprod(sweep(svout$u, 2, svout$d, `*`), svout$v)
        Y22new <- tcrossprod(sweep(svout$u[1:num_cov, , drop = FALSE], 2, svout$d, `*`),
                             svout$v[(m + 1):p, , drop = FALSE])
        Ycomp[1:num_cov, (m + 1):p] <- Y22new
        err <- sum((Y22old - Y22new) ^ 2)
        iter_index <- iter_index + 1
        ## cat(err, "\n")
    }
    return(Y22new)
}


#' Reproduce RUV2, RUV3, and RUV4 with RUVimpute.
#'
#' This is a proof of concept function to inegrate
#' \code{\link{ruvimpute}} with \code{\link[ruv]{RUV2}},
#' \code{\link{ruv3}}, \code{\link[ruv]{RUV4}}.
#'
#' @inheritParams hard_impute
#' @param impute_type Which version of RUV should we reproduce?
#'     \code{"ruv2"}, \code{"ruv3"}, or \code{"ruv4"}?
#'
#' @export
#'
#' @author David Gerard
impute_ruv_reproduce <- function(Y21, Y31, Y32, k,
                                 impute_type = c("ruv2", "ruv3", "ruv4")) {
    impute_type <- match.arg(impute_type)
    m <- ncol(Y21)
    num_cov <- nrow(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)
    degrees_freedom <- nrow(Y31) - k

    if (impute_type == "ruv2") {
        pcout <- pca_naive(Y = rbind(Y21, Y31), r = k)
        Z3hat <- pcout$Z[(num_cov + 1):n, , drop = FALSE]
        Z2hat <- pcout$Z[1:num_cov, , drop = FALSE]
        alpha2hat <- solve(crossprod(Z3hat)) %*% crossprod(Z3hat, Y32)
        Y22hat <- Z2hat %*% alpha2hat
    } else if (impute_type == "ruv3") {
        pcout <- pca_naive(Y = Y31, r = k)
        Z3hat <- pcout$Z
        alpha1hat <- t(pcout$alpha)
        Z2hat <- tcrossprod(Y21, alpha1hat) %*% solve(tcrossprod(alpha1hat))
        alpha2hat <- solve(crossprod(Z3hat)) %*% crossprod(Z3hat, Y32)
        Y22hat <- Z2hat %*% alpha2hat
    } else if (impute_type == "ruv4") {
        pcout <- pca_naive(Y = cbind(Y31, Y32), r = k)
        alpha1hat <- t(pcout$alpha[1:m, , drop = FALSE])
        alpha2hat <- t(pcout$alpha[(m + 1):p, , drop = FALSE])
        Z2hat <- tcrossprod(Y21, alpha1hat) %*% solve(tcrossprod(alpha1hat))
        Y22hat <- Z2hat %*% alpha2hat
    }

    return(Y22hat)
}


## mice_wrapper <- function(Y) {
##     library(mice)
##     mout <- mice::mice(data = Y)
##     mcomplete <- mice::complete(mout)
##     return(mcomplete)
## }

## Throws error in CRAN checks when used.
## flashr_wrapper <- function(Y, max_rank) {
##     if (!requireNamespace("flashr", quietly = TRUE)) {
##         stop("Sorry, flashr needs to be installed to use flashr_wrapper.")
##     }
##     trash <- utils::capture.output(gout <- flashr::greedy(Y = Y, K = max_rank,
##                                                    flash_para = list(partype = "var_col")))
##     Yhat <- gout$l %*% t(gout$f)
##     return(Yhat)
## }
