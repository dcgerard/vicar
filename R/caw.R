## Confounder Adjustment with Weights

#' Iterative procedure for confounder correction with a procedure that
#' returns lfdrs.
#'
#' Quoth the Raven "Caw, caw!"
#'
#' @inheritParams vruv4
#' @param weight_func The function that returns the weights (or
#'     lfdr's). Many forms of input are allowed. See
#'     \code{weight_func_input} for details.
#' @param weight_args Additional arguments to pass to
#'     \code{weight_func}.
#' @param scale_var A logical. Should we scale the variance
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param weight_init A character. How should we initialize the
#'     weights? The options are to initialize in the all-null setting
#'     (\code{"all_null"}), draw the weights randomly form iid
#'     uniforms (\code{"random"}), or run an iteration of
#'     \code{weight_func} prior to the first round of estimating the
#'     confounders (\code{"limma"}). This last step uses limma-ebayes.
#' @param weight_func_input The form of input for
#'     \code{weight_func}. Right now only \code{"summary2"} is
#'     supported, but I intend to support all of the following in the
#'     future. If \code{weight_func_input = "summary1"} then the
#'     function only takes p-values as input (called
#'     \code{pvalues}). If \code{weight_func_input = "summary2"}, then
#'     the function only takes a vector of effect estimates
#'     \code{betahat}, a vector of standard errors \code{sebetahat},
#'     and a vector of degrees of freedom \code{degrees_freedom}. If
#'     \code{weight = "summary3"}, then the input is a matrix of
#'     effects \code{betamat}, an array of covariances
#'     \code{cov_array} where the each \code{cov_array[,,i]} is the
#'     covariance of the elements of \code{betamat[i, ]}, and a vector
#'     of degrees of freedom \code{dfvec}. If \code{weight_func_input
#'     = "full"}, then the input is just a response matrix \code{Y}
#'     and a covariate matrix \code{X}.
#' @param min_scale The minimum estimate for the variance inflation term.
#'
#' @author David Gerard
#'
caw <- function(Y, X, k = NULL,
                cov_of_interest = ncol(X),
                limmashrink = TRUE,
                weight_func = ash_wrap,
                weight_args = list(),
                fa_func = pca_naive,
                fa_args = list(),
                scale_var = TRUE,
                include_intercept = TRUE,
                weight_init = c("all_null", "random", "limma"),
                weight_func_input = c("summary2"),
                degrees_freedom = NULL,
                min_scale = 0.8) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(X), nrow(Y))
    assertthat::assert_that(all(abs(cov_of_interest - round(cov_of_interest)) < 10 ^ -14))
    assertthat::assert_that(all(cov_of_interest >= 1 & cov_of_interest <= ncol(X)))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y))
    assertthat::assert_that(is.null(fa_args$r))
    assertthat::assert_that(is.function(fa_func))
    assertthat::assert_that(is.logical(scale_var))
    assertthat::assert_that(is.function(weight_func))
    assertthat::assert_that(is.list(weight_args))

    if (length(cov_of_interest) > 1) {
        if (weight_func_input == "summary1" | weight_func_input == "summary2") {
            stop("the weight_func_input specified only allows for one covariate of interest.")
        }
    }

    weight_init <- match.arg(weight_init)

    ## only need to call this once
    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink, fa_func = fa_func,
                               fa_args = fa_args, do_factor = TRUE)
    k <- rotate_out$k

    ## Deal with degrees of freedom
    if (!is.null(rotate_out$prior_df) & is.null(degrees_freedom)) {
        degrees_freedom <- rotate_out$prior_df + nrow(X) - ncol(X) - k
    } else if (is.null(degrees_freedom)) {
        degrees_freedom <- nrow(X) - ncol(X) - k
    }
    assertthat::assert_that(length(degrees_freedom) == 1 | length(degrees_freedom) == ncol(Y))
    assertthat::assert_that(all(degrees_freedom > 0))
    rotate_out$degrees_freedom <- degrees_freedom


    p <- ncol(Y)

    if (weight_init == "all_null") {
        w_init <- rep(1, p)
    } else if (weight_init == "random") {
        w_init <- stats::runif(p)
    } else if (weight_init == "limma") {
        if (weight_func_input == "summary1") {
            lmout <- limma::lmFit(object = t(Y), design = X)
            ebout <- limma::eBayes(lmout)
            weight_args$pvalues <- ebout$p.value[, cov_of_interest, drop = TRUE]
            w_init <- do.call(weight_func, args = weight_args)
        } else if (weight_func_input == "summary2") {
            lmout <- limma::lmFit(object = t(Y), design = X)
            ebout <- limma::eBayes(lmout)
            weight_args_temp <- weight_args
            weight_args_temp$sebetahat <-
                sqrt(lmout$cov.coefficients[cov_of_interest, cov_of_interest,
                                            drop = TRUE] * ebout$s2.post)
            weight_args_temp$betahat <- lmout$coefficients[, cov_of_interest, drop = TRUE]
            weight_args_temp$degrees_freedom <- ebout$df.total
            w_init <- do.call(what = weight_func, args = weight_args_temp)
        } else if (weight_func_input == "summary3") {
            lmout <- limma::lmFit(object = t(Y), design = X)
            ebout <- limma::eBayes(lmout)
            weight_args$cov_array <- outer(lmout$cov.coefficients[cov_of_interest,
                                                                  cov_of_interest, drop = FALSE],
                                           ebout$s2.post)
            weight_args$dfvec   <- ebout$df.total
            weight_args$betamat <- lmout$coefficients[, cov_of_interest, drop = FALSE]
            w_init <- do.call(what = weight_func, args = weight_args)
        } else if (weight_func_input == "full") {
            weight_args$Y <- Y
            weight_args$X <- X
            w_init <- do.call(what = weight_func, args = weight_args)
        }
    }

    w_current <- w_init
    Ws <- w_current / rotate_out$sig_diag
    Z2hat <- rotate_out$Y2 %*% diag(Ws) %*% rotate_out$alpha %*%
        solve(crossprod(rotate_out$alpha, diag(Ws) %*% rotate_out$alpha))
    scale_val <- 1

    wzs_vec_init <- c(w_current, Z2hat, scale_val)
    sqout <- SQUAREM::squarem(par = wzs_vec_init, fixptfn = fix_caw_wrapper,
                              rotate_out = rotate_out,
                              weight_func = weight_func,
                              weight_args = weight_args, min_scale = min_scale)

    w_current <- sqout$par[1:p]
    Z2hat     <- matrix(sqout$par[(p + 1):(p + prod(dim(Z2hat)))], ncol = ncol(Z2hat))
    scale_val <- sqout$par[length(sqout$par)]

    ## one last run to get output ---------------------------------------------------
    fout <- fix_caw(w_current = w_current, Z2hat = Z2hat, scale_val = scale_val,
                    rotate_out = rotate_out, weight_func = weight_func,
                    weight_args = weight_args, min_scale = min_scale, return_val = "full")

    return(fout)
}

#' A wrapper for \code{\link{fix_caw}} so that I can use SQUAREM.
#'
#' @inheritParams fix_caw
#' @param wzs_vec A vector of numerics. The first p of which are the
#'     \code{w_current}, the next \code{length(cov_of_interest) * k}
#'     of which are \code{Z2hat}, the last of which is
#'     \code{scale_val}.
#'
#' @author David Gerard
fix_caw_wrapper <- function(wzs_vec, rotate_out,
                            weight_func = ash_wrap,
                            weight_args = list(), min_scale = 0.8) {
    ncov <- ncol(rotate_out$X)
    assertthat::are_equal(length(wzs_vec), ncov * k + p + 1)

    p    <- ncol(rotate_out$Y2)
    ncov <- nrow(rotate_out$Y2)
    k    <- ncol(rotate_out$alpha)
    w_current <- wzs_vec[1:p]
    Z2hat     <- matrix(wzs_vec[(p + 1):(ncov * k + p)], ncol = k)
    scale_val <- wzs_vec[length(wzs_vec)]
    fout <- fix_caw(w_current = w_current, Z2hat = Z2hat, scale_val = scale_val,
                    rotate_out = rotate_out, weight_func = weight_func,
                    weight_args = weight_args, min_scale = min_scale)
    wzs_vec_new <- c(fout$w_current, fout$Z2hat, fout$scale_val)
    return(wzs_vec_new)
}

#' A fixed point iteration for CAW.
#'
#' @param w_current A vector of weights. The probability given Z2hat
#'     that each coefficient is null. These should all be between 0
#'     and 1.
#' @param Z2hat A matrix. The current estimates of the unobserved
#'     confounders.
#' @param scale_val A positive numeric. The variance inflation term.
#' @param rotate_out The output from \code{\link{rotate_model}} with
#'     \code{degrees_freedom} added to the list.
#' @param weight_func The function that calculates the weights given
#'     \code{w_current}. So far, this only supports functions that
#'     take as input three quantities: \code{betahat},
#'     \code{sebetahat}, and \code{degrees_freedom}
#' @param weight_args Additional arguments to pass to
#'     \code{weight_func}.
#' @param min_scale A scalar. The minimum value that \code{scale_val}
#'     is allowed to take. In practice, you wouldn't think that
#'     shrinking the variance is a good idea. But a value of 1 seems
#'     to stringent. So I set this to 0.8 to give you some wiggle
#'     room. Most of the time \code{scale_val} is usually estimated
#'     greater than 1.
#' @param return_val Should we return just the parameters of interest
#'     (\code{"params"}), or the full shebang (\code{"full"})? Setting
#'     \code{return_val = "params"} is used mostly for the
#'     optimization and returns only \code{w_current}, \code{Z2hat},
#'     and \code{scale_val}.
#'
#' @author David Gerard
fix_caw <- function(w_current, Z2hat, scale_val, rotate_out,
                    weight_func = ash_wrap, weight_args = list(), min_scale = 0.8,
                    return_val = c("params", "full")) {

    ## check inpute
    assertthat::assert_that(all(w_current >= 0 & w_current <= 1))
    assertthat::assert_that(is.matrix(Z2hat))
    assertthat::are_equal(ncol(Z2hat), ncol(rotate_out$alpha))
    assertthat::assert_that(is.function(weight_func))
    assertthat::assert_that(is.list(weight_args))
    assertthat::assert_that(is.null(weight_args$betahat))
    assertthat::assert_that(is.null(weight_args$sebetahat))
    assertthat::assert_that(is.null(weight_args$degrees_freedom))
    assertthat::assert_that(min_scale > 0)
    return_val <- match.arg(return_val)

    ## Update Z2hat
    Ws <- w_current / rotate_out$sig_diag
    Z2hat <- rotate_out$Y2 %*% diag(Ws) %*% rotate_out$alpha %*%
        solve(crossprod(rotate_out$alpha, diag(Ws) %*% rotate_out$alpha))

    ## get estimates
    diff_val <- (rotate_out$Y2 - tcrossprod(Z2hat, rotate_out$alpha))
    betahat <- solve(rotate_out$R22) %*% diff_val

    ## get scale val
    ## scale_val <- 1
    ## scale_val <- mean(diff_val ^ 2 %*% diag(Ws))
    scale_val <- max(mean(diff_val ^ 2 %*% diag(Ws)), min_scale)

    ## Generate Z
    if (!is.null(rotate_out$Y1)) {
        R12 <- rotate_out$R12
        R11 <- rotate_out$R11
        Q   <- rotate_out$Q
        beta1_ols <- solve(R11) %*% (rotate_out$Y1 - R12 %*% rotate_out$betahat_ols)
        resid_top <- rotate_out$Y1 - R12 %*% betahat - R11 %*% beta1_ols
        ## resid_top2 <- rotate_out$R12 %*% (betahat - solve(rotate_out$R22) %*% rotate_out$Y2)
        Z1  <- solve(t(rotate_out$alpha) %*% diag(1 / rotate_out$sig_diag) %*%
                     rotate_out$alpha) %*%
            t(rotate_out$alpha) %*% diag(1 / rotate_out$sig_diag) %*% t(resid_top)
        Zhat <- Q %*% rbind(t(Z1), Z2hat, rotate_out$Z3)
    } else {
        Q   <- rotate_out$Q
        Zhat <- Q %*% rbind(Z2hat, rotate_out$Z3)
    }

    XZ <- cbind(rotate_out$X, Zhat)
    ncolx <- ncol(rotate_out$X)
    mult_matrix <- solve(t(XZ) %*% XZ)[ncolx, ncolx, drop = FALSE]

    ## new sebetahat
    sebetahat <- sqrt(scale_val * c(mult_matrix) * rotate_out$sig_diag)

    weight_args$betahat <- c(betahat)
    weight_args$sebetahat <- c(sebetahat)
    weight_args$degrees_freedom <- rotate_out$degrees_freedom

    w_current <- do.call(what = weight_func, args = weight_args)

    if (return_val == "params") {
        return(list(w_current = w_current, Z2hat = Z2hat, scale_val = scale_val))
    } else if (return_val == "full") {
        return_list <- list()
        return_list$w_current       <- w_current
        return_list$Z2hat           <- Z2hat
        return_list$scale_val       <- scale_val
        return_list$betahat         <- betahat
        return_list$sebetahat       <- sebetahat
        return_list$Zhat            <- Zhat
        return_list$mult_matrix     <- mult_matrix
        return_list$degrees_freedom <- rotate_out$degrees_freedom

        weight_args$return_val <- "full"
        fout <- do.call(what = weight_func, args = weight_args)
        return_list$fout <- fout
        return(return_list)
    }
}

qvalue_wrap <- function(pvalues, return_val = c("params", "full")) {
    return_val <- match.arg(return_val)
    if (requireNamespace("qvalue", quietly = TRUE)) {
        if (return_val == "params") {
            return(qvalue::qvalue(p = pvalues)$lfdr)
        } else if (return_val == "full") {
            return(qvalue::qvalue(p = pvalues))
        }
    } else {
        stop("qvalue not installed")
    }
}

#' Wrapper for \code{\link[ashr]{ash.workhorse}} for use in \code{\link{caw}}.
#'
#' @param betahat The current effect size estimates
#' @param sebetahat The current standard error estimates.
#' @param degrees_freedom The degrees of freedom.
#' @param return_val Should we return just the lfdr's
#'     (\code{"params"}) or the entire output from
#'     \code{\link[ashr]{ash.workhorse}} (\code{"full"})?
#'
#' @author David Gerard
ash_wrap <- function(betahat, sebetahat, degrees_freedom,
                     return_val = c("params", "full")) {
    return_val <- match.arg(return_val)
    betahat <- as.vector(betahat)
    sebetahat <- as.vector(sebetahat)
    degrees_freedom <- as.vector(degrees_freedom)
    assertthat::are_equal(length(betahat), length(sebetahat))
    args <- list()
    args$betahat   <- betahat
    args$sebetahat <- sebetahat
    args$df        <- degrees_freedom
    args$outputlevel <- 2
    ashout <- do.call(what = ashr::ash.workhorse, args = args)
    if (return_val == "params") {
        return(ashout$result$lfdr)
    } else if (return_val == "full") {
        return(ashout)
    }
}

ash_arrays <- function(betamat, cov_array, dfvec) {
    assertthat::assert_that(is.matrix(betamat))
    assertthat::assert_that(is.array(cov_array))
    assertthat::are_equal(ncol(betamat), 1)
    assertthat::are_equal(dim(cov_array)[1], 1)
    assertthat::are_equal(dim(cov_array)[2], 1)
    assertthat::are_equal(nrow(betamat), dim(cov_array)[3])
    assertthat::are_equal(nrow(betamat), length(dfvec))
    assertthat::assert_that(is.list(args))

    args$betahat     <- c(betamat)
    args$sebetahat   <- c(cov_array)
    args$df          <- dfvec
    args$outputlevel <- 2
    ashout <- do.call(what = ashr::ash.workhorse, args = args)
    return(ashout$result$lfdr)
}


ash_full <- function(Y, X, cov_of_interest = ncol(X)) {
    lmout <- limma::lmFit(object = t(Y), design = X)
    ebout <- limma::eBayes(lmout)
    sebetahat <- lmout$cov.coefficients[cov_of_interest, cov_of_interest,
                                        drop = TRUE] * ebout$s2.post
    betahat <- lmout$coefficients[, cov_of_interest, drop = TRUE]
    df <- ebout$df.total
    ashout <- ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, df = df)
    return(ashout$result$lfdr)
}
