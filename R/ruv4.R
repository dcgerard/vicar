#' Calibrated RUV4 where the control genes are used to estimate hidden
#' confounders and variance inflation parameter.
#'
#' This function will perform a variant of Removing Unwanted Variation
#' 4-step (RUV4) (Gagnon-Bartsch et al, 2013) where the control genes
#' are used not only to estimate the hidden confounders, but to
#' estimate a variance inflation parameter. This variance inflation
#' step is akin to the "empirical null" approach of Efron (2004).
#'
#' The model is \deqn{Y = XB + ZA + E,} where \eqn{Y} is a matrix of
#' responses (e.g. log-transformed gene expression levels), \eqn{X} is
#' a matrix of covariates, \eqn{B} is a matrix of coefficients,
#' \eqn{Z} is a matrix of unobserved confounders, \eqn{A} is a matrix
#' of unobserved coefficients of the unobserved confounders, and
#' \eqn{E} is the noise matrix where the elements are independent
#' Gaussian and each column shares a common variance. The rows of
#' \eqn{Y} are the observations (e.g. individuals) and the columns of
#' \eqn{Y} are the response variables (e.g. genes).
#'
#' This model is fit using a two-step approach proposed in
#' Gagnon-Bartsch et al (2013) and described in Wang et al (2015),
#' modified to include estimating a variance inflation parameter. In
#' the current implementation, only the coefficients of one covariate
#' can be estimated using ASH. The rest are regressed out using OLS.
#'
#' @param Y A matrix of numerics. These are the response variables
#'     where each column has its own variance. In a gene expression
#'     study, the rows are the individuals and the columns are the
#'     genes.
#' @param X A matrix of numerics. The covariates of interest.
#' @param k A non-negative integer.The number of unobserved
#'     confounders. If not specified and the R package sva is
#'     installed, then this function will estimate the number of
#'     hidden confounders using the methods of Buja and Eyuboglu
#'     (1992).
#' @param cov_of_interest A vector of positive integers. The column
#'     numbers of the covariates in X whose coefficients you want to
#'     apply ASH to.
#' @param ctl A vector of logicals of length \code{ncol(Y)}. If
#'     position i is \code{TRUE} then position i is considered a
#'     negative control. If \code{ctl = NULL} (the default) then ASH
#'     will be run on the OLS estimates and corresponding standard
#'     errors.
#' @param include_intercept A logical. If \code{TRUE}, then it will
#'     check \code{X} to see if it has an intercept term. If not, then
#'     it will add an intercept term. If \code{FALSE}, then \code{X}
#'     will be unchanged.
#' @param gls A logical. Should we use generalized least squares
#'     (\code{TRUE}) or ordinary least squares (\code{FALSE}) for
#'     estimating the confounders? The OLS version is equivalent to
#'     using RUV4 to estimate the confounders.
#' @param likelihood Either \code{"normal"} or \code{"t"}. If
#'     \code{likelihood = "t"}, then the user may provide the degrees
#'     of freedom via \code{degrees_freedom}.
#' @param degrees_freedom if \code{likelihood = "t"}, then this is the
#'     user-defined degrees of freedom for that distribution. If
#'     \code{degrees_freedom} is \code{NULL} then the degrees of
#'     freedom will be the sample size minus the number of covariates
#'     minus \code{k}.
#' @param limmashrink A logical. Should we apply hierarchical
#'     shrinkage to the variances (\code{TRUE}) or not (\code{FALSE})?
#' @param fa_func A factor analysis function. The function must have
#'     as inputs a numeric matrix \code{Y} and a rank (numeric scalar)
#'     \code{r}. It must output a numeric matrix \code{alpha} and a
#'     numeric vector \code{sig_diag}. \code{alpha} is the estimate of
#'     the coefficients of the unobserved confounders, so it must be
#'     an \code{r} by \code{ncol(Y)} matrix. \code{sig_diag} is the
#'     estimate of the column-wise variances so it must be of length
#'     \code{ncol(Y)}. The default is the function \code{pca_naive}
#'     that just uses the first \code{r} singular vectors as the
#'     estimate of \code{alpha}. The estimated variances are just the
#'     column-wise mean square.
#' @param fa_args A list. Additional arguments you want to pass to
#'     fa_func.
#'
#'
#' @return A list whose elements are:
#'
#'     \code{multiplier} A numeric. The estimated variance inflation parameter.
#'
#'     \code{betahat_ols} A vector of numerics. The ordinary least
#'     squares estimates of the coefficients of the covariate of
#'     interest. This is when not including the estimated confounding
#'     variables.
#'
#'     \code{sebetahat_ols} A vector of positive numerics. The
#'     pre-inflation standard errors of \code{ruv$betahat} (NOT
#'     \code{ruv$betahat_ols}).
#'
#'     \code{betahat} A vector of numerics. The ordinary least squares
#'     estimates of the coefficients of the covariate of interest WHEN
#'     YOU ALSO INCLUDE THE ESTIMATES OF THE UNOBSERVED CONFOUNDERS.
#'
#'     \code{sebetahat} A vector of positive numerics. This is equal
#'     to sqrt(ruv$sebethat_ols * ruv$multiplier). This is the
#'     post-inflation adjusted standard errors for \code{ruv$betahat}.
#'
#'     \code{tstats} A vector of numerics. The t-statistics for
#'     testing against the null hypothesis of the coefficient of the
#'     covariate of interest being zero. This is after estimating the
#'     variance inflation parameter but before the posthoc-inflation.
#'
#'     \code{pvalues} A vector of numerics. The p-values of said test
#'     above.
#'
#'     \code{alphahat} A matrix of numerics. The estimates of the
#'     coefficients of the hidden confounders.
#'
#'     \code{sigma2} A vector of positive numerics. The estimates of
#'     the variances.
#'
#'     \code{fnorm_x} A numeric. This is the diagonal element of
#'     \code{t(X) \%*\% X} that corresponds to the covariate of
#'     interest. Returned mostly for debugging reasons and may be
#'     removed in the future.
#'
#'     \code{Z1} A matrix of numerics of length 1. This is the
#'     estimated confounders (after a rotation). Not useful on it's
#'     own and is mostly returned for debugging purposes. It may be
#'     removed in the future.
#'
#' @export
#'
#' @author David Gerard
#'
#' @references Gagnon-Bartsch, J., Laurent Jacob, and Terence
#'     P. Speed. "Removing unwanted variation from high dimensional
#'     data with negative controls."
#'     Berkeley: Department of Statistics. University of California
#'     (2013).
#'
#'     Andreas Buja and Nermin
#'     Eyuboglu. "Remarks on parallel analysis." Multivariate behavior
#'     research, 27(4):509-540, 1992.
#'
#'     Bradley Efron
#'     "Large-Scale Simultaneous Hypothesis Testing: The Choice of a Null
#'     Hypothesis",
#'     Journal of the American Statistical Association, 99:465,
#'     96-104, 2004.
#'
#'     Wang, J., Zhao, Q., Hastie, T., & Owen, A. B
#'     "Confounder Adjustment in Multiple Hypotheses Testing."
#'     arXiv preprint arXiv:1508.04178 (2015).
#'
#'
vicarius_ruv4 <- function(Y, X, ctl, k = NULL,
                          cov_of_interest = ncol(X),
                          likelihood = c("t", "normal"),
                          limmashrink = FALSE, degrees_freedom = NULL,
                          include_intercept = TRUE, gls = TRUE,
                          fa_func = pca_naive, fa_args = list()) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(all(cov_of_interest >= 1 & cov_of_interest <= ncol(X)))
    assertthat::assert_that(is.logical(gls))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y))
    assertthat::assert_that(is.null(fa_args$r))
    assertthat::assert_that(is.function(fa_func))

    likelihood <- match.arg(likelihood)
    if (likelihood == "normal" & !is.null(degrees_freedom)) {
        message("likelihood = \"normal\" but degrees_freedom not NULL. Ignoring degrees_freedom.")
    }
    if (is.null(degrees_freedom)) {
        degrees_freedom <- nrow(X) - ncol(X) - k
    }
    assertthat::assert_that(length(degrees_freedom) == 1 | length(degrees_freedom) == ncol(Y))
    assertthat::assert_that(all(degrees_freedom > 0))

    ## RUN THE ROTATED MODEL HERE -------------------------------------------
    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink, fa_func = fa_func,
                               fa_args = fa_args)

    betahat_ols     <- rotate_out$betahat_ols

    ## RUN RUV4 HERE ----------------------------------------------------------
    if (length(cov_of_interest) == 1) {
        alpha_scaled    <- rotate_out$alpha / rotate_out$Rsub[1, 1]
        sig_diag_scaled <- rotate_out$sig_diag / rotate_out$Rsub[1, 1] ^ 2
        k               <- rotate_out$k
        ruv4_out <- cruv4(betahat_ols = t(betahat_ols), alpha_scaled = alpha_scaled,
                          sig_diag_scaled = sig_diag_scaled, ctl = ctl,
                          degrees_freedom = degrees_freedom,
                          gls = gls, likelihood = likelihood)

        ruv4_out$alphahat <- rotate_out$alpha
        ruv4_out$sigma2   <- rotate_out$sig_diag
        ruv4_out$fnorm_x  <- rotate_out$Rsub[1, 1]
    } else {
        Y2       <- rotate_out$Y2
        alpha    <- rotate_out$alpha
        sig_diag <- rotate_out$sig_diag
        Rsub     <- rotate_out$Rsub
        ruv4_out <- cruv4_multicov(Y2 = t(Y2), alpha = alpha,
                                   sig_diag = sig_diag, ctl = ctl,
                                   Rsub = Rsub,
                                   degrees_freedom = degrees_freedom,
                                   gls = gls, likelihood = likelihood)
    }

    return(ruv4_out)
}

cruv4_multicov <- function(Y2, alpha, sig_diag, ctl, Rsub, degrees_freedom,
                           gls = TRUE, likelihood = "t") {
    assertthat::assert_that(is.matrix(Y2))
    assertthat::assert_that(is.matrix(alpha))
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::assert_that(is.matrix(Rsub))
    assertthat::assert_that(all(sig_diag > 0))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(length(degrees_freedom) == 1 | length(degrees_freedom) == ncol(Y2))
    assertthat::assert_that(all(degrees_freedom > 0))
    assertthat::are_equal(nrow(Y2), nrow(alpha))

    k <- ncol(alpha)

    ## Use control genes to jointly estimate Z1 and variance scaling parameter.
    Yc <- Y2[ctl, , drop = FALSE]
    if (k != 0) {
        alphac       <- alpha[ctl, , drop = FALSE]
        sig_diag_inv <- 1 / sig_diag[ctl]
        if (gls) {
            Z1 <- crossprod(solve(crossprod(alphac, sig_diag_inv * alphac)),
                            crossprod(alphac, sig_diag_inv * Yc))
        } else {
            Z1 <- crossprod(solve(crossprod(alphac, alphac)),
                            crossprod(alphac, Yc))
        }
        resid_mat <- Yc - alphac %*% Z1
        betahat <- (Y2 - alpha %*% Z1) %*% solve(t(Rsub))
    } else {
        Z1        <- NULL
        resid_mat <- Yc
        betahat   <- Y2 %*% solve(Rsub)
    }

    ## Gaussian MLE of variance inflation parameter.
    multiplier <- mean(resid_mat ^ 2 / sig_diag[ctl])

    ## T-likelihood regression and estimate variance inflation
    ## parameter using control genes.
    if (likelihood == "t") {
        if (k != 0) {
            if (!gls) {
                message("gls = FALSE not supported for t-likelihood, using gls = TRUE.")
            }
            alphac <- alpha[ctl, , drop = FALSE]
            tout <- tregress_em(Y = Yc, alpha = alphac,
                                sig_diag = sig_diag[ctl],
                                nu = degrees_freedom, lambda_init = multiplier,
                                Z_init = Z1)
            Z1         <- tout$Z
            betahat    <- (Y2 - alpha %*% Z1) %*% solve(Rsub)
            multiplier <- tout$lambda
        } else if (k == 0) {
            tout <- stats::optim(par = multiplier, fn = tregress_obj_wrapper,
                                 Z = matrix(0, nrow = 1, ncol = ncol(Yc)), Y = Yc,
                                 alpha = matrix(0, nrow = nrow(Yc), ncol = 1),
                                 sig_diag = sig_diag[ctl],
                                 nu = degrees_freedom, method = "Brent",
                                 lower = 0, upper = 10,
                                 control = list(fnscale = -1))
            multiplier <- tout$par
        }
    }

    ## Output values
    ruv4_out <- list()
    mult_matrix       <- solve(Rsub %*% t(Rsub))
    sig_diag_adjusted <- sig_diag * multiplier
    sebetahat         <- sqrt(outer(sig_diag_adjusted, diag(mult_matrix), FUN = "*"))
    sebetahat_ols     <- sqrt(outer(sig_diag, diag(mult_matrix), FUN = "*"))
    betahat_ols       <- Y2 %*% solve(t(Rsub))


    ruv4_out$sigma2            <- sig_diag
    ruv4_out$multiplier        <- multiplier
    ruv4_out$sigma2_adjusted   <- sig_diag_adjusted
    ruv4_out$mult_matrix       <- mult_matrix
    ruv4_out$betahat_ols       <- betahat_ols
    ruv4_out$sebetahat_ols     <- sebetahat_ols
    ruv4_out$betahat           <- betahat
    ruv4_out$sebetahat         <- sebetahat
    ruv4_out$tstats            <- betahat / sebetahat
    ruv4_out$pvalues           <- 2 * (stats::pt(q = -abs(ruv4_out$tstats),
                                                 df = degrees_freedom))
    ruv4_out$Z1                <- Z1
    ruv4_out$alphahat          <- alpha
    ruv4_out$Rsub              <- Rsub
    return(ruv4_out)
}


#' RUV4's second step.
#'
#' @param betahat_ols A matrix of numerics with one column. The ols
#'     estimates of beta.
#' @param alpha_scaled A matrix of numerics. The estimated
#'     coefficients of the unobserved confounders.
#' @param sig_diag_scaled A vector of positive numerics. The estimated
#'     standard errors of \code{betahat_ols}.
#' @param degrees_freedom A positive numeric. The degrees of freedom
#'     of the t-likelihood if using it.
#' @inheritParams vicarius_ruv4
#'
#' @author David Gerard
#'
cruv4 <- function(betahat_ols, alpha_scaled, sig_diag_scaled, ctl, degrees_freedom,
                  gls = TRUE, likelihood = "t") {

    assertthat::assert_that(is.matrix(betahat_ols))
    assertthat::assert_that(is.matrix(alpha_scaled))
    assertthat::assert_that(is.vector(sig_diag_scaled))
    assertthat::assert_that(all(sig_diag_scaled > 0))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(length(degrees_freedom) == 1 | length(degrees_freedom) == ncol(betahat_ols))
    assertthat::assert_that(all(degrees_freedom > 0))
    assertthat::are_equal(nrow(betahat_ols), nrow(alpha_scaled))

    k <- ncol(alpha_scaled)

    ## Use control genes to jointly estimate Z1 and variance scaling parameter.
    Yc <- betahat_ols[ctl, , drop = FALSE]
    if (k != 0) {
        alphac       <- alpha_scaled[ctl, , drop = FALSE]
        sig_diag_inv <- 1 / sig_diag_scaled[ctl]
        if (gls) {
            Z1 <- crossprod(solve(crossprod(alphac, sig_diag_inv * alphac)),
                            crossprod(alphac, sig_diag_inv * Yc))
        } else {
            Z1 <- crossprod(solve(crossprod(alphac, alphac)),
                            crossprod(alphac, Yc))
        }
        resid_mat <- Yc - alphac %*% Z1
        betahat <- betahat_ols - alpha_scaled %*% Z1
    } else {
        Z1        <- NULL
        resid_mat <- Yc
        betahat   <- betahat_ols
    }

    ## Gaussian MLE of variance inflation parameter.
    multiplier <- mean(resid_mat ^ 2 / sig_diag_scaled[ctl])

    ## T-likelihood regression and estimate variance inflation
    ## parameter using control genes.
    if (likelihood == "t") {
        if (k != 0) {
            if (!gls) {
                message("gls = FALSE not supported for t-likelihood, using gls = TRUE.")
            }
            alphac <- alpha_scaled[ctl, , drop = FALSE]
            tout <- tregress_em(Y = Yc, alpha = alphac,
                                sig_diag = sig_diag_scaled[ctl],
                                nu = degrees_freedom, lambda_init = multiplier,
                                Z_init = Z1)
            Z1         <- tout$Z
            betahat    <- betahat_ols - alpha_scaled %*% Z1
            multiplier <- tout$lambda
        } else if (k == 0) {
            tout <- stats::optim(par = multiplier, fn = tregress_obj_wrapper,
                                 Z = matrix(0, nrow = 1, ncol = 1), Y = Yc,
                                 alpha = matrix(0, nrow = nrow(Yc), ncol = 1),
                                 sig_diag = sig_diag_scaled[ctl],
                                 nu = degrees_freedom, method = "Brent",
                                 lower = 0, upper = 10,
                                 control = list(fnscale = -1))
            multiplier <- tout$par
        }
    }

    ## Output values
    sebetahat <- sqrt(sig_diag_scaled * multiplier)

    ruv4_out <- list()
    ruv4_out$multiplier    <- multiplier
    ruv4_out$betahat_ols   <- betahat_ols
    ruv4_out$sebetahat_ols <- sqrt(sig_diag_scaled)
    ruv4_out$betahat       <- betahat
    ruv4_out$sebetahat     <- sebetahat
    ruv4_out$tstats        <- betahat / sebetahat
    ruv4_out$pvalues       <- 2 * (stats::pt(q = -abs(ruv4_out$tstats),
                                             df = degrees_freedom))

    ruv4_out$Z1            <- Z1

    return(ruv4_out)
}


#' QR rotation to independent models.
#'
#' @inheritParams vicarius_ruv4
#'
#' @author David Gerard
#'
rotate_model <- function(Y, X, k, cov_of_interest = ncol(X),
                         include_intercept = TRUE,
                         limmashrink = FALSE, fa_func = pca_naive,
                         fa_args = list(), do_ols = FALSE) {

    if (include_intercept) {
        X_scaled <- apply(X, 2, function(x) {
            x / sqrt(sum(x ^ 2))
        })
        int_term <- rep(1, length = nrow(X)) / sqrt(nrow(X))

        any_int <- any(colSums( (int_term - X_scaled) ^ 2) < 10 ^ -14)
        if (!any_int) {
            X <- cbind(X, rep(1, length = nrow(X)))
        }
    }

    if (is.null(k)) {
        if (requireNamespace("sva", quietly = TRUE)) {
            message("Number of confounders not provided so being estimated with package sva.")
            k <- sva::num.sv(dat = t(Y), mod = X)
        } else {
            stop("If sva is not installed, then k needs to be provided. To install sva, run in R\n   source(\"https://bioconductor.org/biocLite.R\")\n   biocLite(\"sva\")")
        }
    }

    assertthat::assert_that(k + ncol(X) < nrow(X))

    if (k >= sum(ctl) & !do_ols) {
        stop("k is larger than the number of control genes so model not identified.\nReduce k or increase the number of control genes.\nYou can also try out succotashr. To install succotashr, run in R:\n    install.packages(\"devtools\")\n    devtools::install_github(\"dcgerard/succotashr\")")
    }

    ## Place desired covariate as last covariate
    X <- X[, c( (1:ncol(X))[-cov_of_interest], cov_of_interest), drop = FALSE]
    cov_of_interest <- (ncol(X) - length(cov_of_interest) + 1):ncol(X)

    qr_x <- qr_ident(X)
    ## multiply by sign so that it matches with beta_hat_ols
    Q <- qr_x$Q
    R <- qr_x$R
    ## discard first q-1 rows.
    Y_tilde <- crossprod(Q, Y)

    y2start_index <- ncol(X) - length(cov_of_interest) + 1
    y3start_index <- ncol(X) + 1
    Y2 <- Y_tilde[y2start_index:(y3start_index - 1), , drop = FALSE]

    Y3 <- Y_tilde[y3start_index:nrow(Y_tilde), , drop = FALSE]


    ## Factor analysis using all but first row of Y_tilde
    fa_args$Y <- Y3
    fa_args$r <- k
    fa_out    <- do.call(what = fa_func, args = fa_args)
    alpha     <- fa_out$alpha
    sig_diag  <- fa_out$sig_diag

    ## make sure the user didn't screw up the factor analysis.
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::are_equal(length(sig_diag), ncol(Y))
    assertthat::assert_that(all(sig_diag > 0))
    if (k != 0) {
        assertthat::assert_that(is.matrix(alpha))
        assertthat::are_equal(ncol(alpha), k)
        assertthat::are_equal(nrow(alpha), ncol(Y))
    } else {
        assertthat::assert_that(is.null(alpha))
    }

    ## Shrink variances if desired.
    if (requireNamespace("limma", quietly = TRUE) & limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag,
                                       df = nrow(X) - ncol(X) - k)
        sig_diag <- limma_out$var.post
    } else if (!requireNamespace("limma", quietly = TRUE) & limmashrink) {
        stop("limmashrink = TRUE but limma not installed. To install limma, run in R:\n    source(\"https://bioconductor.org/biocLite.R\")\n    biocLite(\"limma\")")
    }

    ## absorb fnorm(X) into Y_tilde[1,], alpha, and sig_diag -----------------
    ## since dealt with sign earlier
    Rsub <- R[cov_of_interest, cov_of_interest, drop = FALSE]
    ## this is betahat from OLS, called Y1 in CATE.
    betahat_ols <- solve(Rsub) %*% Y2

    ## create list for returns
    return_list             <- list()
    return_list$betahat_ols <- betahat_ols
    return_list$alpha       <- alpha
    return_list$sig_diag    <- sig_diag
    return_list$Rsub        <- Rsub
    return_list$Y2          <- Y2
    return_list$Y3          <- Y3
    return_list$k           <- k
    return_list$X           <- X

    return(return_list)
}


#' An identified QR decomposition.
#'
#' This just re-jiggers the output from \code{\link[base]{qr}} so that
#' the upper triangular matrix has positive diagonal elements. This one will only work if the column dimension is less than or equal to the row dimension.
#'
#' @param X A matrix.
qr_ident <- function (X)
{
    assertthat::assert_that(ncol(X) <= nrow(X))
    assertthat::assert_that(is.matrix(X))
    qr_X <- qr(X)
    R <- qr.R(qr_X, complete = TRUE)
    Q <- qr.Q(qr_X, complete = TRUE)
    sign_vec <- sign(diag(R))
    Q[, 1:ncol(X)] <- t(t(Q)[1:ncol(X), , drop = FALSE] * sign_vec)
    R[1:ncol(X), ] <- sign_vec * R[1:ncol(X), , drop = FALSE]
    return(list(R = R, Q = Q))
}


#' Basic PCA.
#'
#' Glorified truncated SVD. The variance estimates are just the
#' column-wise mean-squares, except I divide by \code{nrow(Y) - r}
#' rather than by \code{nrow(Y)}. I allow \code{r = 0}.
#'
#'
#' @param Y A matrix of numerics. The data.
#' @param r the rank.
#'
#' @export
#'
#' @author David Gerard
pca_naive <- function (Y, r) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::are_equal(length(r), 1)
    assertthat::assert_that(r >= 0 & r < min(dim(Y)))

    if (r == 0) {
        Gamma <- NULL
        Z <- NULL
        Sigma <- apply(Y, 2, function(x) mean(x ^ 2))
    } else {
        svd_Y <- svd(Y)
        Gamma <- svd_Y$v[, 1:r, drop = FALSE] %*% diag(svd_Y$d[1:r], r, r) /
            sqrt(nrow(Y))
        Z <- sqrt(nrow(Y)) * svd_Y$u[, 1:r, drop = FALSE]
        Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) sum(x ^ 2)) / (nrow(Y) - r)
    }
    return(list(alpha = Gamma, Z = Z, sig_diag = Sigma))
}

#' EM algorithm to find regression coefficients using t-likelihood
#' when variances are known up to scale.
#'
#' When using a non-standard t-distribution for your regression
#' likelihood, there is a simple latent variable representation that
#' allows us to develop an EM algorithm. This is a very similar
#' procedure to Lange et al (1989).
#'
#'
#' @param Y A matrix of numerics with one column. The response
#'     variables.
#' @param alpha A matrix of numerics, the covariates. It must be that
#'     \code{nrow(Y)} is equal to \code{nrow(alpha)}.
#' @param sig_diag A vector of numerics. The variances of the elements
#'     in \code{Y}, but only assumed to be known up to a scaling
#'     factor.
#' @param nu A positive numeric scalar. The known degrees of freedom
#'     of the t-distribution.
#' @param lambda_init A positive numeric scalar. The initial value of
#'     the variance inflation parameter. Defaults to the estimate
#'     under Gaussian errors.
#' @param Z_init A matrix of numerics with one column. The initial
#'     value of the coefficients of \code{alpha}. Defaults to the
#'     estimate under Gaussian errors.
#' @param control_args A list of control arguments for the EM
#'     algorithm that is passed to SQUAREM.
#'
#' @author David Gerard
#'
#' @references Lange, Kenneth L., Roderick JA Little, and Jeremy MG
#'     Taylor. "Robust statistical modeling using the t distribution."
#'     Journal of the American Statistical Association 84.408 (1989):
#'     881-896.
#'
#' @seealso \code{\link{tregress_obj}} for the objective function that
#'     this function maximizes, \code{\link{tregress_fix}} for the
#'     fixed point iteration that this function calls.
tregress_em <- function(Y, alpha, sig_diag, nu, lambda_init = NULL,
                         Z_init = NULL, control_args = list()) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(alpha))
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::assert_that(all(sig_diag > 0))
    assertthat::assert_that(length(nu) == 1 | length(nu) == nrow(Y))
    assertthat::assert_that(all(nu > 0))
    assertthat::assert_that(is.list(control_args))
    assertthat::are_equal(nrow(Y), nrow(alpha))
    assertthat::are_equal(nrow(Y), length(sig_diag))

    sig_diag_inv <- 1 / sig_diag
    if (is.null(Z_init)) {
        Z_init <- crossprod(solve(crossprod(alpha, sig_diag_inv * alpha)),
                            crossprod(alpha, sig_diag_inv * Y))
    }

    resid_init <- Y - alpha %*% Z_init

    if (is.null(lambda_init)) {
        lambda_init <- mean(resid_init ^ 2 * sig_diag_inv)
    }

    assertthat::assert_that(is.matrix(Z_init))
    assertthat::are_equal(ncol(alpha), nrow(Z_init))
    assertthat::are_equal(length(lambda_init), 1)
    assertthat::assert_that(lambda_init > 0)

    zlambda <- c(c(Z_init), lambda_init)

    sqout <- SQUAREM::squarem(par = zlambda, fixptfn = tregress_fix,
                              objfn = tregress_obj, Y = Y, alpha = alpha,
                              sig_diag = sig_diag, nu = nu,
                              control = control_args)

    Z <- matrix(sqout$par[-length(sqout$par)], nrow = ncol(alpha))
    lambda <- sqout$par[length(sqout$par)]

    return(list(Z = Z, lambda = lambda))
}

#' Fixed point iteration for EM algorithm for regression with
#' t-errors where the variance is known up to scale.
#'
#' @inheritParams tregress_em
#' @param zlambda A vector containing the current estimates of the
#'     coefficients (Z) and the variance inflation parameter
#'     (lambda). The last element is lambda and all other elements are
#'     Z.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{tregress_em}} where this function is called,
#'     \code{\link{tregress_obj}} for the objective function that this
#'     fixed-point iteration increases.
#'
tregress_fix <- function(zlambda, Y, alpha, sig_diag, nu) {

    assertthat::assert_that(is.vector(zlambda))
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(alpha))
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::are_equal(ncol(alpha), length(zlambda) - 1)
    assertthat::are_equal(nrow(alpha), nrow(Y))
    assertthat::assert_that(length(nu) == 1 | length(nu) == nrow(Y))
    assertthat::assert_that(all(nu > 0))

    Z <- matrix(zlambda[-length(zlambda)], nrow = ncol(alpha))
    lambda <- zlambda[length(zlambda)]

    resid_vec <- Y - alpha %*% Z

    wvec <- (nu + 1) / (resid_vec ^ 2 / (lambda * sig_diag) + nu)

    wsig <- c(wvec / sig_diag)

    Znew <- crossprod(solve(crossprod(alpha, wsig * alpha)),
                            crossprod(alpha, wsig * Y))

    resid_new <- Y - alpha %*% Znew

    lambda_new <- mean(resid_new ^ 2 * wsig)

    zlambda_new <- c(c(Znew), lambda_new)

    return(zlambda_new)

}

#' The likelihood for regression with t-errors where the variance is
#' known up to scale.
#'
#' @inheritParams tregress_fix
#'
#' @author David Gerard
#'
#' @seealso \code{\link{tregress_em}} where this function is called,
#'     \code{\link{tregress_fix}} for the fixed point iteration that
#'     increases this objective function.
tregress_obj <- function(zlambda, Y, alpha, sig_diag, nu) {

    assertthat::assert_that(is.vector(zlambda))
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(alpha))
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::are_equal(ncol(alpha), length(zlambda) - 1)
    assertthat::are_equal(nrow(alpha), nrow(Y))
    assertthat::assert_that(length(nu) == 1 | length(nu) == nrow(Y))
    assertthat::assert_that(all(nu > 0))

    Z <- matrix(zlambda[-length(zlambda)], nrow = ncol(alpha))
    lambda <- zlambda[length(zlambda)]

    standard_t<- (Y - alpha %*% Z) / sqrt(lambda * sig_diag)

    llike <- sum(log(stats::dt(x = standard_t, df = nu))) -
        length(sig_diag) * log(lambda) / 2 - sum(log(sig_diag)) / 2

    return(llike)
}

#' Wrapper for \code{\link{tregress_obj}}.
#'
#' This is mostly so I can run \code{stats::optim} with just
#' \code{lambda}.
#'
#' @author David Gerard
#'
#' @inheritParams tregress_obj
#' @param lambda A positive numeric. The variance inflation parameter.
#' @param Z A matrix of numerics with one column. The coefficients.
tregress_obj_wrapper <- function(lambda, Z, Y, alpha, sig_diag, nu) {
    zlambda <- c(c(Z), lambda)
    llike <- tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha,
                          sig_diag = sig_diag, nu = nu)
    return(llike)
}
