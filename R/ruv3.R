#' Removing Unwanted Variation 3.
#'
#' This can be viewed both as a version of RUV2 and a version of RUV4.
#'
#' @inheritParams vruv4
#'
#' @return \code{betahat} The estimates of the coefficients of
#'     interest. The values corresponding to control genes are 0.
#'
#'     \code{sebetahat_unadjusted} The unadjusted standard errors of
#'     \code{betahat}. The values corresponding to control genes are
#'     \code{NA}.
#'
#'     \code{tstats_unadajusted} The t-statistics corresponding to the
#'     coefficients of interest. These use \code{sebetahat_unadjusted}
#'     as the standard errors. The values corresponding to control
#'     genes are \code{NA}.
#'
#'     \code{pvalues_unadajusted} The p-values using said statistics
#'     above.
#'
#'     \code{sebetahat_adjusted} The unadjusted standard errors of
#'     \code{betahat}. This equals \code{sebetahat_unadjusted *
#'     multiplier}. The values corresponding to control genes are
#'     \code{NA}.
#'
#'     \code{tstats_adjusted} The t-statistics corresponding to the
#'     coefficients of interest. These use \code{sebetahat_adjusted}
#'     as the standard errors. The values corresponding to control
#'     genes are \code{NA}.
#'
#'     \code{pvalues_unadajusted} The p-values using said statistics
#'     above.
#'
#'     \code{betahat_ols} The ordinary least squares (OLS) estimates
#'     for all of the coefficients.
#'
#'     \code{sebetahat_ols} The standard errors from OLS regression.
#'
#'     \code{tstats_ols} The t-statistics from OLS regression.
#'
#'     \code{pvalues_ols} The p-values from OLS regression.
#'
#'     \code{sigma2_unadjusted} The unadjusted variance estimates.
#'
#'     \code{sigma2_adjusted} The adjusted variance estimates. This is
#'     equal to \code{sigma2_unadjusted * multiplier}.
#'
#'     \code{Zhat} The estimates of the confounders.
#'
#'     \code{alphahat} The esitmates of the coefficients of the confounders.
#'
#'     \code{multiplier} The estimate of the variance inflation parameter.
#'
#'     \code{mult_matrix} The unscaled covariance of \code{betahat}
#'     after including the confounders.
#'
#'     \code{mult_matrix_ols} The OLS version of \code{mult_matrix}.
#'
#'     \code{degrees_freedom} The degrees of freedom used when
#'     calculating the p-values.
#'
#'     \code{debuglist} A list of elements that aren't really useful
#'     except for unit testing and debugging.
#'
#' @export
#'
#' @author David Gerard
ruv3 <- function(Y, X, ctl, k = NULL, cov_of_interest = ncol(X),
                 include_intercept = TRUE, limmashrink = TRUE,
                 gls = TRUE, fa_func = pca_naive, fa_args = list()) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(all(abs(cov_of_interest - round(cov_of_interest)) < 10 ^ -14))
    assertthat::assert_that(all(cov_of_interest >= 1 & cov_of_interest <= ncol(X)))
    assertthat::assert_that(is.logical(gls))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y))
    assertthat::assert_that(is.null(fa_args$r))
    assertthat::assert_that(is.function(fa_func))

    ## Rotate model ----------------------------------------------------------
    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink,
                               fa_func = fa_func, fa_args = fa_args,
                               do_factor = FALSE)

    degrees_freedom <- nrow(X) - ncol(X) - k

    Y21 <- rotate_out$Y2[, ctl, drop = FALSE]
    Y22 <- rotate_out$Y2[, !ctl, drop = FALSE]
    Y31 <- rotate_out$Y3[, ctl, drop = FALSE]
    Y32 <- rotate_out$Y3[, !ctl, drop = FALSE]
    R22 <- rotate_out$R22

    ## Factor analysis on Y31 ------------------------------------------------
    fa_args$Y <- Y31
    fa_args$r <- k
    faout <- do.call(what = fa_func, args = fa_args)

    alpha1 <- t(faout$alpha)
    Z3 <- faout$Z
    sig_diag1 <- faout$sig_diag

    ## Regression to get Z2 --------------------------------------------------
    if (gls) {
        if (limmashrink) {
            limma_out1 <- limma::squeezeVar(var = sig_diag1,
                                            df = degrees_freedom)
            sig_diag1_temp <- limma_out1$var.post
        } else {
            sig_diag1_temp <- sig_diag1
        }
        Z2 <- Y21 %*% diag(1 / sig_diag1_temp) %*% t(alpha1) %*%
            solve(alpha1 %*% diag(1 / sig_diag1_temp) %*% t(alpha1))
    } else {
        Z2 <- Y21 %*% t(alpha1) %*%
            solve(alpha1 %*% t(alpha1))
    }

    ## find variance adjustment ----------------------------------------------
    resid_mat <- Y21 - Z2 %*% alpha1
    multiplier <- mean(resid_mat ^ 2 %*% (1 / sig_diag1))

    ## Regression to get alpha2 ----------------------------------------------
    alpha2 <- solve(t(Z3) %*% Z3) %*% t(Z3) %*% Y32
    sig_diag2 <- colSums((Y32 - Z3 %*% alpha2) ^ 2) / degrees_freedom

    ## Get beta2 hat, consolodate estimates into big matrices ----------------
    beta2hat <- solve(R22) %*% (Y22 - Z2 %*% alpha2)
    beta2_ols <- rotate_out$betahat_ols

    betahat_long         <- matrix(0, nrow = nrow(beta2hat), ncol = ncol(Y))
    betahat_long[, !ctl] <- beta2hat
    alpha_long           <- matrix(NA, nrow = nrow(alpha2), ncol = ncol(Y))
    alpha_long[, ctl]    <- alpha1
    alpha_long[, !ctl]   <- alpha2
    sig_diag_long        <- rep(NA, length = ncol(Y))
    sig_diag_long[ctl]   <- sig_diag1
    sig_diag_long[!ctl]  <- sig_diag2

    ## Shrink variances if desired.
    if (limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag_long,
                                       df = degrees_freedom)
        sig_diag_long <- limma_out$var.post
        prior_df <- limma_out$df.prior
        degrees_freedom <- prior_df + degrees_freedom
    }

    ## Get remaining confounders ---------------------------------------------
    if (!is.null(rotate_out$Y1)) {
        R12 <- rotate_out$R12
        R11 <- rotate_out$R11
        Q   <- rotate_out$Q
        beta1_ols <- solve(R11) %*% (rotate_out$Y1 - R12 %*% beta2_ols)
        betahat_ols <- rbind(beta1_ols, beta2_ols)
        resid_top <- rotate_out$Y1 - R12 %*% betahat_long - R11 %*% beta1_ols
        if (gls) {
            Z1  <- resid_top %*% diag(1 / sig_diag_long) %*% t(alpha_long) %*%
                solve(alpha_long %*% diag(1 / sig_diag_long) %*% t(alpha_long))
        } else {
            Z1  <- resid_top %*% t(alpha_long) %*% solve(alpha_long %*% t(alpha_long))
        }
        Zhat <- Q %*% rbind(Z1, Z2, Z3)
    } else {
        Q   <- rotate_out$Q
        Zhat <- Q %*% rbind(Z2, Z3)
        betahat_ols <- beta2_ols
    }

    ## Calculate standard errors ---------------------------------------------
    XZ <- cbind(X, Zhat)
    mult_matrix <- solve(t(XZ) %*% XZ)[cov_of_interest, cov_of_interest, drop = FALSE]
    sebetahat_unadjusted <- sqrt(outer(diag(mult_matrix), sig_diag_long, FUN = "*"))
    sebetahat_unadjusted[, ctl] <- NA
    sebetahat_adjusted <- sebetahat_unadjusted * multiplier

    ## OLS statistics --------------------------------------------------------
    sigma2_ols <- colSums((Y - rotate_out$X %*% betahat_ols) ^ 2) / (nrow(X) - ncol(X))
    mult_matrix_ols <- solve(t(rotate_out$X) %*% rotate_out$X)
    sebetahat_ols <- sqrt(outer(diag(mult_matrix_ols), sigma2_ols, FUN = "*"))
    tstats_ols <- betahat_ols / sebetahat_ols
    pvalues_ols <- 2 * (stats::pt(q = -abs(tstats_ols), df = nrow(X) - ncol(X)))

    ## Make return list ------------------------------------------------------
    return_list <- list()
    return_list$betahat              <- betahat_long
    return_list$sebetahat_unadjusted <- sebetahat_unadjusted
    return_list$tstats_unadjusted    <- betahat_long / sebetahat_unadjusted
    return_list$pvalues_unadjusted   <- 2 * (stats::pt(q = -abs(return_list$tstats_unadjusted),
                                                     df = degrees_freedom))
    return_list$sebetahat_adjusted   <- sebetahat_adjusted
    return_list$tstats_adjusted      <- betahat_long / sebetahat_adjusted
    return_list$pvalues_adjusted     <- 2 * (stats::pt(q = -abs(return_list$tstats_adjusted),
                                                     df = degrees_freedom))
    return_list$betahat_ols          <- betahat_ols
    return_list$sebetahat_ols        <- sebetahat_ols
    return_list$tstats_ols           <- tstats_ols
    return_list$pvalues_ols          <- pvalues_ols
    return_list$simga2_unadjusted    <- sig_diag_long
    return_list$sigma2_adjusted      <- sig_diag_long * multiplier
    return_list$Zhat                 <- Zhat
    return_list$alphahat             <- alpha_long
    return_list$multiplier           <- multiplier
    return_list$mult_matrix          <- mult_matrix
    return_list$mult_matrix_ols      <- mult_matrix_ols
    return_list$degrees_freedom      <- degrees_freedom
    return_list$debuglist            <- list()
    if (!is.null(rotate_out$Y1)) {
        return_list$debuglist$Z1         <- Z1
    }
    return_list$debuglist$Z2         <- Z2
    return_list$debuglist$Z3         <- Z3

    return(return_list)
}
