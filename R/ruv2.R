
#' Calibrated RUV2.
#'
#' This function will perform a variant of Removing Unwanted Variation
#' 2-step (RUV2) (Gagnon-Bartsch et al, 2013), where we include a
#' variance inflation parameter in the second step and estimate it by
#' maximum likelihood.
#'
#' See \code{\link{vruv4}} for a description of the model.
#'
#' The variances of the non-control genes using the adjustment for RUV2 are unchanged. The variances for the control genes are the ones that are inflated. Hence, for this method to adjust the variances of the genes of interest (the non-control genes) you must use hierarchical shrinkage of the variances. That is, you must set \code{limmashrink = TRUE}.
#'
#'
#' @inheritParams vruv4
#' @param gls A logical. Should we estimate the part of the
#'     confounders associated with the nuisance parameters with gls
#'     (\code{TRUE}) or with ols (\code{FALSE}).
#'
#' @author David Gerard
#'
#' @export
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
vruv2 <- function(Y, X, ctl, k = NULL,
                  cov_of_interest = ncol(X),
                  likelihood = c("t", "normal"),
                  limmashrink = TRUE, degrees_freedom = NULL,
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

    if (!limmashrink) {
        warning("If limmashrink = FALSE, then variance inflation in RUV2 will have no effect on the non-control genes.\nIt is recommended that you use limmashrink = TRUE .")
    }

    if (likelihood == "t") {
        stop("t-likelihood not implemented yet for vruv2")
    }


    ## RUN THE ROTATED MODEL HERE -------------------------------------------
    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink, fa_func = fa_func,
                               fa_args = fa_args, do_factor = FALSE)


    Y1  <- rotate_out$Y1
    Y2  <- rotate_out$Y2
    Y3  <- rotate_out$Y3
    R22 <- rotate_out$R22
    k   <- rotate_out$k

    if (k == 0) {
        stop("k = 0, so no point in running vruv2")
    }


    ## Do factor Analysis -----------------------------------------------------
    Y23 <- rbind(Y2, Y3)
    fa_args$Y <- Y23[, ctl]
    fa_args$r <- k
    pcout <- do.call(what = fa_func, args = fa_args)
    sig_diag_ctl <- pcout$sig_diag
    Z23 <- pcout$Z
    alpha_ctl <- pcout$alpha

    ## Do RUV2 here -----------------------------------------------------------
    Z2 <- Z23[1:length(cov_of_interest), , drop = FALSE]
    Z3 <- Z23[(length(cov_of_interest) + 1):nrow(Z23), , drop = FALSE]
    alpha <- matrix(NA, ncol = ncol(Y), nrow = ncol(Z3))

    ## RUV2 method
    alphahat_ols <- solve(t(Z3) %*% Z3) %*% t(Z3) %*% Y3
    sigma2_unadjusted <- colSums((Y3 - Z3 %*% alphahat_ols) ^ 2) / (nrow(X) - ncol(X) - k)
    betahat_ols <- solve(R22) %*% (Y2 - Z2 %*% alphahat_ols)

    ## New method
    alpha[, !ctl] <- alphahat_ols[, !ctl]
    alpha[, ctl]  <- t(alpha_ctl)
    betahat <- solve(R22) %*% (Y2 - Z2 %*% alpha)

    ## variances
    r2 <- colMeans((Y3 - Z3 %*% alpha) ^ 2)
    multiplier <- mean(r2[ctl] / sig_diag_ctl)
    sig_diag <- rep(NA, length = ncol(Y))
    sig_diag[ctl] <- sig_diag_ctl
    sig_diag[!ctl] <- r2[!ctl] / multiplier

    ## MLE to UMVUE inspired inflation
    multiplier <- multiplier * (nrow(X) - ncol(X)) / (nrow(X) - ncol(X) - k)

    ## Shrink variances
    if (limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag,
                                       df = nrow(X) - ncol(X) - k)
        sig_diag <- limma_out$var.post
        prior_df <- limma_out$df.prior

        limma_out_unadjusted <- limma::squeezeVar(var = sigma2_unadjusted,
                                                  df = nrow(X) - ncol(X) - k)
        sigma2_unadjusted <- limma_out_unadjusted$var.post
        prior_df_unadjusted <- limma_out_unadjusted$df.prior
    } else {
        prior_df <- nrow(X) - ncol(X) - k
        prior_df_unadjusted <- nrow(X) - ncol(X) - k
    }



    ## Estimate rest of hidden confounders
    if (!is.null(Y1)) {
        R12 <- rotate_out$R12
        R11 <- rotate_out$R11
        Q   <- rotate_out$Q
        beta1_ols <- solve(R11) %*% (Y1 - R12 %*% solve(R22) %*% Y2)
        resid_top <- Y1 - R12 %*% betahat - R11 %*% beta1_ols
        if (gls) {
            Z1  <- t(solve(alpha %*% diag(1 / sig_diag) %*% t(alpha)) %*%
                     alpha %*% diag(1 / sig_diag) %*% t(resid_top))
        } else {
            Z1  <- t(solve(alpha %*%  t(alpha)) %*% alpha %*% t(resid_top))
        }
        Zhat <- Q %*% rbind(Z1, Z2, Z3)
    } else {
        Zhat <- Q %*% rbind(Z2, Z3)
    }


    ## this is different mult_matrix from vruv4.
    ## This is the way the RUV folks do it and makes more sense.
    ## Need to think about consequences much more carefully.
    XZ <- cbind(X, Zhat)
    mult_matrix <- solve(t(XZ) %*% XZ)[cov_of_interest, cov_of_interest, drop = FALSE]

    sebetahat <- sqrt(outer(diag(mult_matrix), sig_diag * multiplier, FUN = "*"))
    sebetahat_ols <- sqrt(outer(diag(mult_matrix), sigma2_unadjusted, FUN = "*"))

    ruv2_out                   <- list()
    ruv2_out$betahat           <- betahat
    ruv2_out$sebetahat         <- sebetahat
    ruv2_out$betahat_ols       <- betahat_ols
    ruv2_out$sebetahat_ols     <- sebetahat_ols
    ruv2_out$sigma2_unadjusted <- sigma2_unadjusted
    ruv2_out$sigma2_adjusted   <- sig_diag * multiplier
    ruv2_out$tstats            <- betahat / sebetahat
    ruv2_out$pvalues           <- 2 * (stats::pt(q = -abs(ruv2_out$tstats),
                                                 df = prior_df))
    ruv2_out$alphahat          <- alpha
    ruv2_out$Zhat              <- Zhat
    ruv2_out$multiplier        <- multiplier
    ruv2_out$mult_matrix       <- mult_matrix

    return(ruv2_out)
}
