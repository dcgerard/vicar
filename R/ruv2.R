#' Calibrated RUV2.
#'
#' This function will perform a variant of Removing Unwanted Variation
#' 2-step (RUV2) (Gagnon-Bartsch et al, 2013), where we include a
#' variance inflation parameter in the factor analysis.
#'
#' See \code{\link{vruv4}} for a description of the model.
#'
#' You can provide your own factor analysis, but it must include an
#' estimate for the variance inflation parameter. This turns out to be
#' pretty hard. The way I do it now seems to work OK.
#'
#'
#' @inheritParams vruv4
#' @param gls A logical. Should we estimate the part of the
#'     confounders associated with the nuisance parameters with gls
#'     (\code{TRUE}) or with ols (\code{FALSE}).
#' @param fa_func A factor analysis function. It must take parameters:
#'     \code{Y} a data matrix of numerics, \code{r} a positive integer
#'     for the rank, and \code{vr} a positive integer for the number
#'     of the first rows that have a different variance than the last
#'     rows. It must return: \code{alpha} a matrix for the factor
#'     loadings, \code{Z} a matrix for the factors, \code{sig_diag} a
#'     vector of the column-wise variances, and \code{lambda} a
#'     numeric for the variance inflation of the first \code{vr} rows
#'     of \code{Y}. The default function is \code{\link{pca_2step}},
#'     which is the main difference between RUV2 and this version.
#' @param use_factor A logical. Should we use the estimates of
#'     \code{alpha} and \code{sig_diag} from the factor analysis
#'     (\code{TRUE}), or re-estimate these using OLS as RUV2 does it
#'     (\code{FALSE})? Right now it's probably a bad idea to have the
#'     settings \code{use_factor = TRUE, fa_limmashrink = TRUE,
#'     limmashrink = TRUE} since then the variance estimates of the
#'     control genes are being shrunk twice.
#' @param force_check A logical. Are you REALLY sure you want to use
#'     another fa_func (\code{FALSE}) or should I ask you again
#'     (\code{TRUE})?
#' @param fa_limmashrink A logical. Should we shrink the variances
#'     during the factor analysis step (\code{TRUE}) or not
#'     (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @seealso \code{\link{pca_2step}} for the special factor analysis
#'     that results in variance inflation in RUV2.
#'
#' @return A list whose elements are:
#'
#'     \code{betahat} A matrix of numerics. The ordinary least
#'     squares estimates of the coefficients of the covariate of
#'     interest WHEN YOU ALSO INCLUDE THE ESTIMATES OF THE UNOBSERVED
#'     CONFOUNDERS.
#'
#'     \code{sebetahat} A matrix of positive numerics. This is the
#'     post-inflation adjusted standard errors for \code{ruv$betahat}.
#'
#'     \code{tstats} A vector of numerics. The t-statistics for
#'     testing against the null hypothesis of the coefficient of the
#'     covariate of interest being zero. This is after estimating the
#'     variance inflation parameter.
#'
#'     \code{pvalues} A vector of numerics. The p-values of said test
#'     above.
#'
#'     \code{alphahat} A matrix of numerics. The estimates of the
#'     coefficients of the hidden confounders. Only identified up to a
#'     rotation on the rowspace.
#'
#'     \code{Zhat} A matrix of numerics. The estimates of the
#'     confounding variables. Only identified up to a rotation on the
#'     columnspace.
#'
#'     \code{sigma2} A vector of positive numerics. The estimates of
#'     the variances PRIOR to inflation.
#'
#'     \code{sigma2_adjusted} A vector of positive numerics. The
#'     estimates of the variances AFTER to inflation. This is equal to
#'     \code{sigma2 * multiplier}.
#'
#'     \code{multiplier} A numeric. The estimated variance inflation
#'     parameter.
#'
#'     \code{mult_matrix} A matrix of numerics. Equal to
#'     \code{solve(t(cbind(X, Zhat)) \%*\% cbind(X, Zhat))}. One
#'     multiplies \code{sigma2} or \code{simga2_adjused} by the
#'     diagonal elements of \code{mult_matrix} to get the standard
#'     errors of \code{betahat}.
#'
#'     \code{degrees_freedom} The degrees of freedom of the t-
#'     statistics.
#'
#'
#' @references
#' \itemize{
#' \item{Gagnon-Bartsch, J., Laurent Jacob, and Terence
#'     P. Speed. "Removing unwanted variation from high dimensional
#'     data with negative controls."
#'     Berkeley: Department of Statistics. University of California
#'     (2013).}
#' \item{Andreas Buja and Nermin
#'     Eyuboglu. "Remarks on parallel analysis." Multivariate behavior
#'     research, 27(4):509-540, 1992.}
#' \item{Gerard, David, and Matthew Stephens. 2019. "Unifying and Generalizing Methods for Removing Unwanted Variation Based on Negative Controls." \emph{Statistica Sinica}, in press. <\href{https://doi.org/10.5705/ss.202018.0345}{doi:10.5705/ss.202018.0345}>.}
#' }
vruv2 <- function(Y, X, ctl, k = NULL, cov_of_interest = ncol(X),
                  likelihood = c("t", "normal"), limmashrink = TRUE,
                  degrees_freedom = NULL, include_intercept = TRUE,
                  gls = TRUE, fa_func = pca_2step, fa_args = list(),
                  use_factor = FALSE, force_check = TRUE,
                  fa_limmashrink = TRUE) {

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

    if (!identical(fa_func, pca_2step) & !force_check) {
        cat ("vruv2 mostly works because of the special factor analysis used.\nAre you sure you want to use another factor analysis function (Y/N)?\n")
        line <- readline()
        if (line == "N") {
            return()
        } else if (line != "Y") {
            cat("Input Y or N\n")
        }
        while (line != "N" & line != "Y") {
            line <- readline()
            if (line == "N") {
                return()
            } else if (line != "Y") {
                cat("Input Y or N\n")
            }
        }
    }

    if (fa_limmashrink & limmashrink & use_factor) {
        warning("Variances of control genes are being shrunk twice.\nThis is probably bad.\nIt is recommended that you change one of limmashrink, fa_limmashrink, or use_factor.")
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
    fa_args$vr <- nrow(Y2)
    fa_args$likelihood <- likelihood
    fa_args$limmashrink <- fa_limmashrink
    pcout <- do.call(what = fa_func, args = fa_args)
    sig_diag_ctl <- pcout$sig_diag
    Z23 <- pcout$Z
    alpha_ctl <- pcout$alpha
    multiplier <- pcout$lambda

    ## Do RUV2 here -----------------------------------------------------------
    Z2 <- Z23[1:length(cov_of_interest), , drop = FALSE]
    Z3 <- Z23[(length(cov_of_interest) + 1):nrow(Z23), , drop = FALSE]
    alpha <- matrix(NA, ncol = ncol(Y), nrow = ncol(Z3))

    ## RUV2 method -----------------------------------------------------------
    alphahat_ruv2 <- solve(t(Z3) %*% Z3) %*% t(Z3) %*% Y3
    sigma2_ruv2 <- colSums((Y3 - Z3 %*% alphahat_ruv2) ^ 2) / (nrow(X) - ncol(X) - k)
    betahat_ruv2 <- solve(R22) %*% (Y2 - Z2 %*% alphahat_ruv2)

    if (use_factor) {
        alpha[, !ctl] <- alphahat_ruv2[, !ctl]
        alpha[, ctl]  <- t(alpha_ctl)
        betahat <- solve(R22) %*% (Y2 - Z2 %*% alpha)

        ## variances
        r2 <- colSums((Y3 - Z3 %*% alpha) ^ 2) / (nrow(X) - ncol(X) - k)
        sig_diag <- rep(NA, length = ncol(Y))
        sig_diag[ctl] <- sig_diag_ctl
        sig_diag[!ctl] <- r2[!ctl]
        sig_diag <- sig_diag
    } else {
        alpha <- alphahat_ruv2
        betahat <- betahat_ruv2
        sig_diag <- sigma2_ruv2
    }

    ## Shrink standard errors if desired
    if (limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag,
                                       df = nrow(X) - ncol(X) - k)
        sig_diag <- limma_out$var.post
        prior_df <- limma_out$df.prior
        degrees_freedom <- prior_df + nrow(X) - ncol(X) - k
    } else {
        degrees_freedom <- nrow(X) - ncol(X) - k
    }
    df_ruv2 <- nrow(X) - ncol(X) - k

    sig_diag_adjusted <- sig_diag * multiplier



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


    XZ <- cbind(X, Zhat)
    mult_matrix <- solve(t(XZ) %*% XZ)[cov_of_interest, cov_of_interest, drop = FALSE]

    sebetahat <- sqrt(outer(diag(mult_matrix), sig_diag_adjusted, FUN = "*"))
    sebetahat_ruv2 <- sqrt(outer(diag(mult_matrix), sigma2_ruv2, FUN = "*"))

    ruv2_out                   <- list()
    ruv2_out$betahat           <- betahat
    ruv2_out$sebetahat         <- sebetahat
    ## ruv2_out$betahat_ruv2      <- betahat_ruv2
    ## ruv2_out$sebetahat_ruv2    <- sebetahat_ruv2
    ## ruv2_out$sigma2_ruv2       <- sigma2_ruv2
    ruv2_out$sigma2            <- sig_diag
    ruv2_out$sigma2_adjusted   <- sig_diag_adjusted
    ruv2_out$tstats            <- betahat / sebetahat
    ruv2_out$pvalues           <- 2 * (stats::pt(q = -abs(ruv2_out$tstats),
                                                 df = degrees_freedom))
    ruv2_out$degrees_freedom   <- degrees_freedom
    ruv2_out$alphahat          <- alpha
    ruv2_out$Zhat              <- Zhat
    ruv2_out$multiplier        <- multiplier
    ruv2_out$mult_matrix       <- mult_matrix
    ruv2_out$debuglist         <- list()
    ruv2_out$debuglist$Z1      <- Z1
    ruv2_out$debuglist$Z2      <- Z2
    ruv2_out$debuglist$Z3      <- Z3

    return(ruv2_out)
}














#' Old version of Calibrated RUV2 that doesn't really work.
#'
#' This function will perform a variant of Removing Unwanted Variation
#' 2-step (RUV2) (Gagnon-Bartsch et al, 2013), where we include a
#' variance inflation parameter in the second step and estimate it by
#' maximum likelihood.
#'
#' See \code{\link{vruv4}} for a description of the model.
#'
#' The variances of the non-control genes using the adjustment for
#' RUV2 are unchanged. The variances for the control genes are the
#' ones that are inflated. Hence, for this method to adjust the
#' variances of the genes of interest (the non-control genes) you must
#' use hierarchical shrinkage of the variances. That is, you must set
#' \code{limmashrink = TRUE}.
#'
#'
#' @inheritParams vruv4
#' @param gls A logical. Should we estimate the part of the
#'     confounders associated with the nuisance parameters with gls
#'     (\code{TRUE}) or with ols (\code{FALSE}).
#'
#' @author David Gerard
#'
#'
#' @references
#' \itemize{
#' \item{Gagnon-Bartsch, J., Laurent Jacob, and Terence
#'     P. Speed. "Removing unwanted variation from high dimensional
#'     data with negative controls."
#'     Berkeley: Department of Statistics. University of California
#'     (2013).}
#' \item{Andreas Buja and Nermin
#'     Eyuboglu. "Remarks on parallel analysis." Multivariate behavior
#'     research, 27(4):509-540, 1992.}
#' \item{Gerard, David, and Matthew Stephens. 2019. "Unifying and Generalizing Methods for Removing Unwanted Variation Based on Negative Controls." \emph{Statistica Sinica}, in press. <\href{https://doi.org/10.5705/ss.202018.0345}{doi:10.5705/ss.202018.0345}>.}
#' }
vruv2_old <- function(Y, X, ctl, k = NULL,
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
