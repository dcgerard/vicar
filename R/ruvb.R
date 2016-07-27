#' Bayesian version of Removing Unwanted Variation.
#'
#' @inheritParams vruv4
#' @param fa_func A function that takes as input a matrix names
#'     \code{Y} that has missing values and returns a list of matrices
#'     called \code{Yhat} where each matrix is the same dimension of
#'     \code{Y} with the missing values filled in.
#' @param fa_args A list of additional parameters to pass to
#'     \code{fa_func}.
#'
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
ruvb <- function(Y, X, ctl, k = NULL, fa_func, fa_args = list(),
                 cov_of_interest = ncol(X), include_intercept = TRUE) {

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
    assertthat::assert_that(is.function(fa_func))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y21))
    assertthat::assert_that(is.null(fa_args$Y31))
    assertthat::assert_that(is.null(fa_args$Y32))

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

    ncontrols <- sum(ctl)
    ncovariates <- length(cov_of_interest)

    Ytilde <- matrix(NA, nrow = nrow(Y21) + nrow(Y31), ncol = ncol(Y21) + ncol(Y22))
    Ytilde[1:ncovariates, 1:ncontrols] <- Y21
    Ytilde[(ncovariates + 1):nrow(Ytilde), 1:ncontrols] <- Y31
    Ytilde[(ncovariates + 1):nrow(Ytilde), (ncontrols + 1):ncol(Ytilde)] <- Y32

    R22inv <- backsolve(R22, diag(nrow(R22)))



    return()
}


rstiefel_wrapper <- function(Y, num_sv) {

}
