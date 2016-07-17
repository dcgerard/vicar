#' Use control genes to estimate hidden confounders and variance
#' inflation parameter, then run ASH.
#'
#' This function will perform a variant of Removing Unwanted Variation
#' 4-step (RUV4) (Gagnon-Bartsch et al, 2013) where the control genes
#' are used not only to estimate the hidden confounders, but to
#' estimate a variance inflation parameter. This variance inflation
#' step is akin to the "empirical null" approach of Efron
#' (2004). After this procedure, Adaptive SHrinkage (ASH) (Stephens,
#' 2016) is performed on the coefficient estimates and the inflated
#' standard errors.
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
#' modified to include estimating a variance inflation
#' parameter. Rather than use OLS in the second step of this two-step
#' procedure, we estimate the coefficients using Adaptive SHrinkage
#' (ASH) (Stephens, 2016). In the current implementation, only the
#' coefficients of one covariate can be estimated using ASH. The rest
#' are regressed out using OLS.
#'
#' @inheritParams vruv4
#' @param ash_args A list of arguments to pass to ash. See
#'     \code{\link[ashr]{ash.workhorse}} for details.
#' @param cov_of_interest A positive integer. The column number of the
#'     covariate in X whose coefficients you are interested in. The
#'     rest are considered nuiszance parameters and are regressed out
#'     by OLS. \code{ash_ruv4} only works with one covariate of
#'     interest right now.
#'
#' @author David Gerard
#'
#' @export
#'
#' @return Except for the list \code{ruv4}, the values returned are
#'     the exact same as in \code{\link[ashr]{ash.workhorse}}. See that
#'     function for more details. Elements in the \code{ruv4} are the
#'     exact same as returned in \code{\link{vruv4}}.
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
#'     Stephens, Matthew. "False Discovery Rates: A New Deal." bioRxiv
#'     (2016): 038216.
#'
#'     Wang, J., Zhao, Q., Hastie, T., & Owen, A. B
#'     "Confounder Adjustment in Multiple Hypotheses Testing."
#'     arXiv preprint arXiv:1508.04178 (2015).
ash_ruv4 <- function(Y, X, ctl = NULL, k = NULL, cov_of_interest = ncol(X),
                     likelihood = c("t", "normal"), ash_args = list(),
                     limmashrink = TRUE, degrees_freedom = NULL,
                     include_intercept = TRUE, gls = TRUE,
                     fa_func = pca_naive, fa_args = list()) {

    if (!requireNamespace("ashr", quietly = TRUE)) {
        stop("ashr must be installed to run ash_ruv4. To install, run in R:\n    install.packages(\"devtools\")\n    devtools::install_github(\"stephens999/ashr\")")
    }

    likelihood <- match.arg(likelihood)

    if (length(cov_of_interest) != 1) {
        stop("cov_of_interest must contain only one number")
    }

    vout <- vruv4(Y = Y, X = X, ctl = ctl, k = k,
                  cov_of_interest = cov_of_interest,
                  likelihood = likelihood, limmashrink = limmashrink,
                  degrees_freedom = degrees_freedom,
                  include_intercept = include_intercept, gls = gls,
                  fa_func = fa_func, fa_args = fa_args)

    betahat   <- c(vout$betahat)
    sebetahat <- c(vout$sebetahat)
    df        <- vout$degrees_freedom

    ash_args$betahat   <- betahat
    ash_args$sebetahat <- sebetahat
    if (likelihood == "t") {
        ash_args$df <- df
    }
    ashout <- do.call(what = ashr::ash.workhorse, args = ash_args)

    ashout$ruv4 <- vout

    return(ashout)
}





#' Use control genes to estimate hidden confounders and variance
#' inflation parameter, then run ASH.
#'
#' This model is fit using a two-step approach proposed in
#' Gagnon-Bartsch et al (2013) modified to include estimating a
#' variance inflation parameter. Rather than use OLS in the second
#' step of this two-step procedure, we estimate the coefficients using
#' Adaptive SHrinkage (ASH) (Stephens, 2016). In the current
#' implementation, only the coefficients of one covariate can be
#' estimated using ASH. The rest are regressed out using OLS.
#'
#' @inheritParams vruv4
#' @inheritParams vruv2
#' @param ash_args A list of arguments to pass to ash. See
#'     \code{\link[ashr]{ash.workhorse}} for details.
#' @param cov_of_interest A positive integer. The column number of the
#'     covariate in X whose coefficients you are interested in. The
#'     rest are considered nuiszance parameters and are regressed out
#'     by OLS. \code{ash_ruv2} only works with one covariate of
#'     interest right now.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{vruv2}} for the variance inflation in RUV2 and
#'     \code{\link[ashr]{ash.workhorse}} for the adaptive shrinkage
#'     method.
#'
#' @return Except for the list \code{ruv2}, the values returned are
#'     the exact same as in \code{\link[ashr]{ash.workhorse}}. See that
#'     function for more details. Elements in the \code{ruv2} are the
#'     exact same as returned in \code{\link{vruv2}}.
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
ash_ruv2 <- function(Y, X, ctl, k = NULL, cov_of_interest = ncol(X),
                     likelihood = c("t", "normal"), ash_args = list(),
                     limmashrink = TRUE, degrees_freedom = NULL,
                     include_intercept = TRUE, gls = TRUE,
                     fa_func = pca_2step, fa_args = list(),
                     use_factor = FALSE) {

    if (!requireNamespace("ashr", quietly = TRUE)) {
        stop("ashr must be installed to run ash_ruv2. To install, run in R:\n    install.packages(\"devtools\")\n    devtools::install_github(\"stephens999/ashr\")")
    }
    
    likelihood <- match.arg(likelihood)
    
    if (length(cov_of_interest) != 1) {
        stop("cov_of_interest must contain only one number")
    }
    
    vout <- vruv2(Y = Y, X = X, ctl = ctl, k = k,
                  cov_of_interest = cov_of_interest,
                  likelihood = likelihood, limmashrink = limmashrink,
                  degrees_freedom = degrees_freedom,
                  include_intercept = include_intercept, gls = gls,
                  fa_func = fa_func, fa_args = fa_args)

    betahat   <- c(vout$betahat)
    sebetahat <- c(vout$sebetahat)
    df        <- vout$degrees_freedom
    
    ash_args$betahat   <- betahat
    ash_args$sebetahat <- sebetahat
    if (likelihood == "t") {
        ash_args$df <- df
    }
    ashout <- do.call(what = ashr::ash.workhorse, args = ash_args)

    ashout$ruv2 <- vout

    return(ashout)
}




#' Random draw from a mixture of normals.
#'
#' @param n the number of samples to draw.
#' @param mixdense An object of class \code{normalmix}.
#'
#' @author David Gerard
#'
#' @export
rnormalmix <- function(n, mixdense) {
    assertthat::are_equal(class(mixdense), "normalmix")
    k <- length(mixdense$pi)
    which_norm <- sample(1:k, size = n, prob = mixdense$pi, replace = TRUE)
    samp <- stats::rnorm(n = n, mean = mixdense$mean[which_norm], sd = mixdense$sd[which_norm])
    return(samp)
}


#' Draw from a mixture of normals.
#'
#' Draw from a mean zero mixture of normals given mixing proportions,
#' mixing standard deviations, and mixing means.
#'
#' Given mixing proportions \code{pi_vals}, and mixing standard
#' deviations (not variances) \code{sd_seq}, this function will draw
#' \code{p} samples from the mean zero mixture of normals.
#'
#' @param pi_vals A vector length \code{M}. The mixing proportions.
#' @param sd_seq A vector of length \code{M}. The mixing standard
#'     deviations.
#' @param mean_seq A vector of length \code{M}. The mixing means.
#' @param n An integer. The number of samples to draw.
#'
#' @export
#'
#' @return A vector of length \code{p} that are drawn from the mean
#'     zero mixture of normals.
rmixnorm <- function(n, pi_vals, sd_seq, mean_seq = rep(0, length(pi_vals))) {
    assertthat::are_equal(length(pi_vals), length(sd_seq))
    assertthat::are_equal(length(pi_vals), length(mean_seq))
    assertthat::assert_that(is.numeric(n))
    assertthat::assert_that(is.numeric(pi_vals))
    assertthat::assert_that(is.numeric(sd_seq))
    assertthat::assert_that(is.numeric(mean_seq))
    assertthat::are_equal(length(n), 1)
    assertthat::assert_that(abs(n - round(n)) < 10 ^ -14)
    assertthat::assert_that(n >= 1)
    assertthat::are_equal(sum(pi_vals), 1)

    M <- length(pi_vals)
    beta <- rep(NA, length = n)
    which.mix <- sample(1:M, size = n, replace = TRUE, prob = pi_vals)
    for (index in 1:M) {
        current.ind <- which.mix == index
        n_m <- sum(current.ind)
        if (n_m > 0) {
            beta[current.ind] <- stats::rnorm(n = n_m,
                                              mean = mean_seq[index],
                                              sd = sd_seq[index])
        }
    }
    return(beta)
}
