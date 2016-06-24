#' Use the distribution of control genes' p-values to adjust all of
#' the p-values.
#'
#' This function will transform all p-values by the empirical cdf of
#' the p-values for the control genes. This will result in the
#' p-values of the control genes being uniformly distributed about 0
#' and 1.
#'
#' @param pvalues A vector of numbers between 0 and 1. The pvalues.
#' @param ctl A vector of logicals of the same length of
#'     \code{pvalues}. A \code{TRUE} indicates that a p-value is from
#'     a control gene and so should be uniformly distributed.
#'
#' @author David Gerard
#'
#' @export
clean_my_p <- function(pvalues, ctl) {
    assertthat::assert_that(is.numeric(pvalues))
    assertthat::assert_that(is.logical(ctl))
    assertthat::are_equal(length(pvalues), length(ctl))
    ctl_ecdf <- stats::ecdf(pvalues[ctl])
    cleaned_p <- ctl_ecdf(pvalues)
    return(cleaned_p)
}

#' Wrapper for a Kolmogorov-Smirnov test using the p-values from the
#' control genes.
#'
#' @inheritParams clean_my_p
#'
#' @author David Gerard
#'
#' @export
smell_my_p <- function(pvalues, ctl) {
    return(stats::ks.test(pvalues[ctl], "punif"))
}


#' Quantile normalize stats to their theoretical distributions.
#'
#' @param stats A vector of numerics. The statistics.
#' @param ctl A vector of logicals of the same length of
#'     \code{pvalues}. A \code{TRUE} indicates that a p-value is from
#'     a control gene and so should be uniformly distributed.
#' @param qfunc A quantile function for the theoretical null
#'     distribution of the statistics in \code{stats}.
#' @param ... More parameters to pass to \code{qfunc}.
#' @author David Gerard
#'
#' @export
#'
#'
clean_my_stats <- function(stats, ctl, qfunc, ...) {
    assertthat::assert_that(is.numeric(stats))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(is.function(qfunc))
    assertthat::are_equal(length(stats), length(ctl))
    ctl_ecdf <- stats::ecdf(stats[ctl])
    cleaned_p <- ctl_ecdf(stats)
    cleaned_stats <- qfunc(cleaned_p, ...)
    return(cleaned_stats)
}

#' Wrapper for a Kolmogorov-Smirnov test using the statistics from the
#' control genes.
#'
#' @inheritParams clean_my_stats
#' @param cdf_func A cdf function.
#'
#' @author David Gerard
#'
#' @export
smell_my_stats <- function(stats, ctl, cdf_func, ...) {
    ksout <- stats::ks.test(x = stats[ctl], y = cdf_func, ...)
    return(ksout)
}
