#' \code{vicar}: Various Ideas for Confounder Adjustment in Regression.
#'
#' \code{vicar} contains functions for implementing RUV3, a
#' Bayesian version of RUV, and a calibrated version of RUV4. It also
#' contains the functions \code{\link{mouthwash}} and \code{\link{backwash}} for confounder
#' adjustment when one does not have access to control genes.
#'
#' @section Main \code{vicar} functions:
#'
#'     \code{\link{ruv3}}: An implementation of RUV3.
#'
#'     \code{\link{ruvb}}: An implementation of RUV-Bayes.
#'
#'     \code{\link{vruv4}}: Variance inflated version of RUV4/CATE.
#'
#'     \code{\link{mouthwash}}: Empirical Bayesian confounder adjustment without control genes.
#'
#'     \code{\link{backwash}}: Variational Bayesian confounder adjustment without control genes.
#'
#' @section Vignettes (run in R):
#'
#'     \code{utils::vignette("customFA", package = "vicar")} For providing custom factor analyses during
#'     confounder adjustment.
#'
#'     \code{utils::vignette("custom_prior", package = "vicar")} For providing custom priors in \code{\link{ruvb}}.
#'
#'     \code{utils::vignette("sample_analysis", package = "vicar")} For a sample analysis on a simulated
#'     dataset using the vicar functions and other confounder adjustment packages.
#'
#'     If the vignette code does not work, then you probably did not build them during install. See
#'     \url{https://github.com/dcgerard/vicar#vignettes}.
#'
#' @docType package
#' @name vicar
#'
#' @author David Gerard
#'
#' @useDynLib vicar
#' @importFrom Rcpp evalCpp
NULL
#> NULL
