#' \code{vicar}: Various Ideas for Confounder Adjustment in Regression.
#'
#' \code{vicar} contains functions for implementing RUV3 and two
#' versions of what we call RUV5: RUV-Bayes and RUV-impute. It also
#' contains the functions \code{\link{mouthwash}} and \code{\link{backwash}} for confounder
#' adjustment when one does not have access to control genes.
#'
#' @section \code{vicar} functions:
#'
#'     \code{\link{ruv3}}: An implementation of RUV3.
#'
#'     \code{\link{ruvb}}: An implementation of RUV-Bayes.
#'
#'     \code{\link{ruvimpute}}: An implementation of RUV-impute.
#'
#'     \code{\link{vruv4}}: Variance inflated version of CATE.
#'
#'     \code{\link{mouthwash}}: Empirical Bayesian confounder adjustment without control genes.
#'
#'     \code{\link{backwash}}: Variational Bayesian confounder adjustment without control genes.
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
