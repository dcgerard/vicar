#' \code{vicar}: (more than just) Variance Inflation for Confounder Adjustment in
#' Regression.
#'
#' \code{vicar} does more than just variance inflation. This package
#' contains function for implementing RUV3, and two versions of what
#' we call RUV5: RUV-Bayes and RUV-impute.
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
#' @docType package
#' @name vicar
#'
#' @author David Gerard
#'
#' @useDynLib vicar
#' @importFrom Rcpp evalCpp
NULL
#> NULL
