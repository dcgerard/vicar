
#' Tests whether a user-specified factor analysis function
#' is appropriate for use in \code{\link{vruv4}}, \code{\link{ruv3}},
#' \code{\link{mouthwash}}, or \code{\link{backwash}}.
#'
#' This function tests that the output of \code{fa_func} is compatible with the
#' confounder adjustment methods in this package. In general, \code{fa_func} needs to take as input
#' at least a matrix \code{Y}, which is dimensions n by p, and the rank of the factor analysis \code{r}.
#' \code{fa_func} needs to return three elements, a matrix \code{alpha} that is p by r, a matrix \code{Z} that
#' is n by r, and a numeric vector \code{sig_diag} that is of length \code{p}.
#'
#' This function also tests for orthogonal equivariance. See Gerard and Stephens (2016) for a definition
#' of left orthogonal equivariance.
#' You can turn off this test with \code{test_ortho = FALSE}, but any factor analysis that is not
#' left orthogonally equivariant would be dubious as the rotation onto the null space is completely
#' arbitrary.
#'
#' @param fa_func The factor analysis function to be tested.
#' @param fa_args A list of additional arguments to pass to \code{fa_func}.
#' @param test_ortho A logical. Should we test for left orthogonal equivariance (\code{TRUE}),
#'     or not (\code{FALSE})?
#'
#' @return A list with some or all of the following elements:
#'
#'     \code{ok} A logical. A value of \code{TRUE} means the function is OK. A value of
#'     \code{FALSE} means the function is not OK.
#'
#'     \code{why} A string vector giving the reason(s) why the function is not OK if \code{ok = FALSE}.
#'
#'     \code{input} A list of input values to help you debug.
#'
#' @export
#'
#' @author David Gerard
#'
fa_tester <- function(fa_func, fa_args = list(), test_ortho = TRUE) {

  assertthat::assert_that(is.list(fa_args))
  assertthat::assert_that(is.null(fa_args$Y))
  assertthat::assert_that(is.null(fa_args$r))

  if (!is.function(fa_func)) {
    return(list(ok = FALSE, why = "fa_func needs to be a function"))
  }


  ## Generate data -----------------------------------------------------------
  n <- 11
  p <- 23
  r <- 5
  left_mat <- matrix(stats::rnorm(n * r), nrow = n)
  right_mat <- matrix(stats::rnorm(r * p), nrow = r)
  Y <- left_mat %*% right_mat + matrix(stats::rnorm(n * p, sd = 1/2), nrow = n)

  ## test dimensions and classes ---------------------------------------------
  fa_args_1 <- fa_args
  fa_args_1$Y <- Y
  fa_args_1$r <- r
  tryout <- tryCatch( {
    fa_out <- do.call(what = fa_func, args = fa_args_1)
  }, error = function(e) {NULL})
  if (is.null(tryout)) {
    return(list(ok = FALSE, why = "fa_func gave error using input", input = fa_args_1))
  }
  faot <- fa_out_test(fa_out = fa_out, fa_args = fa_args_1)
  if (!faot$ok) {
    return(faot)
  }

  if (test_ortho) {
    ## test orthogonal equivariance --------------------------------------------
    Q <- qr.Q(qr(matrix(stats::rnorm(n ^ 2), nrow = n)))
    Y2 <- Q %*% Y

    fa_args_2 <- fa_args
    fa_args_2$Y <- Y2
    fa_args_2$r <- r
    fa_out_2 <- do.call(what = fa_func, args = fa_args_2)

    meanest1 <- tcrossprod(fa_out$Z, fa_out$alpha)
    meanest2 <- crossprod(Q, tcrossprod(fa_out_2$Z, fa_out_2$alpha))

    if (!all(abs(abs(meanest1) - abs(meanest2)) < .Machine$double.eps * 1000)) {
      return(list(ok = FALSE, why = "fa_func is not left orthogonally equivariant.",
                  input = list(fa_args_1, fa_args_2)))
    }

    if (!all(abs(fa_out$sig_diag - fa_out_2$sig_diag) < .Machine$double.eps * 1000)) {
      return(list(ok = FALSE, why = "variance estimates differ on left rotated data",
                  input = list(fa_args_1, fa_args_2)))
    }
  }


  ## test extreme data (r = 1) -----------------------------------------------
  n <- 11
  p <- 23
  r <- 1
  left_mat <- matrix(stats::rnorm(n * r), nrow = n)
  right_mat <- matrix(stats::rnorm(r * p), nrow = r)
  Y <- left_mat %*% right_mat + matrix(stats::rnorm(n * p, sd = 1/2), nrow = n)
  fa_args_1 <- fa_args
  fa_args_1$Y <- Y
  fa_args_1$r <- r
  tryout <- tryCatch( {
    fa_out <- do.call(what = fa_func, args = fa_args_1)
  }, error = function(e) {NULL})
  if (is.null(tryout)) {
    return(list(ok = FALSE, why = "fa_func gave error using input with r = 1", input = fa_args_1))
  }
  faot <- fa_out_test(fa_out = fa_out, fa_args = fa_args_1)
  if (!faot$ok) {
    return(faot)
  }

  return(list(ok = TRUE))
}

#' Test output form fa_func.
#'
#' @param fa_out Output from fa_func.
#' @param fa_args Input to fa_func.
#'
#' @author David Gerard
#'
fa_out_test <- function(fa_out, fa_args) {
  if (is.null(fa_out$alpha)) {
    return(list(ok = FALSE, why = "alpha is missing from returns", input = fa_args))
  }
  if (is.null(fa_out$Z)) {
    return(list(ok = FALSE, why = "Z is missing from returns", input = fa_args))
  }
  if (is.null(fa_out$sig_diag)) {
    return(list(ok = FALSE, why = "sig_diag is missing from returns", input = fa_args))
  }
  if (!is.matrix(fa_out$alpha)) {
    return(list(ok = FALSE, why = "alpha is not a matrix", input = fa_args))
  }
  if (!is.matrix(fa_out$Z)) {
    return(list(ok = FALSE, why = "Z is not a matrix", input = fa_args))
  }
  if (!is.vector(fa_out$sig_diag, mode = "numeric")) {
    return(list(ok = FALSE, why = "sig_diag is not a numeric vector", input = fa_args))
  }
  if (ncol(fa_out$alpha) != ncol(fa_out$Z)) {
    return(list(ok = FALSE, why = "ncol(alpha) does not equal ncol(Z)", input = fa_args))
  }
  if (nrow(fa_out$alpha) != ncol(fa_args$Y)) {
    return(list(ok = FALSE, why = "nrow(alpha) does not equal ncol(Y)", input = fa_args))
  }
  if (nrow(fa_out$Z) != nrow(fa_args$Y)) {
    return(list(ok = FALSE, why = "nrow(Z) does not equal nrow(Y)", input = fa_args))
  }
  if (length(fa_out$sig_diag) != ncol(fa_args$Y)) {
    return(list(ok = FALSE, why = "length(sig_diag) does not equal ncol(Y)", input = fa_args))
  }
  if (!all(fa_out$sig_diag > -1 * .Machine$double.eps * 1000)) {
    return(list(ok = FALSE, why = "sig_diag has some negative numbers", input = fa_args))
  }
  return(list(ok = TRUE))
}


#' Test whether a user-specified function is compatible with \code{\link{ruvb}}.
#'
#' This function tests whether a user-defined function \code{fa_func} is compatible
#' for use in \code{\link{ruvb}}. Let n be the number of samples, p be the number of genes,
#' m be the number of control genes, and q be the number of covariates.
#' In general, \code{fa_func} needs to take as input
#' \code{Y21} a q by m matrix of numerics, \code{Y31} a (n - q) by m matrix of numerics,
#' \code{Y32} a (n - q) by (p - m) matrix of numerics, and \code{k} a positive integer
#' that is the rank of the factor model if assuming a factor model.
#' \code{fa_func} also needs to return \code{Y22_array} a three-way array whose first two dimensions are
#' k and (p - m).
#'
#' @param fa_func The factor analysis function.
#' @param fa_args Additional arguments to pass to fa_func.
#'
#' @return A list with some or all of the following elements:
#'
#'     \code{ok} A logical. A value of \code{TRUE} means the function is OK. A value of
#'     \code{FALSE} means the function is not OK.
#'
#'     \code{why} A string vector giving the reason(s) why the function is not OK if \code{ok = FALSE}.
#'
#'     \code{input} A list of input values to help you debug.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{ruvb}}
#'
fa_tester_ruvb <- function(fa_func, fa_args = list()) {

  assertthat::assert_that(is.list(fa_args))
  assertthat::assert_that(is.null(fa_args$Y21))
  assertthat::assert_that(is.null(fa_args$Y31))
  assertthat::assert_that(is.null(fa_args$Y32))
  assertthat::assert_that(is.null(fa_args$k))

  if (!is.function(fa_func)) {
    return(list(ok = FALSE, why = "fa_func needs to be a function"))
  }

  ## Generate data -----------------------------------------------------------
  n <- 11
  p <- 23
  r <- 5
  k <- 3
  m <- 7
  left_mat <- matrix(stats::rnorm(n * r), nrow = n)
  right_mat <- matrix(stats::rnorm(r * p), nrow = r)
  Y <- left_mat %*% right_mat + matrix(stats::rnorm(n * p, sd = 1/2), nrow = n)
  Y21 <- Y[1:k, 1:m]
  Y22 <- Y[1:k, (m + 1):p]
  Y31 <- Y[(k + 1):n, 1:m]
  Y32 <- Y[(k + 1):n, (m + 1):p]

  fa_args_1 <- fa_args
  fa_args_1$Y21 <- Y21
  fa_args_1$Y31 <- Y31
  fa_args_1$Y32 <- Y32
  fa_args_1$k   <- k
  tryout <- tryCatch( {
    faout <- do.call(what = fa_func, args = fa_args_1)
  }, error = function(e) {NULL})
  if (is.null(tryout)) {
    return(list(ok = FALSE, why = "fa_func gave error using input", input = fa_args_1))
  }

  if (is.null(faout$Y22_array)) {
    return(list(ok = FALSE, why = "Y22_array is missing from returns", input = fa_args_1))
  }
  if (!is.array(faout$Y22_array)) {
    return(list(ok = FALSE, why = "Y22_array is not an array", input = fa_args_1))
  }
  if (length(dim(faout$Y22_array)) != 3) {
    return(list(ok = FALSE, why = "Y22_array needs to have 3 dimensions", input = fa_args_1))
  }
  if (dim(faout$Y22_array)[1] != nrow(Y21)) {
    return(list(ok = FALSE, why = "first dimension of Y22_array needs to equal nrow(Y21)", input = fa_args_1))
  }
  if (dim(faout$Y22_array)[2] != ncol(Y32)) {
    return(list(ok = FALSE, why = "second dimension of Y22_array needs to equal ncol(Y32)", input = fa_args_1))
  }

  return(list(ok = TRUE))
}
