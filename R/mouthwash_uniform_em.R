################################################################################
## Uniform mixtures ------------------------------------------------------------
################################################################################


#' Fixed point iteration when using t errors and a mixture of uniforms prior.
#'
#' @param pi_vals A vector of non-negative numerics that sums to
#'     one. The current values of the mixing proportions.
#' @param z2 A vector of numerics. The current values of the
#'     unobserved confounders.
#' @param xi A positive numeric. The current value of the variance
#'     inflation factor.
#' @param betahat_ols A vector of numerics. The OLS regression
#'     coefficients.
#' @param S_diag A vector of positive numerics. The standard errors.
#' @param alpha_tilde A matrix of numerics. The estimates of the
#'     coefficients of the confounders.
#' @param a_seq A vector of numerics (usually non-positive). The lower
#'     bounds of the uniform mixing densities.
#' @param b_seq A vector of numerics (usually non-negative). The upper
#'     bounds of the uniform mixing densities.
#' @param lambda_seq A vector of numerics greater than 1. The
#'     penalties on the mixing proportions.
#' @param degrees_freedom A positive scalar or a vector of positive
#'     scalars the length of \code{betahat_ols}. The degrees of
#'     freedom of the t-distribution.
#' @param scale_var A logical. Should we update the variance inflation
#'     parameter (\code{TRUE}) or not (\code{FALSE})
#'
#' @author David Gerard
uniform_mix_fix <- function(pi_vals, z2, xi, betahat_ols, S_diag, alpha_tilde, a_seq, b_seq,
                             lambda_seq, degrees_freedom, scale_var = TRUE) {

    ## Make sure input is correct ----------------------------------------------
    zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)

    M <- length(pi_vals)
    p <- length(betahat_ols)
    k <- length(z2)
    assertthat::are_equal(length(a_seq), M)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(ncol(alpha_tilde), k)

    ## Calculate f_tildes -----------------------------------------------------
    sd_mat <- matrix(rep(sqrt(xi * S_diag), M), ncol = M)
    resid_vec <- betahat_ols - alpha_tilde %*% z2
    left_centered_means  <- outer(c(resid_vec), a_seq, FUN = `-`) / sd_mat
    right_centered_means <- outer(c(resid_vec), b_seq, FUN = `-`) / sd_mat

    denom_mat <- matrix(rep(b_seq - a_seq, p), byrow = TRUE, ncol = M)
    top_mat <- ptdiff_mat(X = left_centered_means, Y = right_centered_means,
                          degrees_freedom = degrees_freedom)

    ## top_mat2 <- stats::pt(q = left_centered_means, df = degrees_freedom) -
    ##     stats::pt(q = right_centered_means, df = degrees_freedom)
    ## assertthat::are_equal(top_mat, top_mat2)

    ftilde_mat <- top_mat / denom_mat

    ftilde_mat[, zero_spot] <- dt_wrap(x = betahat_ols, df = degrees_freedom,
                                       mean = alpha_tilde %*% z2,
                                       sd = sqrt(xi * S_diag))

    pif <- sweep(ftilde_mat, MARGIN = 2, STATS = pi_vals, FUN = `*`)
    qvals <- pif / rowSums(pif)
    assertthat::assert_that(all(abs(rowSums(qvals) - 1) < 10 ^ -14))

    pi_new <- (colSums(qvals) + lambda_seq - 1) / (p - M + sum(lambda_seq))
    assertthat::assert_that(abs(sum(pi_new) - 1) < 10 ^ -14)
    assertthat::assert_that(all(pi_new > -10 ^ -8))

    if (any(pi_new < 0)) {
      pi_new[pi_new < 0] <- 0
      pi_new <- pi_new / sum(pi_new)
    }


    oout <- stats::optim(par = z2, fn = unif_int_obj, gr = unif_int_grad,
                         method = "BFGS",
                         control = list(fnscale = -1, maxit = 10),
                         xi = xi, betahat_ols = betahat_ols, alpha_tilde = alpha_tilde,
                         S_diag = S_diag, a_seq = a_seq, b_seq = b_seq, qvals = qvals,
                         degrees_freedom = degrees_freedom)

    z_new <- oout$par
    xi_new <- xi
    if (scale_var) {
        for (index in 1:10) {
            xi_old <- xi_new
            oout <- stats::optim(par = xi_new, fn = unif_int_obj, method = "Brent",
                                 lower = 10 ^ -14, upper = 10,
                                 control = list(fnscale = -1, maxit = 10),
                                 z2 = z_new, betahat_ols = betahat_ols,
                                 alpha_tilde = alpha_tilde, S_diag = S_diag, a_seq = a_seq,
                                 b_seq = b_seq, qvals = qvals, degrees_freedom = degrees_freedom)
            xi_new <- max(oout$par, 10 ^ -14)
            oout <- stats::optim(par = z_new, fn = unif_int_obj, gr = unif_int_grad,
                                 method = "BFGS",
                                 control = list(fnscale = -1, maxit = 10),
                                 xi = xi_new, betahat_ols = betahat_ols, alpha_tilde = alpha_tilde,
                                 S_diag = S_diag, a_seq = a_seq, b_seq = b_seq, qvals = qvals,
                                 degrees_freedom = degrees_freedom)
            z_new <- oout$par
            if (abs(xi_old / xi_new - 1) < 10 ^ -3) break
        }
    }
    return(list(pi_vals = pi_new, z2 = z_new, xi = xi_new))
}

#' Intermediate objective function in EM.
#'
#' @inheritParams uniform_mix_fix
#' @param qvals A matrix of non-negative numerics whose rows sum to
#'     one. The E-step proportions.
#'
#' @author David Gerard
unif_int_obj <- function(z2, xi, betahat_ols, alpha_tilde, S_diag, a_seq, b_seq, qvals,
                         degrees_freedom) {

    ## Make sure input is correct --------------------------------------------
    M <- length(a_seq)
    p <- length(betahat_ols)
    k <- length(z2)
    zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(ncol(alpha_tilde), k)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(ncol(qvals), M)
    assertthat::are_equal(nrow(qvals), p)

    ## Calculate objective function ------------------------------------------
    sd_mat <- matrix(rep(sqrt(xi * S_diag), M), ncol = M)
    resid_vec <- betahat_ols - alpha_tilde %*% z2
    left_centered_means  <- outer(c(resid_vec), a_seq, FUN = `-`) / sd_mat
    right_centered_means <- outer(c(resid_vec), b_seq, FUN = `-`) / sd_mat

    tdiff_mat <- ptdiff_mat_log(X = left_centered_means, Y = right_centered_means,
                                degrees_freedom = degrees_freedom)
    ## tdiff_mat2 <- stats::pt(q = left_centered_means, df = degrees_freedom) -
    ##     stats::pt(q = right_centered_means, df = degrees_freedom)
    ## assertthat::are_equal(tdiff_mat, tdiff_mat2)

    tdiff_mat[, zero_spot] <- dt_wrap(x = betahat_ols, mean = alpha_tilde %*% z2,
                                      sd = sqrt(xi * S_diag), df = degrees_freedom, log = TRUE)

    ## tdiff_mat2 <- ptdiff_mat(X = left_centered_means, Y = right_centered_means,
    ##                          degrees_freedom = degrees_freedom)
    ## tdiff_mat2[, zero_spot] <- dt_wrap(x = betahat_ols, mean = alpha_tilde %*% z2,
    ##                                   sd = sqrt(xi * S_diag), df = degrees_freedom)
    ## sum(log(tdiff_mat2) * qvals)


    obj_val <- sum(tdiff_mat * qvals)
    return(obj_val)
}

#' Gradient wrt z2 of intermediate function in EM.
#'
#' @inheritParams uniform_mix_fix
#' @param qvals A matrix of non-negative numerics whose rows sum to
#'     one. The E-step proportions.
#'
#' @author David Gerard
unif_int_grad <- function(z2, xi, betahat_ols, alpha_tilde, S_diag, a_seq, b_seq, qvals,
                          degrees_freedom) {

    ## Make sure input is correct --------------------------------------------
    M <- length(a_seq)
    p <- length(betahat_ols)
    k <- length(z2)
    zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(ncol(alpha_tilde), k)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(ncol(qvals), M)
    assertthat::are_equal(nrow(qvals), p)

    ## Calculate gradient ----------------------------------------------------
    sd_mat <- matrix(rep(sqrt(xi * S_diag), M), ncol = M)
    resid_vec <- betahat_ols - alpha_tilde %*% z2
    left_centered_means  <- outer(c(resid_vec), a_seq, FUN = `-`) / sd_mat
    right_centered_means <- outer(c(resid_vec), b_seq, FUN = `-`) / sd_mat

    tdiff_mat <- ptdiff_mat(X = left_centered_means, Y = right_centered_means,
                            degrees_freedom = degrees_freedom)

    ## tdiff_mat2 <- stats::pt(q = left_centered_means, df = degrees_freedom) -
    ##     stats::pt(q = right_centered_means, df = degrees_freedom)
    ## assertthat::are_equal(tdiff_mat, tdiff_mat2)


    dtdiff_mat <- stats::dt(x = right_centered_means, df = degrees_freedom) / sd_mat -
        stats::dt(x = left_centered_means, df = degrees_freedom) / sd_mat

    tratio_mat <- dtdiff_mat / tdiff_mat

    ## stupid hack to get rid of NaN's --------------------------------------
    tratio_mat[tdiff_mat == 0] <- 0
    ## ----------------------------------------------------------------------

    tratio_mat[, zero_spot] <- (degrees_freedom + 1) * resid_vec /
        (degrees_freedom * xi * S_diag + resid_vec ^ 2)

    grad_val <- crossprod(alpha_tilde, rowSums(tratio_mat * qvals))

    return(grad_val)
}

#' Wrapper for \code{\link{uniform_mix_fix}}, mostly so I can use SQUAREM.
#'
#' @inheritParams uniform_mix_fix
#' @param pizxi_vec A vector of numerics. The first M of which are
#'     \code{pi_vals}, the next k of which are \code{z2}, the last
#'     element is \code{xi}.
#'
#' @author David Gerard
uniform_mix_fix_wrapper <- function(pizxi_vec, betahat_ols, S_diag, alpha_tilde, a_seq, b_seq,
                                      lambda_seq, degrees_freedom, scale_var = TRUE) {
    ## Make sure input is correct
    p <- length(betahat_ols)
    M <- length(a_seq)
    k <- ncol(alpha_tilde)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::assert_that(is.logical(scale_var))
    assertthat::are_equal(length(pizxi_vec), M + k + 1)

    pi_vals <- pizxi_vec[1:M]
    z2      <- pizxi_vec[(M + 1):(M + k)]
    xi      <- pizxi_vec[length(pizxi_vec)]
    fout <- uniform_mix_fix(pi_vals = pi_vals, z2 = z2, xi = xi,
                            betahat_ols = betahat_ols, S_diag = S_diag,
                            alpha_tilde = alpha_tilde, a_seq = a_seq,
                            b_seq = b_seq, lambda_seq = lambda_seq,
                            degrees_freedom = degrees_freedom,
                            scale_var = scale_var)
    pizxi_new <- c(fout$pi_vals, fout$z2, fout$xi)
    return(pizxi_new)
}

#' Log-likelihood when using t errors and mixture of uniforms prior.
#'
#' @inheritParams uniform_mix_fix
#' @param var_inflate_pen The penalty to apply on the variance inflation parameter.
#'     Defaults to 0, but should be something non-zero when \code{alpha = 1}
#'     and \code{scale_var = TRUE}.
#'
#' @author David Gerard
uniform_mix_llike <- function(pi_vals, z2, xi, betahat_ols, S_diag, alpha_tilde, a_seq, b_seq,
                              lambda_seq, degrees_freedom, var_inflate_pen = 0) {

    ## Make sure input is correct ----------------------------------------------
    zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)

    M <- length(pi_vals)
    p <- length(betahat_ols)
    k <- length(z2)
    assertthat::are_equal(length(a_seq), M)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(ncol(alpha_tilde), k)

    ## make sure hard boundary is observed -----------------------------------------
    if (any(pi_vals < -10 ^ -14) | xi < 10 ^ -14) {
      return(-Inf)
    }

    ## Calculate f_tildes -----------------------------------------------------
    sd_mat <- matrix(rep(sqrt(xi * S_diag), M), ncol = M)
    resid_vec <- betahat_ols - alpha_tilde %*% z2
    left_centered_means  <- outer(c(resid_vec), a_seq, FUN = `-`) / sd_mat
    right_centered_means <- outer(c(resid_vec), b_seq, FUN = `-`) / sd_mat

    denom_mat <- matrix(rep(b_seq - a_seq, p), byrow = TRUE, ncol = M)

    top_mat <- ptdiff_mat(X = left_centered_means, Y = right_centered_means,
                            degrees_freedom = degrees_freedom)

    ## top_mat2 <- stats::pt(q = left_centered_means, df = degrees_freedom) -
    ##     stats::pt(q = right_centered_means, df = degrees_freedom)
    ## assertthat::are_equal(top_mat, top_mat2)

    ftilde_mat <- top_mat / denom_mat

    ftilde_mat[, zero_spot] <- dt_wrap(x = betahat_ols, df = degrees_freedom,
                                       mean = alpha_tilde %*% z2,
                                       sd = sqrt(xi * S_diag))

    llike <- sum(log(rowSums(sweep(ftilde_mat, MARGIN = 2, STATS = pi_vals, FUN = `*`))))

    if (all(lambda_seq == 1)) {
        pen <- 0
    } else {
        pen <- sum(log(pi_vals[lambda_seq > 1]) * (lambda_seq[lambda_seq > 1] - 1))
    }

    vpen <- -var_inflate_pen / xi

    return(llike + pen + vpen)
}

#' Wrapper for \code{\link{uniform_mix_llike}}, mostly for the SQUAREM package.
#'
#' This returns the negative log-likelihood for use in SQUAREM.
#'
#' @inheritParams uniform_mix_fix
#' @param pizxi_vec A vector of numerics. The first M of which are
#'     \code{pi_vals}, the next k of which are \code{z2}, the last
#'     element is \code{xi}.
#'
#' @author David Gerard
uniform_mix_llike_wrapper <- function(pizxi_vec, betahat_ols, S_diag, alpha_tilde, a_seq, b_seq,
                                      lambda_seq, degrees_freedom, scale_var = TRUE) {
    ## Make sure input is correct
    p <- length(betahat_ols)
    M <- length(a_seq)
    k <- ncol(alpha_tilde)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::assert_that(is.logical(scale_var))
    assertthat::are_equal(length(pizxi_vec), M + k + 1)

    pi_vals <- pizxi_vec[1:M]
    z2      <- pizxi_vec[(M + 1):(M + k)]
    xi      <- pizxi_vec[length(pizxi_vec)]
    llike <- uniform_mix_llike(pi_vals = pi_vals, z2 = z2, xi = xi,
                             betahat_ols = betahat_ols, S_diag = S_diag,
                             alpha_tilde = alpha_tilde, a_seq = a_seq,
                             b_seq = b_seq, lambda_seq = lambda_seq,
                             degrees_freedom = degrees_freedom)
    return(-1 * llike)
}


#' More stable way to calculate differences in T cdf's.
#'
#' If \eqn{T()} is the cdf of a t-distribution if
#' \code{degrees_freedom} degrees of freedom, then this function will
#' calculate \eqn{T(x) - T(y)} or \eqn{T(-y) - T(-x)} (which are are
#' equivalent), depending on which one is more numerically stable.
#'
#' @param x A quantile of the t.
#' @param y The other quantile of the t.
#' @param degrees_freedom The degrees of freedom of the t.
#'
#' @author David Gerard
ptdiff <- function(x, y, degrees_freedom) {
    if (abs(x) > abs(y)) {
        diff <- stats::pt(x, df = degrees_freedom) - stats::pt(y, df = degrees_freedom)
    } else {
        diff <- stats::pt(-1 * y, df = degrees_freedom) - stats::pt(-1 * x, df = degrees_freedom)
    }
    return(diff)
}

#' More stable way to calculate differences in T cdfs when input is a matrix.
#'
#' @param X A matrix of quantiles of the t.
#' @param Y Another matrix of quantiles of the t.
#' @param degrees_freedom The degrees of freedom of the t.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{ptdiff}}.
ptdiff_mat <- function(X, Y, degrees_freedom) {

    assertthat::are_equal(dim(X), dim(Y))

    which_switch <- abs(Y) > abs(X)

    diff <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
    diff[which_switch] <- stats::pt(-Y[which_switch], df = degrees_freedom) -
        stats::pt(-X[which_switch], df = degrees_freedom)
    diff[!which_switch] <- stats::pt(X[!which_switch], df = degrees_freedom) -
        stats::pt(Y[!which_switch], df = degrees_freedom)
    return(diff)
}


#' Log version of \code{\link{ptdiff_mat}}.
#'
#' @inheritParams ptdiff_mat
#'
#' @author David Gerard
#'
#' @seealso \code{\link{ptdiff_mat}}.
ptdiff_mat_log <- function(X, Y, degrees_freedom) {
    assertthat::are_equal(dim(X), dim(Y))
    which_switch <- abs(Y) > abs(X)
    ldiff_mat <- matrix(NA, nrow = nrow(X), ncol = ncol(X))

    a <- stats::pt(X[!which_switch], df = degrees_freedom, log.p = TRUE)
    b <- stats::pt(Y[!which_switch], df = degrees_freedom, log.p = TRUE)
    ldiff_mat[!which_switch] <- a + log(1 - exp(b - a))

    a <- stats::pt(-X[which_switch], df = degrees_freedom, log.p = TRUE)
    b <- stats::pt(-Y[which_switch], df = degrees_freedom, log.p = TRUE)
    ldiff_mat[which_switch] <- b + log(1 - exp(a - b))
    return(ldiff_mat)
}
