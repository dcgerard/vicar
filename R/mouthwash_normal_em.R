###############################################################################
## Normal Mixtures ------------------------------------------------------------
###############################################################################

#' Wrapper for \code{\link{normal_mix_llike}} so that I can use it in SQUAREM.
#'
#' This returns the negative log-liklihood for SQUAREM to use.
#'
#' @inheritParams normal_mix_llike
#' @inheritParams mouthwash_second_step
#' @param scale_var Not used. Only included in the arguments for
#'     SQUAREM application.
#' @param pizxi_vec A vector whose first M elements are pi_vals, next
#'     k elements are z2, and last element is xi.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{normal_mix_llike}}
normal_mix_llike_wrapper <- function(pizxi_vec, betahat_ols, S_diag, alpha_tilde, tau2_seq,
                                     lambda_seq, scale_var = TRUE, var_inflate_pen = 0) {

    ## Make sure input is correct
    p <- length(betahat_ols)
    M <- length(tau2_seq)
    k <- ncol(alpha_tilde)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::assert_that(is.logical(scale_var))
    assertthat::are_equal(length(pizxi_vec), M + k + 1)

    pi_vals <- pizxi_vec[1:M]
    z2      <- pizxi_vec[(M + 1):(M + k)]
    xi      <- pizxi_vec[length(pizxi_vec)]
    nout <- normal_mix_llike(pi_vals = pi_vals, z2 = z2, xi = xi, betahat_ols = betahat_ols,
                             S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                             lambda_seq = lambda_seq, var_inflate_pen = var_inflate_pen)
    return(-1 * nout)
}

#' Penalized MOUTHWASH likelihood when using a mixture of normals.
#'
#' @inheritParams mouthwash_second_step
#' @param z2 The current value of the unobserved confounders
#'     corresponding to the covariates of interest.
#' @param xi The current value of the variance inflation parameter.
#' @param tau2_seq The grid of variances. This is the same thing as
#'     grid_seq in \code{\link{mouthwash_second_step}} if
#'     \code{likelihood = "normal"}.
#' @param pi_vals The current values of the mixing proportions.
#'
#' @author David Gerard
normal_mix_llike <- function(pi_vals, z2, xi, betahat_ols, S_diag, alpha_tilde, tau2_seq,
                             lambda_seq, var_inflate_pen = 0) {
    ## Make sure input is correct ------------------------------------------------------
    M <- length(tau2_seq)
    zero_spot <- which(abs(tau2_seq) < 10 ^ -12)
    assertthat::are_equal(length(zero_spot), 1)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(length(pi_vals), M)
    assertthat::assert_that(all(lambda_seq >= 1))
    assertthat::assert_that(all(tau2_seq >= 0))
    assertthat::are_equal(length(betahat_ols), nrow(alpha_tilde))
    assertthat::are_equal(length(z2), ncol(alpha_tilde))
    assertthat::assert_that(xi > 0)
    assertthat::are_equal(length(S_diag), nrow(alpha_tilde))

    ## make sure hard boundary is observed -----------------------------------------
    if (any(pi_vals < -10 ^ -12) | xi < 10 ^ -12) {
      return(-Inf)
    }

    ## Mixing variances and means ------------------------------------------------------
    mix_var  <- outer(xi * S_diag, tau2_seq, FUN = `+`) ## p by M
    mix_mean <- matrix(rep(alpha_tilde %*% z2, M), ncol = M)
    mix_obs  <- matrix(rep(betahat_ols, M), ncol = M)

    ## Log densities ------------------------------------------------------------------
    ldmix <- stats::dnorm(x = mix_obs, mean = mix_mean, sd = sqrt(mix_var), log = TRUE)
    ldmax <- apply(ldmix, 1, max)

    obs_like <- rowSums(sweep(x = exp(ldmix - ldmax), MARGIN = 2, STATS = pi_vals, FUN = `*`))
    llike <- sum(log(obs_like) + ldmax)

    if (all(lambda_seq == 1)) {
        pen <- 0
    } else {
        pen <- sum(log(pi_vals[lambda_seq > 1]) * (lambda_seq[lambda_seq > 1] - 1))
    }

    vpen <- -var_inflate_pen / xi

    return(llike + pen + vpen)
}

#' Wrapper for \code{\link{normal_mix_fix}} so that I can use SQUAREM.
#'
#' @inheritParams normal_mix_fix
#' @inheritParams mouthwash_second_step
#' @inheritParams normal_mix_llike
#' @param pizxi_vec A vector whose first M elements are pi_vals, next
#'     k elements are z2, and last element is xi.
#'
#' @seealso \code{\link{normal_mix_fix}}
#'
#' @author David Gerard
normal_mix_fix_wrapper <- function(pizxi_vec, betahat_ols, S_diag, alpha_tilde, tau2_seq,
                                   lambda_seq, scale_var = TRUE, var_inflate_pen = 0) {

    ## Make sure input is correct
    p <- length(betahat_ols)
    M <- length(tau2_seq)
    k <- ncol(alpha_tilde)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::assert_that(is.logical(scale_var))
    assertthat::are_equal(length(pizxi_vec), M + k + 1)

    pi_vals <- pizxi_vec[1:M]
    z2      <- pizxi_vec[(M + 1):(M + k)]
    xi      <- pizxi_vec[length(pizxi_vec)]
    nout <- normal_mix_fix(pi_vals = pi_vals, z2 = z2, xi = xi, betahat_ols = betahat_ols,
                           S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                           lambda_seq = lambda_seq, scale_var = scale_var,
                           var_inflate_pen = var_inflate_pen)

    pizxi_new <- c(nout$pi_vals, nout$z2, nout$xi)
    return(pizxi_new)
}

#' A fixed point iteration for updating the mixing proportions and the
#' confounders associated with the covariates of interest when using a
#' mixture of normals prior.
#'
#' @inheritParams normal_mix_llike
#' @param scale_var Should we optimize over a variance inflation
#'     parameter (\code{TRUE}) or not (\code{FALSE})?
#'
#' @return A list with the following elements.
#'
#'     \code{pi_vals}: The update for \code{pi_vals}.
#'
#'     \code{z2}: The update for \code{z2}.
#'
#'     \code{xi}: The update for \code{xi}.
#'
#' @author David Gerard
normal_mix_fix <- function(pi_vals, z2, xi, betahat_ols, S_diag, alpha_tilde, tau2_seq,
                           lambda_seq, scale_var = TRUE, var_inflate_pen = 0) {
    ## Make sure input is correct ------------------------------------------------------
    M <- length(tau2_seq)
    p <- length(betahat_ols)
    zero_spot <- which(abs(tau2_seq) < 10 ^ -8)
    assertthat::are_equal(length(zero_spot), 1)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(length(pi_vals), M)
    assertthat::assert_that(all(lambda_seq >= 1))
    assertthat::assert_that(all(tau2_seq >= 0))
    assertthat::are_equal(p, nrow(alpha_tilde))
    assertthat::are_equal(length(z2), ncol(alpha_tilde))
    assertthat::assert_that(xi > 0)
    assertthat::are_equal(length(S_diag), p)

    ## Mixing variances and means ------------------------------------------------------
    mix_var  <- outer(xi * S_diag, tau2_seq, FUN = `+`) ## p by M
    mix_mean <- matrix(rep(alpha_tilde %*% z2, M), ncol = M)
    mix_obs  <- matrix(rep(betahat_ols, M), ncol = M)

    ## Log densities ------------------------------------------------------------------
    ldmix <- stats::dnorm(x = mix_obs, mean = mix_mean, sd = sqrt(mix_var), log = TRUE)
    ldmax <- apply(ldmix, 1, max)

    denom <- log(rowSums(sweep(x = exp(ldmix - ldmax), MARGIN = 2, STATS = pi_vals, FUN = `*`))) + ldmax

    lqvals_lesspi <- ldmix - denom
    qvals <- sweep(exp(lqvals_lesspi), MARGIN = 2, STATS = pi_vals, FUN = `*`)
    assertthat::assert_that(all(abs(rowSums(qvals) - 1) < 10 ^ -8))

    ## dtemp <- stats::dnorm(x = mix_obs, mean = mix_mean, sd = sqrt(mix_var)) %*% diag(pi_vals)
    ## assertthat::are_equal(dtemp / rowSums(dtemp), qvals)

    pi_new <- (colSums(qvals) + lambda_seq - 1) / (p - M + sum(lambda_seq))
    assertthat::assert_that(abs(sum(pi_new) - 1) < 10 ^ -8)
    assertthat::assert_that(all(pi_new > -10 ^ -8))

    if (any(pi_new < 0)) {
        pi_new[pi_new < 0] <- 0
        pi_new <- pi_new / sum(pi_new)
    }

    theta_diag <- rowSums(qvals / mix_var) / 2


    z_new <- solve(t(alpha_tilde) %*% diag(theta_diag) %*% alpha_tilde) %*%
        t(alpha_tilde) %*% diag(theta_diag) %*% betahat_ols


    if (scale_var) {
        ## run a few iterations to jointly update xi and z
        xi_new <- xi
        for (iter_index in 1:10) {
            resid_vec <- betahat_ols - alpha_tilde %*% z_new
            xi_old <- xi_new
            oout <- stats::optim(par = xi_old, fn = brent_obj_norm, method = "Brent",
                                 qvals = qvals, S_diag = S_diag, tau2_seq = tau2_seq,
                                 resid_vec = resid_vec, lower = 10 ^ -14, upper = 10,
                                 control = list(fnscale = -1, maxit = 10),
                                 var_inflate_pen = var_inflate_pen)
            xi_new <- max(oout$par, 10 ^ -14)

            mix_var  <- outer(xi_new * S_diag, tau2_seq, FUN = `+`)
            theta_diag <- rowSums(qvals / mix_var) / 2
            z_new <- solve(t(alpha_tilde) %*% diag(theta_diag) %*% alpha_tilde) %*%
                t(alpha_tilde) %*% diag(theta_diag) %*% betahat_ols

            if (abs(xi_new / xi_old - 1) < 10 ^ -3) break
        }
    } else {
        xi_new <- xi
    }
    return(list(pi_vals = pi_new, z2 = z_new, xi = xi_new))
}

#' Objective function for updating variance inflation parameter during
#' EM for mixtures of normals prior.
#'
#' @inheritParams normal_mix_fix
#' @param qvals The probabilities from the E step. A matrix.
#' @param resid_vec The residuals from the OLS estimates and the bias
#'     term induced by the confounders.
#'
#' @author David Gerard
brent_obj_norm <- function(xi, qvals, S_diag, tau2_seq, resid_vec, var_inflate_pen = 0) {

    ## Make sure input is correct.
    p <- length(S_diag)
    M <- length(tau2_seq)
    assertthat::are_equal(length(xi), 1)
    assertthat::are_equal(nrow(qvals), p)
    assertthat::are_equal(ncol(qvals), M)
    assertthat::are_equal(length(resid_vec), p)
    assertthat::assert_that(xi > 0)

    mix_var  <- outer(xi * S_diag, tau2_seq, FUN = `+`) ## p by M
    theta_diag <- rowSums(qvals / mix_var)
    obj <- sum(theta_diag * resid_vec ^ 2) + sum(qvals * log(mix_var))
    obj <- -obj / 2
    vpen <- -var_inflate_pen / xi
    return(obj + vpen)
}
