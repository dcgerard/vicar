## Functions for performing MOUTHWASH. This is a re-coding from
## succotashr because I'm slightly better at programming now.


#' MOUTHWASH: Maximize Over Unobservables To Help With Adaptive SHrinkage.
#'
#'
#'
#' @inheritParams vruv4
#' @param grid_seq The grid for the mixing distribution. If
#'     \code{mixing_dist = "uniform"} or \code{"+uniform"}, then these
#'     should be the non-zero limits of the uniform distributions. If
#'     \code{mixing_dist = "sym_uniform"}, then these should be the
#'     right limits of the uniform distributions. If \code{mixing_dist
#'     = "normal"}, then these should be the variances of the mixing
#'     normal distributions.
#' @param pi_init A numeric vector. These are the initial values of
#'     the mixing proportions.
#' @param lambda_seq A numeric vector with elements all greater than
#'     or equal to 1. These are the tuning parameters for the mixing
#'     proportions. This can only be specified if \code{grid_seq} is
#'     also specified.
#' @param mixing_dist A character. Should we use a mixture of uniforms
#'     (\code{"uniform"}), a mixture of uniforms with minimum at 0
#'     (\code{"+uniform"}), a mixture of uniforms symmetric at 0
#'     (\code{"sym_uniform"}), or a mixture of normals
#'     (\code{"normal"})?
#' @param lambda_type A character. Should we apply a penalty on zero
#'     (\code{"zero_conc"}) or no penalty (\code{"uniform"}). Not used
#'     if \code{lambda_seq} is not \code{NULL}.
#' @param lambda0 A numeric greater than or equal to 1. The penalty on
#'     zero if \code{lambda_type = "zero_conc"}.
#' @param pi_init_type How should we initialize the mixing
#'     proportions? By concentrating on zero (\code{"zero_conc"}), by
#'     equal weights on all mixing distributions (\code{"uniform"}),
#'     or by sampling uniformly on the simplex (\code{"random"})?
#' @param scale_var A logical. Should we estimate a variance inflation
#'     parameter (\code{TRUE}) or not (\code{FALSE})?
#'
#' @export
#'
#' @author David Gerard
mouthwash <- function(Y, X, k = NULL, cov_of_interest = ncol(X),
                      include_intercept = TRUE, limmashrink = TRUE,
                      fa_func = pca_naive, fa_args = list(),
                      likelihood = c("t", "normal"),
                      mixing_dist = c("uniform", "+uniform", "sym_uniform", "normal"),
                      lambda_type = c("zero_conc", "uniform"),
                      pi_init_type = c("zero_conc", "uniform", "random"),
                      degrees_freedom = NULL,
                      pi_init = NULL, grid_seq = NULL, lambda_seq = NULL,
                      lambda0 = 10,
                      scale_var = TRUE) {

    ## Make sure input is correct --------------------------------------------------------
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(is.numeric(cov_of_interest))
    assertthat::are_equal(length(cov_of_interest), 1)
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.function(fa_func))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(lambda0 >= 1)

    likelihood   <- match.arg(likelihood)
    mixing_dist  <- match.arg(mixing_dist)
    pi_init_type <- match.arg(pi_init_type)
    lambda_type  <- match.arg(lambda_type)

    if (likelihood == "normal" & mixing_dist == "normal") {
        stop("normal mixtures not implemented for t-likelihood yet (or likely ever).")
    }

    ## Rotate ---------------------------------------------------------------------------
    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink, fa_func = fa_func,
                               fa_args = fa_args, do_factor = TRUE)


    ## Deal with degrees of freedom -----------------------------------------------------
    if (likelihood == "normal") {
        if (!is.null(degrees_freedom)) {
            message("likelihood = \"normal\" but degrees_freedom not NULL. Setting degrees_freedom to Inf")
        }
        degrees_freedom <- Inf
    }
    if (!is.null(rotate_out$prior_df) & is.null(degrees_freedom) & likelihood == "t") {
        degrees_freedom <- rotate_out$prior_df + nrow(X) - ncol(X) - k
        if (degrees_freedom == Inf) {
            message("limma estimated df = Inf . Changing likelihood to \"normal\".")
            likelihood <- "normal"
        }
    } else if (is.null(degrees_freedom) & likelihood == "t") {
        degrees_freedom <- nrow(X) - ncol(X) - k
    }
    assertthat::assert_that(length(degrees_freedom) == 1 | length(degrees_freedom) == ncol(Y))
    assertthat::assert_that(all(degrees_freedom > 0))


    ## rescale alpha and sig_diag by R22 to get data for second step ------------------
    alpha_tilde <- rotate_out$alpha / c(rotate_out$R22)
    S_diag      <- rotate_out$sig_diag / c(rotate_out$R22 ^ 2)
    betahat_ols <- rotate_out$betahat_ols

    ## Set grid and penalties ---------------------------------------------------------
    if (!is.null(lambda_seq) & is.null(grid_seq)) {
        stop("lambda_seq specified but grid_seq is NULL")
    } else if (!is.null(lambda_seq) & !is.null(grid_seq)) {
        M <- length(lambda_seq)
        assertthat::are_equal(length(lambda_seq), length(grid_seq))
    }

    if (is.null(grid_seq)) {
        grid_vals <- get_grid_var(betahat_ols = betahat_ols, S_diag = S_diag)
        if (mixing_dist == "normal") {
            grid_seq <- sign(grid_vals$tau2_seq) * sqrt(abs(grid_vals$tau2_seq))
        } else {
            grid_seq <- grid_vals$tau2_seq
        }
        M <- length(grid_seq)
    }

    zero_spot <- which(abs(grid_seq) < 10 ^ -14) ## the location of the zero
    assertthat::are_equal(length(zero_spot), 1)

    if (is.null(lambda_seq)) {
        if (lambda_type == "uniform") {
            lambda_seq <- rep(1, M)
        } else if (lambda_type == "zero_conc") {
            lambda_seq <- rep(1, M)
            lambda_seq[zero_spot] <- lambda0
        }
    }

    result_out <- mouthwash_second_step(betahat_ols = betahat_ols, S_diag = S_diag,
                                        alpha_tilde = alpha_tilde, grid_seq = grid_seq,
                                        lambda_seq = lambda_seq, mixing_dist = mixing_dist,
                                        likelihood = likelihood, pi_init_type = pi_init_type)

}

#' The second step of MOUTHWASH.
#'
#' @inheritParams mouthwash
#' @param betahat_ols A vector of numerics. The OLS estimates of the
#'     coefficients of interest.
#' @param S_diag A vector of positive numerics. The estimated standard
#'     errors.
#' @param alpha_tilde A matrix. The number of rows should be equal the
#'     length of betahat_ols. The number of columns should equal the
#'     number of hidden confounders.
#' @param grid_seq The grid of variances (if \code{mixing_dist =
#'     "normal"}) or the grid of uniform bounds. See
#'     \code{\link{mouthwash}} for details.
#'
#' @author David Gerard
#'
#' @export
#'
mouthwash_second_step <- function(betahat_ols, S_diag, alpha_tilde, grid_seq, lambda_seq,
                                  mixing_dist = c("uniform", "+uniform", "sym_uniform", "normal"),
                                  likelihood = c("t", "normal"),
                                  pi_init_type = c("zero_conc", "uniform", "random"),
                                  scale_var = TRUE) {

    ## Make sure input is correct -------------------------------------------------
    M <- length(grid_seq)
    k <- ncol(alpha_tilde)
    zero_spot <- which(abs(grid_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)

    assertthat::are_equal(length(betahat_ols), nrow(alpha_tilde))
    assertthat::are_equal(length(S_diag), length(betahat_ols))
    assertthat::are_equal(length(lambda_seq), M)

    mixing_dist  <- match.arg(mixing_dist)
    likelihood   <- match.arg(likelihood)
    pi_init_type <- match.arg(pi_init_type)

    ## initialize parameters and run EM --------------------------------------------
    z2_init <- matrix(stats::rnorm(k), ncol = 1)
    pi_init <- initialize_mixing_prop(pi_init_type = pi_init_type, zero_spot = zero_spot, M = M)
    pizxi_init <- c(pi_init, z2_init, 1)

    if (likelihood == "normal" & mixing_dist == "normal") {
        sqout <- SQUAREM::squarem(par = pizxi_init, fixptfn = normal_mix_fix_wrapper,
                                  objfn = normal_mix_llike_wrapper, betahat_ols = betahat_ols,
                                  S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = grid_seq,
                                  lambda_seq = lambda_seq, scale_var = scale_var,
                                  control = list(tol = 10 ^ -4))
    }

}

#' Function for initializing mixing proportions.
#'
#' @inheritParams mouthwash
#' @param M The number of mixing distributions.
#' @param zero_spot The location of the zero.
#'
#' @author David Gerard
initialize_mixing_prop <- function(pi_init_type = c("zero_conc", "uniform", "random"),
                                   zero_spot, M) {
    pi_init_type <- match.arg(pi_init_type)
    assertthat::assert_that(M > 0)

    if (pi_init_type == "zero_conc") {
        pi_init <- rep(0.1 / (M - 1), length = M)
        pi_init[zero_spot] <- 0.9
    } else if (pi_init_type == "uniform") {
        pi_init <- rep(1 / M, length = M)
    } else if (pi_init_type == "random") {
        pi_init <- stats::runif(n = M)
        pi_init <- pi_init / sum(pi_init)
    }
    return(pi_init)
}

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
                                     lambda_seq, scale_var = TRUE) {

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
                             lambda_seq = lambda_seq)
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
                             lambda_seq) {
    ## Make sure input is correct ------------------------------------------------------
    M <- length(tau2_seq)
    zero_spot <- which(abs(tau2_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(length(pi_vals), M)
    assertthat::assert_that(all(lambda_seq >= 1))
    assertthat::assert_that(all(tau2_seq >= 0))
    assertthat::are_equal(length(betahat_ols), nrow(alpha_tilde))
    assertthat::are_equal(length(z2), ncol(alpha_tilde))
    assertthat::assert_that(xi > 0)
    assertthat::are_equal(length(S_diag), nrow(alpha_tilde))

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

    return(llike + pen)
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
                                   lambda_seq, scale_var = TRUE) {

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
                           lambda_seq = lambda_seq, scale_var = scale_var)

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
                           lambda_seq, scale_var = TRUE) {
    ## Make sure input is correct ------------------------------------------------------
    M <- length(tau2_seq)
    p <- length(betahat_ols)
    zero_spot <- which(abs(tau2_seq) < 10 ^ -14)
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
    assertthat::assert_that(all(abs(rowSums(qvals) - 1) < 10 ^ -14))

    ## dtemp <- stats::dnorm(x = mix_obs, mean = mix_mean, sd = sqrt(mix_var)) %*% diag(pi_vals)
    ## assertthat::are_equal(dtemp / rowSums(dtemp), qvals)

    pi_new <- (colSums(qvals) + lambda_seq - 1) / (p - M + sum(lambda_seq))
    assertthat::assert_that(abs(sum(pi_new) - 1) < 10 ^ -14)

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
                                 resid_vec = resid_vec, lower = 0, upper = 10,
                                 control = list(fnscale = -1, maxit = 10))
            xi_new <- oout$par

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
brent_obj_norm <- function(xi, qvals, S_diag, tau2_seq, resid_vec) {

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
    return(obj)
}



#' Default way to set grid of variances.
#'
#' @param betahat_ols The OLS estimates of the coefficients of interest.
#' @param S_diag The standard error estimates.
#'
#' @author David Gerard
get_grid_var <- function(betahat_ols, S_diag) {

    ## Check input ----------------------------------------------------
    assertthat::are_equal(length(betahat_ols), length(S_diag))
    assertthat::assert_that(all(S_diag > 0))

    ## default grid to be same as in ASH ------------------------------
    tau2_min <- min(S_diag) / 100
    tau2_max <- 4 * max(betahat_ols ^ 2 - S_diag)
    if (tau2_max < 0) {
        tau2_max <- 64 * tau2_min
    }
    tau2_current <- tau2_min
    tau2_seq <- c(0, tau2_current)
    mult_fact <- 2
    while (tau2_current <= tau2_max) {
        tau2_current <- tau2_current * mult_fact
        tau2_seq <- c(tau2_seq, tau2_current)
    }
    M <- length(tau2_seq)
    return(list(tau2_seq = tau2_seq, M = M))
}


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
#' @param alpha_tilde A matrix of numerics. The estimats of the
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
    top_mat <- stats::pt(q = left_centered_means, df = degrees_freedom) -
        stats::pt(q = right_centered_means, df = degrees_freedom)
    ftilde_mat <- top_mat / denom_mat

    ftilde_mat[, zero_spot] <- dt_wrap(x = betahat_ols, df = degrees_freedom,
                                       mean = alpha_tilde %*% z2,
                                       sd = sqrt(xi * S_diag))

    pif <- sweep(ftilde_mat, MARGIN = 2, STATS = pi_vals, FUN = `*`)
    qvals <- pif / rowSums(pif)
    assertthat::assert_that(all(abs(rowSums(qvals) - 1) < 10 ^ -14))

    pi_new <- (colSums(qvals) + lambda_seq - 1) / (p - M + sum(lambda_seq))
    assertthat::assert_that(abs(sum(pi_new) - 1) < 10 ^ -14)


    oout <- stats::optim(par = z2, fn = unif_int_obj, gr = unif_int_grad,
                         method = "L-BFGS-B",
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
                                 lower = 0, upper = 10,
                                 control = list(fnscale = -1, maxit = 10),
                                 z2 = z_new, betahat_ols = betahat_ols,
                                 alpha_tilde = alpha_tilde, S_diag = S_diag, a_seq = a_seq,
                                 b_seq = b_seq, qvals = qvals, degrees_freedom = degrees_freedom)
            xi_new <- oout$par
            oout <- stats::optim(par = z_new, fn = unif_int_obj, gr = unif_int_grad,
                                 method = "L-BFGS-B",
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
    tdiff_mat <- stats::pt(q = left_centered_means, df = degrees_freedom) -
        stats::pt(q = right_centered_means, df = degrees_freedom)

    tdiff_mat[, zero_spot] <- dt_wrap(x = betahat_ols, mean = alpha_tilde %*% z2,
                                      sd = sqrt(xi * S_diag), df = degrees_freedom)

    obj_val <- sum(log(tdiff_mat) * qvals)
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
    tdiff_mat <- stats::pt(q = left_centered_means, df = degrees_freedom) -
        stats::pt(q = right_centered_means, df = degrees_freedom)

    dtdiff_mat <- stats::dt(x = right_centered_means, df = degrees_freedom) / sd_mat -
        stats::dt(x = left_centered_means, df = degrees_freedom) / sd_mat

    tratio_mat <- dtdiff_mat / tdiff_mat

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
#'
#' @author David Gerard
uniform_mix_llike <- function(pi_vals, z2, xi, betahat_ols, S_diag, alpha_tilde, a_seq, b_seq,
                              lambda_seq, degrees_freedom) {

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
    right_centered_means <- outer(c(resid_vec), b_seq, FUN = `-`)

    denom_mat <- matrix(rep(b_seq - a_seq, p), byrow = TRUE, ncol = M)
    top_mat <- stats::pt(q = left_centered_means, df = degrees_freedom) -
        stats::pt(q = right_centered_means, df = degrees_freedom)
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

    return(llike + pen)
}

#' Wrapper for \code{\link{uniform_mix_llike}}, mostly for the SQUAREM package.
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
    return(llike)
}


#' Wrapper for dt with a non-zero mena and non-1 scale parameter.
#'
#' @param x A numeric. Where to evaluate the density function.
#' @param df A positive numeric. The degrees of freedom.
#' @param mean A numeric. The mean of the t. Defaults to 0.
#' @param sd A positive numeric. The standard deviation of the
#'     t. Defaults to 1.
#' @param log A logical. Should we return the log-density
#'     (\code{TRUE}) or the density (\code{FALSE})?
dt_wrap <- function(x, df, mean = 0, sd = 1, log = FALSE) {
    x_new <- (x - mean) / sd
    if (!log) {
        dval <- stats::dt(x_new, df = df) / sd
    } else {
        dval <- stats::dt(x_new, df = df, log = TRUE) - log(sd)
    }
    return(dval)
}

#' Wrapper for pt with a non-zero mena and non-1 scale parameter.
#'
#' @param x A numeric. Where to evaluate the cdf.
#' @inheritParams dt_wrap
#'
pt_wrap <- function(x, df, mean = 0, sd = 1) {
    x_new <- (x - mean) / sd
    pval <- stats::pt(x_new, df = df)
    return(pval)
}
