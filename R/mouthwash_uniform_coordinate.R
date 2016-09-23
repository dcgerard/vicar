###############################################################################
## Coordinate ascent for MOUTHWASH with uniform mixtures ----------------------
###############################################################################

#' Coordinate ascent for optimizing t likelihood with uniform mixtures.
#'
#' Set \code{degrees_freedom = Inf} if you want to use a normal likelihood.
#'
#' @inheritParams mouthwash_second_step
#' @param tol The tolerance for the stopping condition.
#' @param maxit The maximum number of iterations to run in the
#'     optimization scheme.
#' @param pi_init The intial values of the mixing proportions. A
#'     numeric vector of positive values that sum to one.
#' @param z_init The initial values of z2. A numeric vector.
#' @param xi_init The initial (or known) value of the variance
#'     inflation term.
#' @param plot_update A logical. Should I plot the the path of the
#'     log-likelihood (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
mouthwash_coordinate <- function(pi_init, z_init, xi_init, betahat_ols, S_diag,
                                 alpha_tilde, a_seq, b_seq, lambda_seq,
                                 degrees_freedom, scale_var = TRUE,
                                 tol = 10 ^ -6, maxit = 100, plot_update = FALSE) {

    ## Make sure input is correct ---------------------------------------------
    M <- length(a_seq)
    k <- length(z_init)
    p <- length(betahat_ols)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(ncol(alpha_tilde), k)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(sum(pi_init), 1)
    assertthat::assert_that(all(b_seq - a_seq > -10 ^ -14))
    assertthat::assert_that(all(lambda_seq - 1 > -10 ^ -14))
    assertthat::assert_that(xi_init > 0)

    zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)



    pi_new <- pi_init
    z_new  <- z_init
    xi_new <- xi_init

    llike_new <- uniform_mix_llike(pi_vals = pi_new, z2 = z_new, xi = xi_new,
                                   betahat_ols = betahat_ols, S_diag = S_diag,
                                   alpha_tilde = alpha_tilde, a_seq = a_seq,
                                   b_seq = b_seq, lambda_seq = lambda_seq,
                                   degrees_freedom = degrees_freedom)

    llike_vec <- llike_new
    err <- tol + 1
    iter_index <- 1
    while (err > tol & iter_index < maxit) {
        llike_old <- llike_new
        pi_old    <- pi_new
        z_old     <- z_new
        xi_old    <- xi_new

        ## Update z with quasi-Newton ----------------------------------------------------
        optim_out <- stats::optim(par = z_new, fn = uniform_mix_llike,
                                  gr = mouthwash_z_grad, method = "BFGS",
                                  pi_vals = pi_new, xi = xi_new, betahat_ols = betahat_ols,
                                  S_diag = S_diag, alpha_tilde = alpha_tilde, a_seq = a_seq,
                                  b_seq = b_seq, lambda_seq = lambda_seq,
                                  degrees_freedom = degrees_freedom,
                                  control = list(fnscale = -1))

        z_new <- optim_out$par

        ## update xi with Brent if scale_var = TRUE----------------------------------------
        if (scale_var) {
            optim_out <- stats::optim(par = xi_new, fn = uniform_mix_llike, method = "Brent",
                                      lower = 10 ^ -14, upper = 10,
                                      pi_vals = pi_new, z2 = z_new, betahat_ols = betahat_ols,
                                      S_diag = S_diag, alpha_tilde = alpha_tilde, a_seq = a_seq,
                                      b_seq = b_seq, lambda_seq = lambda_seq,
                                      degrees_freedom = degrees_freedom,
                                      control = list(fnscale = -1))

            xi_new <- optim_out$par
        }

        ## update pi_vals with ashr ------------------------------------------------------
        ash_data <- ashr::set_data(betahat = c(betahat_ols - alpha_tilde %*% z_new),
                                   sebetahat = c(sqrt(xi_new * S_diag)),
                                   lik = ashr::t_lik(degrees_freedom),
                                   alpha = 0)
        ash_g <- ashr::unimix(pi = pi_new, a = a_seq, b = b_seq)

        if (requireNamespace(package = "Rmosek", quietly = TRUE)) {
            optmethod <- "mixIP"
        } else {
            optmethod <- "mixEM"
        }
        control_default <- list(K = 1, method = 3, square = TRUE,
                                step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                                tol = 1e-07, maxiter = 500, trace = FALSE)

        ash_out <- ashr:::estimate_mixprop(data = ash_data, g = ash_g, prior = lambda_seq,
                                           optmethod = optmethod, control = list())

        pi_new <- ash_out$optreturn$pihat

        ## stoping criteria
        llike_new <- uniform_mix_llike(pi_vals = pi_new, z2 = z_new, xi = xi_new,
                                   betahat_ols = betahat_ols, S_diag = S_diag,
                                   alpha_tilde = alpha_tilde, a_seq = a_seq,
                                   b_seq = b_seq, lambda_seq = lambda_seq,
                                   degrees_freedom = degrees_freedom)
        llike_vec <- c(llike_vec, llike_new)

        assertthat::assert_that((llike_new -llike_old) > - 10 ^ -14)


        err <- abs(llike_new / llike_old - 1)
        iter_index <- iter_index + 1

        if (plot_update) {
            graphics::plot(llike_vec, type = "l", xlab = "Iteration", ylab = "log-likelihood")
            graphics::mtext(paste0("Diff: ", format(err, digits = 2)))
        }
    }

    return(list(pi_vals = pi_new, z2 = z_new, xi = xi_new))
}


#' Gradient wrt z of \code{\link{uniform_mix_llike}}.
#'
#' @inheritParams uniform_mix_fix
#'
#' @author David Gerard
#'
#' @seealso \code{\link{uniform_mix_llike}} for the objective function
#'     and \code{\link{mouthwash_coordinate}} for where this function
#'     is used.
mouthwash_z_grad <- function(pi_vals, z2, xi, betahat_ols, S_diag, alpha_tilde, a_seq, b_seq,
                             lambda_seq, degrees_freedom){

    ## Make sure input is OK --------------------------------------------------
    M <- length(a_seq)
    k <- length(z2)
    p <- length(betahat_ols)
    assertthat::are_equal(length(b_seq), M)
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::are_equal(length(S_diag), p)
    assertthat::are_equal(nrow(alpha_tilde), p)
    assertthat::are_equal(ncol(alpha_tilde), k)
    assertthat::are_equal(length(pi_vals), M)
    zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    assertthat::are_equal(length(zero_spot), 1)

    az <- alpha_tilde %*% z2
    resid_vec <- c(betahat_ols - az)

    ## calculate ftildes (likelihoods) --------------------------------------
    sd_mat <- matrix(rep(sqrt(xi * S_diag), M), ncol = M)
    left_means  <- outer(resid_vec, a_seq, FUN = "-") / sd_mat
    right_means <- outer(resid_vec, b_seq, FUN = "-") / sd_mat

    denom_mat <- matrix(rep(b_seq - a_seq, times = p), byrow = TRUE, ncol = M)


    ftilde_mat <- ptdiff_mat(left_means, right_means, degrees_freedom = degrees_freedom) /
        denom_mat

    dense_vec <- dt_wrap(x = betahat_ols, df = degrees_freedom,
                        mean = az, sd = sqrt(xi * S_diag))
    ftilde_mat[, zero_spot] <- dense_vec

    ## calculate fbars (derivatives) ---------------------------------------
    tdiff_mat <- stats::dt(right_means, df = degrees_freedom) / sd_mat -
        stats::dt(left_means, df = degrees_freedom) / sd_mat

    fbar_mat <- tdiff_mat / denom_mat

    fbar_mat[, zero_spot] <- dense_vec * resid_vec * (degrees_freedom + 1) /
        (degrees_freedom * xi * S_diag + resid_vec ^ 2)

    ## calcualte weights for rows of alpha_tilde ---------------------------
    weight_vec <- rowSums(sweep(fbar_mat, MARGIN = 2, STATS = pi_vals, FUN = `*`)) /
        rowSums(sweep(ftilde_mat, MARGIN = 2, STATS = pi_vals, FUN = `*`))

    ## now gradient --------------------------------------------------------
    grad_final <- crossprod(alpha_tilde, weight_vec)

    return(grad_final)
}
