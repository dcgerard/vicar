context("MOUTHWASH functions")


test_that("mouthwash works ok", {
    set.seed(68)
    n <- 11
    p <- 19
    k <- 3
    cov_of_interest <- 2
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    include_intercept <- FALSE
    limmashrink       <- TRUE
    fa_func           <- pca_naive
    fa_args           <- list()
    likelihood        <- "t"
    mixing_dist       <- "uniform"
    degrees_freedom   <- NULL
    pi_init           <- NULL
    tau_seq           <- NULL
    lambda_seq        <- NULL
    lambda0           <- 10

}
)

test_that("normal_mix_llike and normal_mix_fix works ok", {
    set.seed(124895)
    p <- 103
    k <- 3
    S_diag <- stats::rchisq(p, 5)
    alpha_tilde <- matrix(stats::rnorm(k * p), nrow = p)
    z <- matrix(stats::rnorm(k), ncol = 1)
    beta <- matrix(stats::rnorm(p), ncol = 1)
    betahat_ols <- beta + alpha_tilde %*% z + rnorm(p, mean = 0, sd = sqrt(S_diag))

    M             <- 23
    grid_seq      <- seq(0, 10, length = M)
    lambda_seq    <- rep(1, M)
    lambda_seq[1] <- 10
    pi_vals <- rep(1 / M, length = M)
    xi <- 1

    fval <- normal_mix_llike(betahat_ols = betahat_ols, alpha_tilde = alpha_tilde, z2 = z, xi = xi,
                             S_diag = S_diag, tau2_seq = grid_seq, lambda_seq = lambda_seq,
                             pi_vals = pi_vals)

    mix_var  <- outer(xi * S_diag, grid_seq, FUN = `+`) ## p by M
    mix_mean <- matrix(rep(alpha_tilde %*% z, M), ncol = M)
    mix_obs  <- matrix(rep(betahat_ols, M), ncol = M)
    eval <- sum(log(rowSums(stats::dnorm(x = mix_obs, mean = mix_mean, sd = sqrt(mix_var),
                                         log = FALSE) %*% diag(pi_vals))))
    pen <- -log(M) * 9

    expect_equal(fval, eval + pen)


    z2 <- z
    itermax <- 20
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- normal_mix_llike(pi_vals = pi_vals, z2 = z2, xi = xi,
                                     betahat_ols = betahat_ols,
                                     S_diag = S_diag, alpha_tilde = alpha_tilde,
                                     tau2_seq = grid_seq,
                                     lambda_seq = lambda_seq)

    for (index in 2:itermax) {

        fout <- normal_mix_fix(pi_vals = pi_vals, z2 = z2, xi = xi, betahat_ols = betahat_ols,
                               S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = grid_seq,
                               lambda_seq = lambda_seq, scale_var = TRUE)
        pi_vals <- fout$pi_vals
        z2 <- fout$z2
        xi <- fout$xi
        llike_vec[index] <- normal_mix_llike(pi_vals = pi_vals, z2 = z2, xi = xi,
                                             betahat_ols = betahat_ols,
                                             S_diag = S_diag, alpha_tilde = alpha_tilde,
                                             tau2_seq = grid_seq,
                                             lambda_seq = lambda_seq)
    }

    expect_true(all(llike_vec[1:(itermax - 1)] <= llike_vec[2:itermax]))


    pi_vals <- rep(1 / M, length = M)
    z2 <- z
    xi <- 1
    pizxi_vec <- c(pi_vals, z2, xi)
    llike_vec2 <- rep(NA, length = itermax)
    llike_vec2[1] <- normal_mix_llike_wrapper(pizxi_vec = pizxi_vec,
                                              betahat_ols = betahat_ols,
                                              S_diag = S_diag,
                                              alpha_tilde = alpha_tilde,
                                              tau2_seq = grid_seq,
                                              lambda_seq = lambda_seq,
                                              scale_var = TRUE)
    for (index in 2:itermax) {
        pizxi_vec <- normal_mix_fix_wrapper(pizxi_vec = pizxi_vec,
                                            betahat_ols = betahat_ols,
                                            S_diag = S_diag,
                                            alpha_tilde = alpha_tilde,
                                            tau2_seq = grid_seq,
                                            lambda_seq = lambda_seq,
                                            scale_var = TRUE)
        llike_vec2[index] <- normal_mix_llike_wrapper(pizxi_vec = pizxi_vec,
                                                      betahat_ols = betahat_ols,
                                                      S_diag = S_diag,
                                                      alpha_tilde = alpha_tilde,
                                                      tau2_seq = grid_seq,
                                                      lambda_seq = lambda_seq,
                                                      scale_var = TRUE)
    }

    expect_equal(llike_vec, -1 * llike_vec2)
}
)
