context("MOUTHWASH functions")


test_that("mouthwash works ok", {
    set.seed(68)
    n <- 20
    p <- 102
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

    betahat_ols <- (solve(t(X) %*% X) %*% t(X) %*% Y)[3, ]

    start_unif <- proc.time()
    mout <- mouthwash(Y = Y, X = X, k = 1, degrees_freedom = 2,
                      likelihood = "t", mixing_dist = "uniform",
                      lambda0 = 10)

    ## mold <- readRDS("old_mouth.Rds")
    ## expect_equal(mout, mold)

    tot_unif <- proc.time() - start_unif

    start_normal <- proc.time()
    mout <- mouthwash(Y = Y, X = X, k = 1, mixing_dist = "normal", likelihood = "normal")
    tot_normal <- proc.time() - start_normal
    mout$pi0

    expect_error(mouthwash(Y = Y, X = X, k = 1, degrees_freedom = 2,
                           likelihood = "t", mixing_dist = "normal"))

    mout2 <- mouthwash(Y = Y, X = X, k = 1, mixing_dist = "normal",
                      likelihood = "normal", sprop = 0.5)

    mout3 <- mouthwash(Y = Y, X = X, k = 1, mixing_dist = "normal",
                       likelihood = "normal", sprop = 0.1)

    expect_equal(mout$result$betahat, mout2$result$betahat, tol = 10^-5)
    expect_equal(mout$result$betahat, mout3$result$betahat, tol = 10^-5)


    mse1 <- sum((mout$result$PosteriorMean - mout2$result$PosteriorMean) ^ 2)
    mse2 <- sum((mout$result$PosteriorMean - mout3$result$PosteriorMean) ^ 2)
    expect_true(mse1 > mse2)
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


test_that("uniform_mix_llike and uniform_mix_fix work", {
    set.seed(991)

    p <- 103
    k <- 3
    S_diag <- stats::rchisq(p, 5)
    alpha_tilde <- matrix(stats::rnorm(k * p), nrow = p)
    z2 <- matrix(stats::rnorm(k), ncol = 1)
    beta <- matrix(stats::rnorm(p), ncol = 1)
    betahat_ols <- beta + alpha_tilde %*% z2 + rnorm(p, mean = 0, sd = sqrt(S_diag))

    M             <- 23
    a_seq         <- seq(-10, 0, length = M)
    b_seq         <- seq(10, 0, length = M)
    lambda_seq    <- rep(1, M)
    lambda_seq[length(lambda_seq)] <- 10
    pi_vals <- rep(1 / M, length = M)
    xi <- 1
    degrees_freedom <- 3
    scale_var <- TRUE
    pi_init_type <- "zero_conc"

    itermax <- 20
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- uniform_mix_llike(pi_vals = pi_vals, z2 = z2, xi = xi,
                                      betahat_ols = betahat_ols, S_diag = S_diag,
                                      alpha_tilde = alpha_tilde, a_seq = a_seq, b_seq = b_seq,
                                      lambda_seq = lambda_seq, degrees_freedom = degrees_freedom)
    xi_new <- xi
    z_new  <- z2
    pi_new <- pi_vals
    for (iter_index in 2:itermax) {
        uout <- uniform_mix_fix(pi_vals = pi_new, z2 = z_new, xi = xi_new,
                                betahat_ols = betahat_ols, S_diag = S_diag,
                                alpha_tilde = alpha_tilde, a_seq = a_seq, b_seq = b_seq,
                                lambda_seq = lambda_seq, degrees_freedom = degrees_freedom,
                                scale_var = TRUE)
        xi_new <- uout$xi
        z_new  <- uout$z2
        pi_new <- uout$pi_vals
        llike_vec[iter_index] <- uniform_mix_llike(pi_vals = pi_new, z2 = z_new, xi = xi_new,
                                                   betahat_ols = betahat_ols, S_diag = S_diag,
                                                   alpha_tilde = alpha_tilde, a_seq = a_seq,
                                                   b_seq = b_seq, lambda_seq = lambda_seq,
                                                   degrees_freedom = degrees_freedom)
    }

    expect_true(all(llike_vec[1:(itermax - 1)] <= llike_vec[2:itermax]))


    pizxi_vec <- c(pi_vals, z2, xi)
    llike_vec2 <- rep(NA, length = itermax)
    llike_vec2[1] <- uniform_mix_llike_wrapper(pizxi_vec = pizxi_vec, betahat_ols = betahat_ols,
                                               S_diag = S_diag, alpha_tilde = alpha_tilde,
                                               a_seq = a_seq, b_seq = b_seq,
                                               lambda_seq = lambda_seq,
                                               degrees_freedom = degrees_freedom)
    for (iter_index in 2:itermax) {
        pizxi_vec <- uniform_mix_fix_wrapper(pizxi_vec = pizxi_vec, betahat_ols = betahat_ols,
                                             S_diag = S_diag, alpha_tilde = alpha_tilde,
                                             a_seq = a_seq, b_seq = b_seq,
                                             lambda_seq = lambda_seq,
                                             degrees_freedom = degrees_freedom,
                                             scale_var = TRUE)
        llike_vec2[iter_index] <- uniform_mix_llike_wrapper(pizxi_vec = pizxi_vec,
                                                            betahat_ols = betahat_ols,
                                                            S_diag = S_diag,
                                                            alpha_tilde = alpha_tilde,
                                                            a_seq = a_seq, b_seq = b_seq,
                                                            lambda_seq = lambda_seq,
                                                            degrees_freedom = degrees_freedom,
                                                            scale_var = TRUE)
    }

    expect_equal(-1 * llike_vec2, llike_vec)

}
)


test_that("ptdiff_mat and ptdiff_mat_log work ok", {
    set.seed(525)
    n <- 10
    p <- 1000
    X <- matrix(stats::rchisq(n * p, df = 1), nrow = n)
    Y <- matrix(-1 * stats::rchisq(n * p, df = 1), nrow = n)

    pt1 <- log(ptdiff_mat(X, Y, degrees_freedom = 3))
    pt2 <- ptdiff_mat_log(X, Y, degrees_freedom = 3)
    expect_equal(pt1, pt2)

}
)



test_that("mouthwash_coord works", {
    set.seed(724)

    p <- 103
    k <- 3
    S_diag <- stats::rchisq(p, 5)
    alpha_tilde <- matrix(stats::rnorm(k * p), nrow = p)
    z2 <- matrix(stats::rnorm(k), ncol = 1)
    beta <- matrix(stats::rnorm(p), ncol = 1)
    betahat_ols <- beta + alpha_tilde %*% z2 + rnorm(p, mean = 0, sd = sqrt(S_diag))

    M             <- 23
    a_seq         <- seq(-10, 0, length = M)
    b_seq         <- seq(10, 0, length = M)
    lambda_seq    <- rep(1, M)
    lambda_seq[length(lambda_seq)] <- 10
    pi_init <- rep(1 / M, length = M)
    xi_init <- 1
    degrees_freedom <- 3
    scale_var <- TRUE
    tol <- 10 ^ -6
    maxit <- 100
    z_init <- z2

    mout <- mouthwash_coordinate(pi_init = pi_init, z_init = z_init, xi_init = xi_init,
                                 betahat_ols = betahat_ols,
                                 S_diag = S_diag, alpha_tilde = alpha_tilde, a_seq = a_seq,
                                 b_seq = b_seq, lambda_seq = lambda_seq,
                                 degrees_freedom = degrees_freedom, scale_var = scale_var)
}
)


test_that("mouthwash_z_grad returns similar results for df large and df = Inf", {
  set.seed(916)

  p <- 103
  k <- 3
  S_diag <- stats::rchisq(p, 5)
  alpha_tilde <- matrix(stats::rnorm(k * p), nrow = p)
  z2 <- matrix(stats::rnorm(k), ncol = 1)
  beta <- matrix(stats::rnorm(p), ncol = 1)
  betahat_ols <- beta + alpha_tilde %*% z2 + rnorm(p, mean = 0, sd = sqrt(S_diag))

  M             <- 23
  a_seq         <- seq(-10, 0, length = M)
  b_seq         <- seq(10, 0, length = M)
  lambda_seq    <- rep(1, M)
  lambda_seq[length(lambda_seq)] <- 10
  pi_vals <- rep(1 / M, length = M)
  xi <- 1
  degrees_freedom <- 3
  scale_var <- TRUE

  mout_inf <- mouthwash_z_grad(pi_vals = pi_vals, z2 = z2, xi = xi,
                               betahat_ols = betahat_ols, S_diag = S_diag,
                               alpha_tilde = alpha_tilde, a_seq = a_seq, b_seq = b_seq,
                               lambda_seq = lambda_seq, degrees_freedom = Inf)

  mout_10000 <- mouthwash_z_grad(pi_vals = pi_vals, z2 = z2, xi = xi,
                                 betahat_ols = betahat_ols, S_diag = S_diag,
                                 alpha_tilde = alpha_tilde, a_seq = a_seq, b_seq = b_seq,
                                 lambda_seq = lambda_seq, degrees_freedom = 10000)

  expect_equal(mout_inf, mout_10000, tolerance = 10 ^ -4)
}
)
