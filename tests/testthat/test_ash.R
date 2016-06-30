library(vicar)
context("ASH Wrappers")

test_that("ash_ruv4 works", {
    set.seed(21)
    n <- 11
    p <- 113
    k <- 2
    q <- 3

    pi_vals <- c(0.5, 0.3, 0.2)
    sd_seq  <- c(0, 1, 2)

    X <- cbind(rep(1, n), sample(c(0, 1), size = n, replace = TRUE))
    beta <- matrix(NA, nrow = k, ncol = p)
    beta[1, ] <- stats::rnorm(p)
    beta[2, ] <- rmixnorm(n = p, pi_vals = pi_vals, sd_seq = sd_seq)
    alpha <- matrix(stats::rnorm(q * p), ncol = p)
    Z <- matrix(stats::rnorm(n * q), nrow = n)
    sig_diag <- stats::rchisq(p, 5) / 5
    E <- matrix(rnorm(n * p), nrow = n) %*% diag(sqrt(sig_diag))
    Y <- X %*% beta + Z %*% alpha + E

    which_null <- beta[2, ] == 0
    ctl <- which_null
    ncontrol <- 31
    ctl[ctl][sample(1:sum(ctl), size = sum(ctl) - ncontrol)] <- FALSE

    vout <- vruv4(Y = Y, X = X, ctl = ctl, k = q, likelihood = "normal")
    vashout <- ash_ruv4(Y = Y, X = X, ctl = ctl, k = q, likelihood = "normal")

    expect_equal(vout, vashout$ruv4)

}
)
