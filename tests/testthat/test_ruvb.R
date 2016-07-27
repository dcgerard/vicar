library(vicar)
context("ruvb")

test_that("See if RUV3 returns correct OLS estimates and matches vruv2", {
    set.seed(730)
    n <- 17
    p <- 103
    k <- 5
    q <- 3

    X <- matrix(rnorm(n * q), nrow = n)
    beta <- matrix(rnorm(q * p), nrow = q)
    Z <- matrix(rnorm(n * k), nrow = n)
    alpha <- matrix(rnorm(k * p), nrow = k)
    E <- matrix(rnorm(n * p), nrow = n)
    Y <- X %*% beta + Z %*% alpha + E
    ctl <- rep(FALSE, length = p)
    ctl[1:13] <- TRUE
    include_intercept <- FALSE
    cov_of_interest <- 2:3
    fa_args <- list()



}
)
