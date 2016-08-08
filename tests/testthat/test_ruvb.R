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


test_that("bsvd works OK", {
    set.seed(881)
    n <- 11
    p <- 37
    ncontrols <- 13
    k <- 5
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

    bsvd_out <- bsvd(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k,
                     nsamp = 10, keep = 1, print_update = FALSE)

    expect_equal(dim(bsvd_out$Y22_array)[1:2], c(ncovs, p - ncontrols))
}
)
