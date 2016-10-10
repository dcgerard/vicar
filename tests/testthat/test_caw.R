context("CAW")

test_that("caw works", {
    n <- 11
    p <- 37
    q <- 2
    k <- 2

    ## data
    X     <- matrix(rnorm(n * q), nrow = n)
    beta  <- matrix(rnorm(q * p), nrow = q)
    Z     <- matrix(rnorm(n * k), nrow = n)
    alpha <- matrix(rnorm(k * p), nrow = k)
    E     <- matrix(rnorm(n * p), nrow = n)
    Y     <- X %*% beta + Z %*% alpha + E

    ## params
    cov_of_interest   <- 1
    limmashrink       <- TRUE
    weight_func       <- ash_wrap
    weight_args       <- list()
    fa_func           <- pca_naive
    fa_args           <- list()
    scale_var         <- TRUE
    include_intercept <- FALSE


}
)
