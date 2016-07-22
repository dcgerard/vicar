library(vicar)
context("ruv3")

test_that("See if RUV3 returns correct OLS estimates and matches vruv2", {
    set.seed(21)
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
    limmashrink = FALSE
    fa_func <- pca_naive
    fa_args <- list()
    likelihood <- "normal"
    gls <- FALSE
    cov_of_interest <- 2:3

    ruv3out <- ruv3(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest,
                    include_intercept = include_intercept,
                    limmashrink = limmashrink, gls = TRUE)

    lmout <- limma::lmFit(t(Y), X)

    expect_equal(c(t(lmout$coefficients)), c(ruv3out$betahat_ols))
    expect_equal(c(lmout$stdev.unscaled * lmout$sigma), c(t(ruv3out$sebetahat_ols)))

    ruv2out <- vruv2(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest,
                     include_intercept = include_intercept, limmashrink = limmashrink,
                     gls = TRUE, likelihood = "normal", use_factor = FALSE,
                     fa_limmashrink = FALSE)


    expect_equal(ruv3out$alphahat, ruv2out$alphahat)
    expect_equal(ruv2out$debuglist$Z3, ruv3out$debuglist$Z3)
    expect_equal(ruv2out$debuglist$Z2, ruv3out$debuglist$Z2)
    expect_equal(ruv3out$betahat[, !ctl], ruv2out$betahat[, !ctl])
    expect_equal(ruv2out$sebetahat[, !ctl], ruv2out$sebetahat[, !ctl])

    ## plot(ruv3out$betahat, ruv2out$betahat)
    ## abline(0, 1)
}
)


test_that("ruvimpute works ok", {
    set.seed(24)
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
    impute_func <- softimpute_wrapper
    impute_args <- list()
    do_variance <- TRUE

    sout <- ruvimpute(Y = Y, X = X, ctl = ctl,
                      impute_func = impute_func,
                      cov_of_interest = cov_of_interest)
    sout2 <- ruvimpute(Y = Y, X = X, ctl = ctl,
                       impute_func = flashr_wrapper,
                       impute_args = list(max_rank = n - k - q - 1),
                       cov_of_interest = cov_of_interest)

    ruv3out <- ruv3(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest)
    ruv4out <- vruv4(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest)
    ruv2out <- vruv2(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest)
    plot(sout$beta2hat, sout2$beta2hat)

    mean((beta[cov_of_interest, !ctl] - ruv3out$betahat[, !ctl]) ^ 2)
    mean((beta[cov_of_interest, !ctl] - ruv2out$betahat[, !ctl]) ^ 2)
    mean((beta[cov_of_interest, !ctl] - t(ruv4out$betahat[!ctl, ])) ^ 2)
    mean((beta[cov_of_interest, !ctl] - sout2$beta2hat) ^ 2)
    mean((beta[cov_of_interest, !ctl] - sout$beta2hat) ^ 2)

}
)
