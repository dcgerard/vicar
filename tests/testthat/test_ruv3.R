library(vicar)
context("ruv3")

test_that("See if RUV3 returns correct OLS estimates", {
    set.seed(19)
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
                    limmashrink = limmashrink, gls = FALSE)

    lmout <- limma::lmFit(t(Y), X)

    expect_equal(c(t(lmout$coefficients)), c(ruv3out$betahat_ols))
    expect_equal(c(lmout$stdev.unscaled * lmout$sigma), c(t(ruv3out$sebetahat_ols)))


    ## ruv2out <- ruv::RUV2(Y = Y, X = X[, cov_of_interest, drop = FALSE], ctl = ctl,
    ##                      k = q, Z = X[, -cov_of_interest, drop = FALSE])
    ## ruv4out <- ruv::RUV4(Y = Y, X = X[, cov_of_interest, drop = FALSE], ctl = ctl,
    ##                      k = q, Z = X[, -cov_of_interest, drop = FALSE])
    ## vruv4out <- vruv4(Y = Y, X = X, ctl = ctl, k = q, cov_of_interest = cov_of_interest,
    ##                   likelihood = likelihood, include_intercept = include_intercept,
    ##                   limmashrink = limmashrink, gls = TRUE)

    ## plot(ruv3out$betahat[, !ctl], ruv2out$betahat[, !ctl])
    ## plot(ruv3out$betahat[, !ctl], ruv4out$betahat[, !ctl])
    ## plot(ruv3out$betahat[, !ctl], t(vruv4out$betahat[!ctl, ]))
    ## plot(ruv2out$betahat[, !ctl], ruv4out$betahat[, !ctl])

    ## mean((ruv3out$betahat[, !ctl] - beta[cov_of_interest, !ctl]) ^ 2)
    ## mean((ruv2out$betahat[, !ctl] - beta[cov_of_interest, !ctl]) ^ 2)
    ## mean((ruv4out$betahat[, !ctl] - beta[cov_of_interest, !ctl]) ^ 2)
    ## mean((t(vruv4out$betahat[!ctl, ]) - beta[cov_of_interest, !ctl]) ^ 2)

}
)
