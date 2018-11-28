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


test_that("ruvimpute returns RUV2, RUV3, and RUV4", {
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
    cov_of_interest <- 2

    imp2 <- ruvimpute(Y = Y, X = X, ctl = ctl, k = k,
                      impute_func = impute_ruv_reproduce,
                      impute_args = list(impute_type = "ruv2"),
                      cov_of_interest = cov_of_interest,
                      include_intercept = FALSE)
    imp3 <- ruvimpute(Y = Y, X = X, ctl = ctl, k = k,
                      impute_func = impute_ruv_reproduce,
                      impute_args = list(impute_type = "ruv3"),
                      cov_of_interest = cov_of_interest,
                      include_intercept = FALSE)
    imp4 <- ruvimpute(Y = Y, X = X, ctl = ctl, k = k,
                      impute_func = impute_ruv_reproduce,
                      impute_args = list(impute_type = "ruv4"),
                      cov_of_interest = cov_of_interest,
                      include_intercept = FALSE)

    impout_miss <- ruvimpute(Y = Y, X = X, ctl = ctl, k = k,
                      impute_func = em_miss,
                      cov_of_interest = cov_of_interest,
                      include_intercept = FALSE)

    ruv3out <- ruv3(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest,
                    include_intercept = FALSE, gls = FALSE)
    ruv4out <- ruv::RUV4(Y = Y, X = X[, cov_of_interest, drop = FALSE], ctl = ctl,
                         k = k, Z = X[, -cov_of_interest, drop = FALSE])
    ruv2out <- ruv::RUV2(Y = Y, X = X[, cov_of_interest, drop = FALSE], ctl = ctl,
                         k = k, Z = X[, -cov_of_interest, drop = FALSE])

    my_ruv4out <- vruv4(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest,
                        likelihood = "normal", limmashrink = FALSE, include_intercept = FALSE,
                        gls = FALSE)

    if (packageVersion("ruv") == "0.9.6") {
      expect_equal(c(imp2$beta2hat), c(ruv2out$betahat[, !ctl]))
      expect_equal(c(imp3$beta2hat), c(ruv3out$betahat[, !ctl]))
      expect_equal(c(imp4$beta2hat), c(ruv4out$betahat[, !ctl]))
      expect_equal(c(t(my_ruv4out$betahat[!ctl, ])), c(ruv4out$betahat[, !ctl]))
    }


    ## plot(sout$beta2hat, sout2$beta2hat)

    # mean((beta[cov_of_interest, !ctl] - ruv3out$betahat[, !ctl]) ^ 2)
    # mean((beta[cov_of_interest, !ctl] - ruv2out$betahat[, !ctl]) ^ 2)
    # mean((beta[cov_of_interest, !ctl] - t(ruv4out$betahat[, !ctl ])) ^ 2)
    # mean((beta[cov_of_interest, !ctl] - impout_miss$beta2hat) ^ 2)
    ## mean((beta[cov_of_interest, !ctl] - sout2$beta2hat) ^ 2)


    ruv3out <- ruv3(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = 1:ncol(X),
                    include_intercept = FALSE, gls = FALSE)


}
)


test_that("Using em_miss in RUVimpute gets variance estimates", {
    set.seed(998)
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
    cov_of_interest <- 2:3

    emout <- ruvem(Y = Y, X = X, ctl = ctl, k = k, cov_of_interest = cov_of_interest)

}
)
