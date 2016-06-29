library(vicar)
context("RUV2")


test_that("vruv2 works", {
    set.seed(944)
    n <- 13
    p <- 73
    k <- 5
    cov_of_interest <- c(1, 2)
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:41] <- 0
    ctl <- rep(FALSE, length = p)
    ctl[1:23] <- TRUE
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 1
    fa_args <- list()
    fa_func <- pca_naive
    limmashrink <- FALSE
    include_intercept <- FALSE
    gls <- TRUE
    likelihood <- "normal"

    expect_warning(vout <- vruv2(Y = Y, X = X, ctl = ctl, k = num_sv,
                                 cov_of_interest = cov_of_interest,
                                 likelihood = likelihood,
                                 include_intercept = FALSE,
                                 limmashrink = FALSE))

    expect_equal(vout$betahat_ols[, !ctl], vout$betahat[, !ctl])
    expect_equal(vout$sigma2_unadjusted[!ctl], vout$sigma2_adjusted[!ctl])


    expect_equal(solve(t(cbind(X, vout$Zhat)) %*%
                       cbind(X, vout$Zhat))[cov_of_interest, cov_of_interest],
                 vout$mult_mat)

    if (requireNamespace("ruv", quietly = TRUE)) {
        ruv2out <- ruv::RUV2(Y = Y, X = X[, cov_of_interest, drop = FALSE],
                             Z = X[, -cov_of_interest, drop = FALSE],
                             ctl = ctl, k = num_sv)
        expect_equal(ruv2out$betahat, vout$betahat_ols)
        expect_equal(ruv2out$sigma2, vout$sigma2_unadjusted)
        expect_equal(1, length(unique(round(c(ruv2out$alpha[, !ctl]) /
                                            c(vout$alpha[, !ctl]), digits = 3))))
    }
}
)


test_that("limmashrink works OK", {
    set.seed(11)
    n <- 13
    p <- 19
    k <- 5
    cov_of_interest <- c(1, 2)
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest[1],] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 1
    fa_args <- list()
    fa_func <- pca_naive
    limmashrink <- FALSE
    include_intercept <- FALSE
    gls <- TRUE
    likelihood <- "normal"

    vout <- vruv2(Y = Y, X = X, ctl = ctl, k = num_sv,
                  cov_of_interest = cov_of_interest,
                  likelihood = likelihood,
                  include_intercept = FALSE,
                  limmashrink = TRUE)
}
)


test_that("outlist.Rd data works", {
    load("outlist.Rd")
    vruv2(Y = outlist$Y, X = outlist$X, ctl = outlist$ctl, k = outlist$num_sv,
          likelihood = "normal")
}
)
