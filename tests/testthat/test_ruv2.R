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
    fa_func <- qmle_ruv2
    limmashrink <- FALSE
    include_intercept <- FALSE
    gls <- TRUE
    likelihood <- "normal"

    expect_warning(vout <- vruv2_old(Y = Y, X = X, ctl = ctl, k = num_sv,
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


test_that("pcaruv2_fix increases likelihood", {
    set.seed(966)
    n <- 11
    p <- 53
    r <- 5
    vr <- 2
    lambda <- 2

    A <- matrix(stats::rnorm(n * r), nrow = n)
    B <- matrix(stats::rnorm(r * p), ncol = p)
    sig_diag <- stats::rchisq(p, df = 4) / 4
    E <- matrix(stats::rnorm(n * p), nrow = n) %*% diag(sqrt(sig_diag))
    E[1:vr, ] <- E[1:vr, ] * lambda

    Y <- A %*% B + E

    ## make sure objective function increases
    resids <- E
    R1 <- resids[1:vr, ]
    R2 <- resids[(vr + 1): n, ]
    r1 <- colSums(R1 ^ 2)
    r2 <- colSums(R2 ^ 2)
    Sigma_init <- r2 / (n - vr)
    lambda_init <- mean(r1 / Sigma_init) / vr
    sig_lambda <- c(Sigma_init, lambda_init)

    itermax <- 20
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- pcaruv2_obj(sig_lambda = sig_lambda,
                                r1 = r1, r2 = r2, n = n,
                                vr = vr)

    for (index in 2:itermax) {
        sig_lambda <- pcaruv2_fix(sig_lambda = sig_lambda, r1 = r1, r2 = r2, n = n, vr = vr)
        llike_vec[index] <- pcaruv2_obj(sig_lambda = sig_lambda,
                                        r1 = r1, r2 = r2, n = n,
                                        vr = vr)
    }
    expect_true(all(llike_vec[1:(length(llike_vec) - 1)] -
                    llike_vec[2:length(llike_vec)] < 10 ^ -13))


}
)


test_that("pca_ruv2 works", {
    n <- 11
    p <- 1011
    r <- 5
    vr <- 2
    lambda <- 2

    A <- matrix(stats::rnorm(n * r), nrow = n)
    B <- matrix(stats::rnorm(r * p), ncol = p)
    sig_diag <- stats::rchisq(p, df = 4) / 4
    E <- matrix(stats::rnorm(n * p), nrow = n) %*% diag(sqrt(sig_diag))
    E[1:vr, ] <- E[1:vr, ] * sqrt(lambda)

    Y <- A %*% B + E

    pc1 <- pca_ruv2(Y = Y, r = r, vr = vr, mle = FALSE)
    pc2 <- pca_ruv2(Y = Y, r = r, vr = vr, mle = TRUE)
    pc3 <- pca_naive(Y = Y, r = r)

    ## plot(pc1$sig_diag, sig_diag)
    ## plot(pc2$sig_diag, sig_diag)
    ## abline(0, 1)
    ## plot(pc3$sig_diag, sig_diag)
    ## abline(0, 1)


    sum((pc3$sig_diag - sig_diag) ^ 2)
    sum((pc2$sig_diag - sig_diag) ^ 2)
    sum((pc1$sig_diag - sig_diag) ^ 2)

    pca_ruv2(E, r = 0, vr = vr)
}
)



test_that("update_sig_alpha increases qmle_obj", {
    set.seed(98)
    n <- 11
    p <- 41
    r <- 5
    vr <- 2
    lambda <- 2

    A <- matrix(stats::rnorm(n * r), nrow = n)
    B <- matrix(stats::rnorm(r * p), ncol = p)
    sig_diag <- stats::rchisq(p, df = 4) / 4
    E <- matrix(stats::rnorm(n * p), nrow = n) %*% diag(sqrt(sig_diag))
    E[1:vr, ] <- E[1:vr, ] * sqrt(lambda)

    Y <- A %*% B + E

    Y1 <- Y[1:vr, ]
    Y2 <- Y[(vr + 1):n, ]

    alpha <- B


    uout1 <- update_sig_alpha(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                              Y1 = Y1, Y2 = Y2)
    uout2 <- update_sig_alpha_basic(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                                    Y1 = Y1, Y2 = Y2)

    expect_equal(uout1, uout2)

    llike1 <- qmle_obj(alpha = alpha, sig_diag = sig_diag,
                       lambda = 2, Y1 = Y1, Y2 = Y2)
    llike2 <- qmle_obj_basic(alpha = alpha, sig_diag = sig_diag,
                       lambda = 2, Y1 = Y1, Y2 = Y2)

    expect_equal(llike1, llike2)

    itermax <- 10
    llike_vec1 <- rep(NA, itermax)
    llike_vec2 <- rep(NA, itermax)
    llike_vec1[1] <- qmle_obj(alpha = alpha, sig_diag = sig_diag,
                              lambda = 2, Y1 = Y1, Y2 = Y2)
    llike_vec2[1] <- qmle_obj_basic(alpha = alpha,
                                    sig_diag = sig_diag, lambda = 2,
                                    Y1 = Y1, Y2 = Y2)
    for (index in 2:itermax) {
        uout <- update_sig_alpha_basic(alpha = alpha, sig_diag = sig_diag,
                                 lambda = lambda, Y1 = Y1, Y2 = Y2)
        alpha <- uout$alpha
        sig_diag <- uout$sig_diag
        llike_vec1[index] <- qmle_obj(alpha = alpha, sig_diag = sig_diag,
                                     lambda = 2, Y1 = Y1, Y2 = Y2)
        llike_vec2[index] <- qmle_obj_basic(alpha = alpha, sig_diag = sig_diag,
                                     lambda = 2, Y1 = Y1, Y2 = Y2)
    }

    expect_true(all(llike_vec2[1:(itermax - 1)] - llike_vec2[2:itermax] < 10 ^ -14))
    expect_true(all(llike_vec1[1:(itermax - 1)] - llike_vec1[2:itermax] < 10 ^ -14))
    expect_equal(llike_vec1, llike_vec2)

    qout <- qmle_ruv2(Y = Y, r = r, vr = vr)

    Zclassic <- Y %*% diag(1 / qout$sig_diag) %*% t(qout$alpha) %*%
        solve(qout$alpha %*% diag(1 / qout$sig_diag) %*% t(qout$alpha))

    expect_equal(Zclassic, qout$Z)
}
)
