library(vicar)
context("QR Rotation and RUV4")

test_that("rotated_model same as ols when no confounders", {
    set.seed(68)
    n <- 11
    p <- 19
    k <- 3
    cov_of_interest <- 2
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    qr_x <- qr_ident(X)
    expect_true(all(qr_x$Q %*% qr_x$R - X < 10 ^ -15))

    num_sv <- 0

    rotate_out <- rotate_model(Y = Y, X = X, k = num_sv,
                               cov_of_interest = cov_of_interest,
                               include_intercept = FALSE,
                               limmashrink = FALSE, do_ols = FALSE)


    xinv <- solve(t(X) %*% X)
    expect_equal(1 / rotate_out$Rsub[1, 1] ^ 2, xinv[cov_of_interest, cov_of_interest])

    betahatmat <- xinv %*% t(X) %*% Y
    expect_equal(betahatmat[cov_of_interest, ], c(rotate_out$betahat_ols))

    resid_mat <- Y - X %*% betahatmat
    sig2est <- colSums(resid_mat ^ 2) / (n - k - num_sv)

    expect_equal(rotate_out$sig_diag, sig2est)

}
)


## should remove this before publishing
test_that("vicarius_ruv4 same as ashr::ash_ruv", {
    set.seed(68)
    n <- 11
    p <- 19
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 2

    ruv4_out <- vicarius_ruv4(Y = Y, X = X, ctl = ctl, k = num_sv,
                              cov_of_interest = cov_of_interest,
                              likelihood = "normal")

    ash_out <- ashr::ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv,
                             cov_of_interest = cov_of_interest,
                             likelihood = "normal", posthoc_inflate = FALSE)

    expect_equal(ash_out$ruv$multiplier, ruv4_out$multiplier)
    expect_equal(ash_out$ruv$alphahat, ruv4_out$alphahat * -1)
    expect_equal(ash_out$ruv$sebetahat_ols, ruv4_out$sebetahat_ols)
    expect_equal(ash_out$ruv$Z1, ruv4_out$Z1 * -1)
}
)


test_that("tregress_em increases likelihood", {
    set.seed(871)
    p  <- 21
    k  <- 5
    nu <- 5

    alpha <- matrix(stats::rnorm(p * k), nrow = p)
    Z     <- matrix(stats::rnorm(k), ncol = 1)
    sig_diag <- stats::rchisq(p, df = 2)
    E <- matrix(stats::rt(p, df = nu), ncol = 1) * sqrt(sig_diag)

    Y <- alpha %*% Z + E

    lambda_init <- 1
    Z_init <- rep(0, length = ncol(alpha))
    zlambda <- c(Z_init, lambda_init)

    itermax <- 20
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha,
                                 sig_diag = sig_diag, nu = nu)

    for (index in 2:itermax) {
        zlambda <- tregress_fix(zlambda = zlambda, Y = Y, alpha = alpha,
                                sig_diag = sig_diag, nu = nu)
        llike_vec[index] <- tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha,
                                         sig_diag = sig_diag, nu = nu)
    }
    expect_true(all(llike_vec[2:length(llike_vec)] >= llike_vec[1:(length(llike_vec) - 1)]))

    oout <- stats::optim(par = c(Z_init, lambda_init), fn = tregress_obj, Y = Y,
                         alpha = alpha, sig_diag = sig_diag, nu = nu,
                         control = list(fnscale = -1, maxit = 5000))

    expect_equal(llike_vec[index], oout$value, tol = 10 ^ -5)

    tregress_obj(zlambda = oout$par, Y = Y, alpha = alpha, sig_diag = sig_diag, nu = nu)
    tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha, sig_diag = sig_diag, nu = nu)

    expect_equal(zlambda, oout$par, tol = 10 ^ -3)

}
)


test_that("cruv4_multicov is same as RUV4 when no gls", {
    set.seed(71)
    n <- 11
    p <- 19
    k <- 4
    cov_of_interest <- c(1, 2)
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest[1],] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 2
    fa_args <- list()
    fa_func <- pca_naive
    limmashrink <- FALSE
    include_intercept <- FALSE
    gls <- TRUE
    likelihood <- "normal"

    vout <- vicarius_ruv4(Y = Y, X = X, ctl = ctl, k = num_sv, include_intercept = FALSE,
                          cov_of_interest = cov_of_interest, likelihood = "normal",
                          gls = FALSE)

    ruvout <- ruv::RUV4(Y = Y, X = X[, cov_of_interest, drop = FALSE],
                        ctl = ctl, k = num_sv,
                        Z = X[, -cov_of_interest, drop = FALSE])

    expect_equal(ruvout$betahat, t(vout$betahat))
    expect_equal(ruvout$sigma2, vout$sigma2)

}
)


test_that("cruv4 and cruv4_multicov give same answers", {
    load("testlist.Rd")
    alpha_scaled    <- testlist$rotate_out$alph / testlist$rotate_out$Rsub[1, 1]
    sig_diag_scaled <- testlist$rotate_out$sig_diag / testlist$rotate_out$Rsub[1, 1] ^ 2
    k               <- testlist$rotate_out$k
    uniout <- cruv4(betahat_ols = t(testlist$betahat_ols), alpha_scaled = alpha_scaled,
                    sig_diag_scaled = sig_diag_scaled, ctl = testlist$ctl,
                    degrees_freedom = testlist$degrees_freedom,
                    gls = testlist$gls, likelihood = "t")
    multiout <- cruv4_multicov(Y2 = t(testlist$rotate_out$Y2), alpha = testlist$rotate_out$alpha,
                               sig_diag = testlist$rotate_out$sig_diag, ctl = testlist$ctl,
                               Rsub = testlist$rotate_out$Rsub,
                               degrees_freedom = testlist$degrees_freedom,
                               gls = testlist$gls, likelihood = "t")

    expect_equal(multiout$betahat, uniout$betahat)
    expect_equal(multiout$multiplier, uniout$multiplier)
    expect_equal(c(multiout$sebetahat), uniout$sebetahat)
    expect_equal(uniout$Z1, multiout$Z1)
    expect_equal(uniout$betahat_ols, multiout$betahat_ols)
    expect_equal(uniout$sebetahat_ols, c(multiout$sebetahat_ols))
    expect_equal(uniout$tstats, multiout$tstats)
    expect_equal(uniout$pvalues, multiout$pvalues)
}
)
