context("EM with missing block")

test_that("em_miss_fix agrees with rubin and thayer when no covariates", {
    set.seed(881)
    n <- 11
    p <- 37
    ncontrols <- 13
    k <- 5
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    alpha_init <- t(pcout$alpha)
    sig_diag_init <- pcout$sig_diag
    alpha_sigma <- c(c(alpha_init), sig_diag_init)
    S2 <- crossprod(Y21)
    S  <- crossprod(cbind(Y31, Y32))



    alpha_sigma_new <- em_miss_fix(alpha_sigma = alpha_sigma, S = S, S2 = NULL, ncovs = 0,
                                       ncontrols = ncontrols, n = n, p = p, k = k)

    rubin_new <- rubin_copy(alpha_sigma = alpha_sigma, S = S, n = n, p = p, k = k)

    expect_equal(rubin_new, alpha_sigma_new)

    alphavec_me    <- alpha_sigma_new[1:(k * p)]
    alphavec_rubin <- rubin_new[1:(k * p)]
    sig_diag_me    <- alpha_sigma_new[((k * p) + 1):length(alpha_sigma_new)]
    sig_diag_rubin <- rubin_new[((k * p) + 1):length(rubin_new)]

    ## library(ggplot2)
    ## qplot(alphavec_me, alphavec_rubin) + geom_abline(slope = 1, intercept = 0)
    ## qplot(sig_diag_me, sig_diag_rubin) + geom_abline(slope = 1, intercept = 0)


    itermax <- 100
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- em_miss_obj(alpha_sigma = alpha_sigma, S = S, S2 = S2, ncovs = ncovs,
                                ncontrols = ncontrols, n = n, p = p, k = k)

    for(index in 2:itermax) {
        alpha_sigma_new <- em_miss_fix(alpha_sigma = alpha_sigma, S = S, S2 = NULL, ncovs = 0,
                                       ncontrols = ncontrols, n = n, p = p, k = k)
        alpha_sigma <- alpha_sigma_new

        llike_vec[index] <- em_miss_obj(alpha_sigma = alpha_sigma, S = S, S2 = NULL, ncovs = 0,
                                        ncontrols = ncontrols, n = n, p = p, k = k)
    }
    ## plot(llike_vec, type = "l")

    expect_true(all(llike_vec[1:(length(llike_vec) - 1)] <= llike_vec[2:length(llike_vec)]))

    ## alphavec_me <- alpha_sigma[1:(k * p)]
    ## sig_diag_me    <- alpha_sigma[((k * p) + 1):length(alpha_sigma)]
    ## cateout <- cate::fa.em(Y = cbind(Y31, Y32), r = k, maxiter = 1000, tol = 10 ^ -12)
    ## qplot(alphavec_me, c(t(cateout$Gamma))) + geom_abline(slope = 1, intercept = 0)
    ## qplot(sig_diag_me, cateout$Sigma) + geom_abline(slope = 1, intercept = 0)
    ## (cateout$Sigma / sig_diag_me)

    ## alpha_sigma_cate <- c(c(t(cateout$Gamma)), cateout$Sigma)
    ## em_miss_obj(alpha_sigma = alpha_sigma_cate, S = S, S2 = NULL, ncovs = 0,
    ##             ncontrols = ncontrols, n = n, p = p, k = k)
    ## em_miss_obj(alpha_sigma = alpha_sigma, S = S, S2 = NULL, ncovs = 0,
    ##             ncontrols = ncontrols, n = n, p = p, k = k)

    ## see if optim gives about the same results.
    lower_limits <- c(rep(-Inf, (k * p)), rep(10^-10, p))
    oout <- optim(par = alpha_sigma_new, fn = em_miss_obj, S = S, S2 = NULL, ncovs = 0,
                  ncontrols = ncontrols, n = n, p = p, k = k,
                  control = list(fnscale = -1),
                  method = "L-BFGS", lower = lower_limits)
    sig_diag_optim <- oout$par[((k * p) + 1):length(alpha_sigma)]
    alpha_optim <- oout$par[1:(k * p)]
    expect_equal(oout$par, alpha_sigma_new, tol = 10^-6)
}
)


test_that("em_miss_fix increase likelihood", {
    set.seed(92)
    n <- 11
    p <- 37
    ncontrols <- 13
    k <- 5
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    alpha_init <- t(pcout$alpha)
    sig_diag_init <- pcout$sig_diag
    alpha_sigma_init <- c(c(alpha_init), sig_diag_init)
    alpha_sigma <- alpha_sigma_init
    S2 <- crossprod(Y21)
    S  <- crossprod(cbind(Y31, Y32))

    ## lower_limits <- c(rep(-Inf, (k * p)), rep(10^-10, p))
    ## oout <- optim(par = alpha_sigma, fn = em_miss_obj, S = S, S2 = S2, ncovs = ncovs,
    ##               ncontrols = ncontrols, n = n, p = p, k = k, control = list(fnscale = -1),
    ##               method = "BFGS", lower = lower_limits)


    itermax <- 100
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- em_miss_obj(alpha_sigma = alpha_sigma, S = S, S2 = S2, ncovs = ncovs,
                                ncontrols = ncontrols, n = n, p = p, k = k)
    for(index in 2:itermax) {
        alpha_sigma_new <- em_miss_fix(alpha_sigma = alpha_sigma, S = S, S2 = S2, ncovs = ncovs,
                                       ncontrols = ncontrols, n = n, p = p, k = k)
        alpha_sigma <- alpha_sigma_new

        llike_vec[index] <- em_miss_obj(alpha_sigma = alpha_sigma, S = S, S2 = S2, ncovs = ncovs,
                                        ncontrols = ncontrols, n = n, p = p, k = k)
    }
    plot(llike_vec, type = "l")

    expect_true(all(llike_vec[1:(length(llike_vec) - 1)] <= llike_vec[2:length(llike_vec)]))

}
)


test_that("em_miss_fix_fast has same update as em_miss_fix", {
    set.seed(912)
    n <- 11
    p <- 79
    ncontrols <- 13
    k <- 5
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    alpha_init <- t(pcout$alpha)
    sig_diag_init <- pcout$sig_diag
    alpha_sigma_init <- c(c(alpha_init), sig_diag_init)
    alpha_sigma <- alpha_sigma_init
    S2 <- crossprod(Y21)
    S  <- crossprod(cbind(Y31, Y32))

    ## Sys.sleep(5)
    ## time.start <- proc.time()
    alpha_sigma_new <- em_miss_fix(alpha_sigma = alpha_sigma, S = S, S2 = S2, ncovs = ncovs,
                                   ncontrols = ncontrols, n = n, p = p, k = k)
    ## time.int <- proc.time() - time.start
    ## time.start <- proc.time()
    alpha_sigma_new_fast <- em_miss_fix_fast(alpha_sigma = alpha_sigma, Y21 = Y21,
                                             Y31 = Y31, Y32 = Y32, k = k)
    ## time.fast <- proc.time() - time.start
    ## plot(alpha_sigma_new, alpha_sigma_new_fast)
    ## abline(0, 1)

    expect_equal(alpha_sigma_new, alpha_sigma_new_fast)

    llike_old <- em_miss_obj(alpha_sigma = alpha_sigma, S = S, S2 = S2, ncovs = ncovs,
                             ncontrols = ncontrols, n = n, p = p, k = k)
    llike_fast <- em_miss_obj_fast(alpha_sigma = alpha_sigma, Y21 = Y21,
                                   Y31 = Y31, Y32 = Y32, k = k)

    expect_equal(llike_old, llike_fast)


    itermax <- 100
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- em_miss_obj_fast(alpha_sigma = alpha_sigma, Y21 = Y21,
                                     Y31 = Y31, Y32 = Y32, k = k)
    for (index in 1:itermax) {
        alpha_sigma <- em_miss_fix_fast(alpha_sigma = alpha_sigma, Y21 = Y21,
                                        Y31 = Y31, Y32 = Y32, k = k)
        llike_vec[index] <- em_miss_obj_fast(alpha_sigma = alpha_sigma, Y21 = Y21,
                                     Y31 = Y31, Y32 = Y32, k = k)
    }

    plot(llike_vec, type = "l")
    expect_true(all(llike_vec[1:(itermax - 1)] <= llike_vec[2:itermax]))

    emout <- em_miss(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k)

    llike_final <- em_miss_obj_fast(alpha_sigma = c(c(emout$alpha), emout$sig_diag), Y21 = Y21,
                                    Y31 = Y31, Y32 = Y32, k = k)

    llike_final
    max(llike_vec)

    expect_true(max(llike_vec) <= llike_final)
}
)


