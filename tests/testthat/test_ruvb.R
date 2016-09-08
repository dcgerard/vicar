library(vicar)
context("ruvb")

test_that("see if ruvb works ok. In particular, if uniform prior returns same results the two ways I calculate it", {
    set.seed(71)
    n <- 11
    p <- 71
    k <- 3
    q <- 2

    X <- matrix(rnorm(n * q), nrow = n)
    beta <- matrix(rnorm(q * p), nrow = q)
    beta[, 1:29] <- 0
    Z <- matrix(rnorm(n * k), nrow = n)
    alpha <- matrix(rnorm(k * p), nrow = k)
    E <- matrix(rnorm(n * p), nrow = n)
    Y <- X %*% beta + Z %*% alpha + E
    ctl <- rep(FALSE, length = p)
    ctl[1:13] <- TRUE
    include_intercept <- FALSE
    cov_of_interest <- 2
    fa_args <- list()
    fa_func <- bfl

    set.seed(73)
    return_list <- ruvb(Y = Y, X = X, ctl = ctl, k = q,
                        cov_of_interest = cov_of_interest,
                        include_intercept = FALSE, prior_fun = NULL,
                        fa_args = list(nsamp = 1000, thin = 1, display_progress = FALSE))

    set.seed(73) ## returns identity
    return_list2 <- ruvb(Y = Y, X = X, ctl = ctl, k = q,
                         cov_of_interest = cov_of_interest,
                         include_intercept = FALSE,
                         prior_fun = function(beta_mat){ 1 },
                         return_log = FALSE,
                         fa_args = list(nsamp = 1000, thin = 1, display_progress = FALSE))



    expect_equal(return_list$lfsr1, return_list2$lfsr1)
    expect_equal(return_list$means, return_list2$means)
    expect_equal(return_list$sd, return_list2$sd, tol = 10 ^ -3)
    expect_equal(return_list$lfsr2, return_list2$lfsr2, tol = 10 ^ -3)

    expect_true(sum(abs(return_list$lower - return_list2$lower)) /
                sum(abs(return_list$lower)) < 0.02)
    expect_true(sum(abs(return_list$upper - return_list2$upper)) /
                sum(abs(return_list$upper)) < 0.02)



    cout <- cate::cate.fit(X.primary = X[, cov_of_interest, drop = FALSE],
                           X.nuis = X[, -cov_of_interest, drop = FALSE],
                           Y = Y, r = q, adj.method = "nc", nc = ctl)
    cate_sebetahat <- c(sqrt(cout$beta.cov.row * cout$beta.cov.col) /
                        sqrt(nrow(X)))

    ## library(ggplot2)
    ## dat <- data.frame(beta = c(beta[cov_of_interest, !ctl]),
    ##                   lower = c(return_list$posterior_lower),
    ##                   upper = c(return_list$posterior_upper),
    ##                   median = c(return_list$posterior_medians),
    ##                   mean = c(return_list$posterior_means),
    ##                   lfsr = c(return_list$lfsr),
    ##                   cate_betahat = c(t(cout$beta[!ctl, ])),
    ##                   cate_lower = c(t(cout$beta[!ctl, ])) - 2 * cate_sebetahat[!ctl],
    ##                   cate_upper = c(t(cout$beta[!ctl, ])) + 2 * cate_sebetahat[!ctl])

    ## qplot(x = 1:nrow(dat), y = median, data = dat) +
    ##     geom_errorbar(mapping = aes(ymin = lower, ymax = upper)) +
    ##     geom_vline(xintercept = sum(c(beta[cov_of_interest, !ctl]) == 0) + 0.5,
    ##                col = 2, lty = 2) +
    ##     geom_point(mapping = aes(y = cate_betahat, x = 1:nrow(dat)), color = I("red")) +
    ##     geom_errorbar(mapping = aes(ymin = cate_lower, ymax = cate_upper), color = I("red"),
    ##                   alpha = I(0.5)) +
    ##     geom_point(mapping = aes(y = beta, x = 1:nrow(dat)), color = I("blue"), alpha = I(0.5))


    ## sum((dat$beta - dat$mean) ^ 2)
    ## sum((dat$beta - dat$cate_betahat) ^ 2)

    ## lfsr_order <- order(dat$lfsr)
    ## qplot(x = 1:nrow(dat), y = dat$lfsr[lfsr_order], color = (dat$beta[lfsr_order] != 0))

}
)


test_that("bfl works OK", {
    set.seed(81)
    n <- 9
    p <- 23
    ncontrols <- 7
    k <- 3
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

    ## bsvd_out <- bsvd(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k,
    ##                  print_update = TRUE, plot_update = TRUE, nsamp = 1000, keep = 1)

    ## nsamp = 10000
    ## burnin = round(nsamp / 4); keep = 20
    ## print_update = TRUE; plot_update = FALSE
    ## hetero_factors = FALSE; rho_0 = 0.1; alpha_0 = 0.1
    ## delta_0 = 0.1; lambda_0 = 0.1; nu_0 = 1; beta_0 = 1
    ## eta_0 = 1; tau_0 = 1

    bfout   <- bfa_wrapper(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k,
                           print_status = 100000)
    bfl_out <- bfl(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k,
                   print_update = FALSE, plot_update = FALSE,
                   nsamp = 100, keep = 1)
    gdout <- gdfa(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k, nsamp = 100,
                  keep = 1, print_update = FALSE)
    expect_equal(dim(bfl_out$Y22_array)[1:2], c(ncovs, p - ncontrols))

    ## hist(bfl_out$xi_mat[, 1])
    ## hist(colMeans(1 / bfl_out$xi_mat))
    ## mean(1 / bfl_out$psi_phi[, 2])
    ## bfl_out$Y22_array


}
)

test_that("bfa_gd_gibbs works ok", {
dat <- readRDS("bfa_gd_examp.Rds")
Y22out <- bfa_gd_gibbs(Linit = dat$Linit, Finit = dat$Finit,
                       xi_init = dat$xi_init, phi_init = dat$phi_init,
                       zeta_init = dat$zeta_init, theta_init = dat$theta_init,
                       kappa_init = dat$kappa_init, Y22init = dat$Y22init,
                       Y21 = dat$Y21, Y31 = dat$Y31, Y32 = dat$Y32,
                       nsamp = dat$nsamp, burnin = dat$burnin, thin = dat$thin,
                       rho_0 = 1, alpha_0 = 1, delta_0 = 1, lambda_0 = 1,
                       nu_0 = 1, beta_0 = 1, eta_0 = 1, tau_0 = 1,
                       hetero_factors = TRUE, display_progress = FALSE)

dat2 <- readRDS("bfa_gd_examp2.Rds")
Y22out <- bfa_gd_gibbs(Linit = dat2$Linit, Finit = dat2$Finit,
                       xi_init = dat2$xi_init, phi_init = dat2$phi_init,
                       zeta_init = dat2$zeta_init, theta_init = dat2$theta_init,
                       kappa_init = dat2$kappa_init, Y22init = dat2$Y22init,
                       Y21 = dat2$Y21, Y31 = dat2$Y31, Y32 = dat2$Y32,
                       nsamp = dat2$nsamp, burnin = dat2$burnin,
                       thin = dat2$thin, rho_0 = dat2$rho_0,
                       alpha_0 = dat2$alpha_0, delta_0 = dat2$delta_0,
                       lambda_0 = dat2$lambda_0, nu_0 = dat2$nu_0,
                       beta_0 = dat2$beta_0, eta_0 = dat2$eta_0,
                       tau_0 = dat2$tau_0, hetero_factors = TRUE,
                       display_progress = FALSE)

Y22out <- bfa_gs_linked_gibbs(Linit = dat2$Linit, Finit = dat2$Finit,
                              xi_init = dat2$xi_init, phi_init = dat2$phi_init,
                              zeta_init = dat2$zeta_init, Y22init = dat2$Y22init,
                              Y21 = dat2$Y21, Y31 = dat2$Y31, Y32 = dat2$Y32,
                              nsamp = dat2$nsamp, burnin = dat2$burnin,
                              thin = dat2$thin, rho_0 = dat2$rho_0,
                              alpha_0 = dat2$alpha_0,
                              beta_0 = dat2$beta_0, eta_0 = dat2$eta_0,
                              tau_0 = dat2$tau_0,
                              display_progress = FALSE)

Y22outr <- bfa_gs_linked_gibbs_r(Linit = dat2$Linit, Finit = dat2$Finit,
                                 xi_init = dat2$xi_init, phi_init = dat2$phi_init,
                                 zeta_init = dat2$zeta_init, Y22init = dat2$Y22init,
                                 Y21 = dat2$Y21, Y31 = dat2$Y31, Y32 = dat2$Y32,
                                 nsamp = dat2$nsamp, burnin = dat2$burnin,
                                 thin = dat2$thin, rho_0 = dat2$rho_0,
                                 alpha_0 = dat2$alpha_0,
                                 beta_0 = dat2$beta_0, eta_0 = dat2$eta_0,
                                 tau_0 = dat2$tau_0,
                                 display_progress = FALSE)
})



test_that("calc_lfsr_g works", {
    y <- -10:10
    g <- rep(1, 20)
    cout <- calc_lfsr_g(y, g)
    expect_equal(cout, 0.5)

    g <- c(rep(0, 5), rep(1, 15))
    cout <- calc_lfsr_g(y, g)
    expect_equal(cout, 1/3)
}
)

test_that("inputting own prior results in same posterior summaries for log and no log", {
    set.seed(71)
    n <- 11
    p <- 37
    k <- 3
    q <- 2

    X <- matrix(rnorm(n * q), nrow = n)
    beta <- matrix(rnorm(q * p), nrow = q)
    beta[, 1:29] <- 0
    Z <- matrix(rnorm(n * k), nrow = n)
    alpha <- matrix(rnorm(k * p), nrow = k)
    E <- matrix(rnorm(n * p), nrow = n)
    Y <- X %*% beta + Z %*% alpha + E
    ctl <- rep(FALSE, length = p)
    ctl[1:13] <- TRUE
    include_intercept <- FALSE
    cov_of_interest <- 2
    fa_args <- list()
    fa_func <- bfl

    set.seed(73)
    return_list1 <- ruvb(Y = Y, X = X, ctl = ctl, k = q,
                         cov_of_interest = cov_of_interest,
                         include_intercept = FALSE,
                         prior_fun = hier_fun, return_log = TRUE,
                         fa_args = list(nsamp = 1000, thin = 1, display_progress = FALSE))

    set.seed(73)
    return_list2 <- ruvb(Y = Y, X = X, ctl = ctl, k = q,
                         cov_of_interest = cov_of_interest,
                         include_intercept = FALSE,
                         prior_fun = hier_fun, return_log = FALSE,
                         fa_args = list(nsamp = 1000, thin = 1, display_progress = FALSE),
                         prior_args = list(return_log = FALSE))

    expect_equal(return_list1, return_list2)
}
)
