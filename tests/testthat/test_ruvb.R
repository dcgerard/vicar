library(vicar)
context("ruvb")

test_that("See if RUV3 returns correct OLS estimates and matches vruv2", {
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

    return_list <- ruvb(Y = Y, X = X, ctl = ctl, k = q, cov_of_interest = cov_of_interest,
                        include_intercept = FALSE,
                        fa_args = list(nsamp = 100, keep = 1, print_update = FALSE))

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
    p <- 11
    ncontrols <- 7
    k <- 1
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

    ## bsvd_out <- bsvd(Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k,
    ##                  print_update = TRUE, plot_update = TRUE, nsamp = 1000, keep = 1)

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
