library(vicar)
context("leapp")

test_that("vicarius_leapp works", {
    skip("no leapp yet")
    set.seed(98)
    n <- 13
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

    leapp_out <- vicarius_leapp(Y = Y, X = X, k = num_sv)

}
)

test_that("leapp_cross works", {
    skip("no leapp yet")
    set.seed(50)
    p <- 1000
    k <- 3
    X <-matrix(rnorm(k * p), nrow = p)
    alpha <- matrix(rnorm(k), ncol = 1)
    E <- matrix(rnorm(p), ncol = 1)

    gamma <- succotashr::draw_beta(pi_vals = c(0.8, 0.1, 0.1),
                                   tau_seq = c(0, 1, 2),
                                   p = p)

    Y <- X %*% alpha + E + gamma


    H <- X %*% solve(t(X) %*% X) %*% t(X)

    ## vout <- vicarius_principis_IPOD(X = X, Y = Y, H = H, sigma = 1)

    lout <- leapp::IPOD(X = X, Y = Y, H = H, method = "soft")
    vout <- cleapp(X = X, Y = Y, H = H, method = "soft")

    sum(vout$gamma == 0)

    sum((lout$gamma - gamma) ^ 2)
    sum((vout$gamma - gamma) ^ 2)
}
)
