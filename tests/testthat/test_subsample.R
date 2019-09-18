context("subsample MOUTHWASH")


test_that("subsample gives same results when TRUE and all genes are subsampled", {
  skip("travis is a little finicky with this one")
  set.seed(1)
  set.seed(69)
  n <- 20
  p <- 102
  k <- 3
  cov_of_interest <- 2
  X <- matrix(stats::rnorm(n * k), nrow = n)
  beta <- matrix(stats::rnorm(k * p), nrow = k)
  beta[, 1:round(p/2)] <- 0
  E <- matrix(stats::rnorm(n * p), nrow = n)
  Y <- X %*% beta + E

  include_intercept <- FALSE
  limmashrink       <- TRUE
  fa_func           <- pca_naive
  fa_args           <- list()
  likelihood        <- "t"
  mixing_dist       <- "uniform"
  degrees_freedom   <- NULL
  pi_init           <- NULL
  tau_seq           <- NULL
  lambda_seq        <- NULL
  lambda0           <- 10

  mout1 <- mouthwash(Y = Y, X = X, k = k, cov_of_interest = 2, include_intercept = FALSE,
                     subsample = TRUE, num_sub = 102)
  mout2 <- mouthwash(Y = Y, X = X, k = k, cov_of_interest = 2, include_intercept = FALSE,
                     subsample = FALSE)
  plot(mout1$result$PosteriorMean, mout2$result$PosteriorMean)
  abline(0, 1)
  expect_true(all(abs(mout1$result$PosteriorMean - mout2$result$PosteriorMean) < 10 ^ -2))
}
)
