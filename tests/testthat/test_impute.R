context("factanal and other imputation functions")


test_that("factanal_wrapper works", {
    set.seed(76)
    n <- 11
    p <- 79
    ncontrols <- 13
    k <- 5
    ncovs <- 3

    Y21 <- matrix(rnorm(ncovs * ncontrols), nrow = ncovs)
    Y31 <- matrix(rnorm((n - ncovs) * ncontrols), nrow = n - ncovs)
    Y32 <- matrix(rnorm((n - ncovs) * (p - ncontrols)), nrow = n - ncovs)

}
)
