context("Hierarchical Prior")

test_that("hier_fun works ok", {
    ## Generate test data
    pstar <- 73
    n <- 3
    beta_mat <- matrix(rnorm(n * pstar), nrow = n)
    shape_param = 1
    rate_param = 1

    hier_fun(beta_mat)

}
)
