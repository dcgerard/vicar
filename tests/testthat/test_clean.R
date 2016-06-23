library(vicar)
context("Cleaning Test Statistics")


test_that("clean_my_p works", {
    p <- 2000
    pvalues <- stats::rbeta(n = p, shape1 = 1/2, shape2 = 2)
    ctl <- rep(FALSE, length = p)
    ctl[1:100] <- TRUE

    ksout <- smell_my_p(pvalues, ctl)
    cleaned_p <- clean_my_p(pvalues, ctl)
    ksout_cleaned <- smell_my_p(cleaned_p, ctl)

    expect_true(ksout$p.value < ksout_cleaned$p.value)
}
)

test_that("clean_my_stats works", {
    set.seed(99)
    p <- 100
    stats <- stats::rbeta(n = p, shape1 = 1/2, shape2 = 2)
    ctl <- rep(FALSE, length = p)
    ctl[1:100] <- TRUE
    kout <- smell_my_stats(stats = stats, ctl = ctl, cdf_func = qt, df = 4)
    cleaned_stats <- clean_my_stats(stats = stats, ctl = ctl, qfunc = qt, df = 4)
}
)
