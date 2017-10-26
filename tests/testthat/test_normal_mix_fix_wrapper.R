context("normal_mix_fix_wrapper works with error prone data")

test_that("normal_mix_fix_wrapper is ok", {
  err_dat <- readRDS("normal_mix_fix_wrapper_err_dat.RDS")
  out <- do.call(what = normal_mix_fix_wrapper, args = err_dat)
}
)
