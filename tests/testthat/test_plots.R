context("Test Plots")

test_that("plot.mouthwash works", {
  mout <- readRDS(file = "norm_mout.RDS")
  plot(mout)
  mout <- readRDS(file = "unif_mout.RDS")
  plot(mout)
  bout <- readRDS(file = "bout.RDS")
  plot(bout)
  ruvbout <- readRDS(file = "ruvbout.RDS")
  ##plot(ruvbout)
}
)
