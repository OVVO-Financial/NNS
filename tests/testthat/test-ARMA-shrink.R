test_that("NNS.ARMA shrink ensembles with the means forecast for several seasonal factors", {
  set.seed(11)
  t <- 1:60
  z <- sin(2 * pi * t / 12) + rnorm(60, sd = 0.2)

  # Single seasonal factor (vectorized path) and several factors (recursive
  # path) must both honor shrink = TRUE; the recursive path used to ignore it.
  for (periods in list(12, c(12, 11))) {
    plain <- NNS.ARMA(z, h = 6, seasonal.factor = periods, method = "lin",
                      plot = FALSE, negative.values = TRUE)
    shrunk <- NNS.ARMA(z, h = 6, seasonal.factor = periods, method = "lin",
                       plot = FALSE, negative.values = TRUE, shrink = TRUE)
    expect_false(isTRUE(all.equal(plain, shrunk)))
  }
})
