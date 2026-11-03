lin_arma <- function(x, ...) NNS.ARMA(x, h = 3, seasonal.factor = 1, method = "lin", plot = FALSE, ...)

test_that("held-out negative values do not change the clamp", {
  expect_equal(lin_arma(c(20:1, -1), training.set = 20), lin_arma(20:1))
  expect_equal(lin_arma(20:1), c(0, 0, 0))
  expect_true(all(lin_arma(c(20:1, -1)) < 0))
})

test_that("NNS.ARMA.optim clamps validation forecasts from a nonnegative training window", {
  x <- c(seq(58, 2, -2), 1, 0.5, 0.2, 0.1, -0.1)
  run <- function(...) NNS.ARMA.optim(x, training.set = 29, seasonal.factor = c(1, 2), ncores = 1,
                                      print.trace = FALSE, obj.fn = expression(mean((predicted - actual)^2)), ...)
  o <- run()
  expect_true(all(o$errors + tail(x, 5) >= 0))
  expect_true(all(o$results >= 0))
  # The out-of-sample forecast is fitted on the whole series, which is negative.
  expect_true(any(run(h = 3)$results < 0))
})
