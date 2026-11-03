mse <- expression(mean((predicted - actual)^2))

run_optim <- function(...) {
  NNS.ARMA.optim(AirPassengers[1:132], training.set = 105, seasonal.factor = c(6, 12, 24),
                 obj.fn = mse, ncores = 1, print.trace = FALSE, ...)
}

test_that("obj.fn scores the returned validation forecast", {
  o <- run_optim()
  expect_true(o$bias.shift != 0)
  expect_equal(o$obj.fn, mean((o$results - AirPassengers[106:132])^2))
})

test_that("obj.fn scores the bias-shifted validation forecast when h is set", {
  o <- run_optim(h = 12)
  actual <- AirPassengers[106:132]
  validation <- pmax(0, o$errors + actual + o$bias.shift)
  expect_equal(o$obj.fn, mean((validation - actual)^2))
})
