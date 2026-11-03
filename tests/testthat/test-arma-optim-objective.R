mse <- expression(mean((predicted - actual)^2))

run_optim <- function(...) {
  NNS.ARMA.optim(AirPassengers[1:132], training.set = 105, seasonal.factor = c(6, 12, 24),
                 obj.fn = mse, ncores = 1, print.trace = FALSE, ...)
}

test_that("obj.fn scores the returned validation forecast", {
  o <- run_optim()
  expect_equal(o$validation.forecast, o$results)
  expect_equal(o$obj.fn, mean((o$results - AirPassengers[106:132])^2))
})

test_that("obj.fn scores the bias-shifted validation forecast when h is set", {
  o <- run_optim(h = 12)
  actual <- AirPassengers[106:132]
  validation <- pmax(0, o$errors + actual + o$bias.shift)
  expect_equal(o$obj.fn, mean((validation - actual)^2))
})

# Force a persistent bias, rather than relying on a random fixture selecting a
# particular gravity shift after other forecasting changes.
test_that("an accepted split-window correction is included in the returned score", {
  testthat::local_mocked_bindings(
    NNS.ARMA = function(variable, h, ...) rep(200, h),
    .package = "NNS"
  )
  o <- run_optim()
  actual <- AirPassengers[106:132]
  expect_true(o$bias.validation$accepted)
  expect_true(o$bias.shift != 0)
  expect_equal(o$raw.results, rep(200, length(actual)))
  expect_equal(o$results, pmax(0, o$raw.results + o$bias.shift))
  expect_equal(o$obj.fn, mean((o$results - actual)^2))
})
