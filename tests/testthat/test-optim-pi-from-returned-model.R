mse <- expression(mean((predicted - actual)^2))

two_cycles <- function() {
  set.seed(9)
  t <- 1:160
  50 + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
}

run_optim <- function(...) NNS.ARMA.optim(..., ncores = 1, print.trace = FALSE)

test_that("errors and bias.shift come from the returned shrink model", {
  x <- two_cycles()
  o <- run_optim(x, training.set = 128, seasonal.factor = c(7, 12), obj.fn = mse,
                 bias.correction = "split")
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = "lin",
                weights = c(0.5, 0.5), shrink = TRUE, plot = FALSE)
  expect_true(o$shrink)
  expect_equal(o$errors, p - tail(x, 32))
  expect_equal(o$bias.validation$estimated.shift,
               -NNS.gravity(o$errors[o$bias.validation$estimate.rows]))
  expect_equal(o$bias.shift, o$bias.validation$applied.shift)
  expect_equal(o$results, p + o$bias.shift)
})

test_that("errors come from the returned nns.regress model", {
  x <- as.numeric(AirPassengers[1:132])
  o <- run_optim(x, training.set = 105, seasonal.factor = c(6, 12, 24), objective = "max")
  smoothed <- NNS.reg(1:105, x[1:105], plot = FALSE, smooth = TRUE)$Fitted.xy$y.hat
  p <- NNS.ARMA(smoothed, training.set = 105, h = 27, seasonal.factor = o$periods,
                method = o$method, weights = o$weights, shrink = TRUE, plot = FALSE)
  expect_true(o$nns.regress)
  expect_equal(o$errors, p - tail(x, 27))
  expect_equal(o$results, p + o$bias.shift)
})

test_that("the interval adds the LPM and UPM VaR of actual - results", {
  x <- two_cycles()
  o <- run_optim(x, training.set = 128, seasonal.factor = c(7, 12), obj.fn = mse)
  r <- tail(x, 32) - o$results
  expect_equal(o$lower.pred.int, o$results + LPM.VaR(0.025, 0, r))
  expect_equal(o$upper.pred.int, o$results + UPM.VaR(0.025, 0, r))
})

test_that("the out-of-sample interval uses the validation residuals", {
  x <- two_cycles()
  o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
  r <- -(o$errors + o$bias.shift)
  expect_equal(o$lower.pred.int, o$results + LPM.VaR(0.025, 0, r))
  expect_equal(o$upper.pred.int, o$results + UPM.VaR(0.025, 0, r))
})

test_that("a nonnegative series takes residuals from the clamped forecasts it returns", {
  # A persistent positive residual centre shifts alternating raw forecasts
  # downward; the lower member must be clamped before the bands are computed.
  testthat::local_mocked_bindings(
    NNS.ARMA = function(variable, h, ...) rep(c(1, 9), length.out = h),
    .package = "NNS"
  )
  v <- rep(1, 80)
  o <- NNS.ARMA.optim(v, seasonal.factor = c(2, 4), training.set = 64, negative.values = FALSE,
                      obj.fn = expression(mean(predicted)), objective = "min",
                      print.trace = FALSE, ncores = 1, plot = FALSE,
                      bias.correction = "split")
  expect_true(any(o$results == 0))
  res <- tail(v, 16) - o$results
  expect_equal(o$upper.pred.int, pmax(0, o$results + UPM.VaR(.025, 0, res)))
  expect_equal(o$lower.pred.int, pmax(0, o$results + LPM.VaR(.025, 0, res)))
  # Residuals of the unclamped forecasts move the bands on this series.
  unclamped <- -(o$errors + o$bias.shift)
  expect_gt(max(abs(c(pmax(0, o$results + UPM.VaR(.025, 0, unclamped)) - o$upper.pred.int,
                      pmax(0, o$results + LPM.VaR(.025, 0, unclamped)) - o$lower.pred.int))), 0.01)
})

