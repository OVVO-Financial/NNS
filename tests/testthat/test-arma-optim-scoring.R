mse <- expression(mean((predicted - actual)^2))

default_score <- function(predicted, actual) eval(eval(formals(NNS.ARMA.optim)$obj.fn))

random_walk <- function(seed) {
  set.seed(seed)
  30 + cumsum(rnorm(100))
}

two_cycles <- function() {
  set.seed(9)
  t <- 1:160
  50 + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
}

run_optim <- function(...) NNS.ARMA.optim(..., ncores = 1, print.trace = FALSE)

test_that("the lin search scores forecasts with negative.values", {
  t <- 1:160
  x <- pmax(0.1, 20 - 0.15 * t + sin(2 * pi * t / 12))
  o <- run_optim(x, training.set = 128, seasonal.factor = 12, lin.only = TRUE,
             negative.values = TRUE, obj.fn = mse)
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = 12, method = "lin",
                negative.values = TRUE, plot = FALSE)
  expect_equal(o$obj.fn, mean((p - tail(x, 32))^2))
})

test_that("linear.approximation = FALSE evaluates the nonlin method", {
  x <- AirPassengers[1:96]
  o <- run_optim(x, training.set = 84, seasonal.factor = c(6, 12), linear.approximation = FALSE)
  p <- NNS.ARMA(x, training.set = 84, h = 12, seasonal.factor = 12, method = "nonlin", plot = FALSE)
  expect_equal(o$method, "nonlin")
  expect_equal(as.numeric(o$periods), 12)
  expect_equal(o$obj.fn, default_score(p, tail(x, 12)))
})

test_that("the both stage is scored with NNS.ARMA(method = \"both\")", {
  x <- random_walk(5)
  o <- run_optim(x, training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  p <- NNS.ARMA(x, training.set = 80, h = 20, seasonal.factor = 2, method = "both", plot = FALSE)
  expect_equal(o$method, "both")
  expect_equal(o$obj.fn, mean((p - tail(x, 20))^2))
})

test_that("equal weights are selected when they score better", {
  o <- run_optim(two_cycles(), h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
  expect_equal(as.numeric(o$periods), c(7, 12))
  expect_equal(o$weights, c(0.5, 0.5))
})

test_that("the out-of-sample forecast applies nns.regress to the full series", {
  x <- random_walk(8)
  o <- run_optim(x, h = 20, training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  smoothed <- NNS.reg(1:100, x, plot = FALSE, smooth = TRUE)$Fitted.xy$y.hat
  p <- NNS.ARMA(smoothed, h = 20, seasonal.factor = o$periods, method = o$method,
                weights = o$weights, shrink = TRUE, plot = FALSE)
  expect_true(o$nns.regress)
  expect_true(o$shrink)
  expect_equal(o$results, p + o$bias.shift)
})

test_that("obj.fn is the equal-weights score when equal weights win (min)", {
  x <- two_cycles()
  o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = o$method,
                weights = c(0.5, 0.5), plot = FALSE)
  expect_equal(o$obj.fn, mean((p - tail(x, 32))^2))
})

test_that("obj.fn is the equal-weights score when equal weights win (max)", {
  x <- two_cycles()
  o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), objective = "max",
             obj.fn = expression(-mean((predicted - actual)^2)))
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = o$method,
                weights = c(0.5, 0.5), plot = FALSE)
  expect_equal(o$obj.fn, -mean((p - tail(x, 32))^2))
})

test_that("nns.regress is not selected by a smooth that saw the validation window", {
  # Smoothing all 100 points scores 19.4 against 29.8 raw; smoothing the 80 training points scores 56.8.
  o <- run_optim(random_walk(4), training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  expect_false(o$nns.regress)
})

test_that("validation results with nns.regress use a smooth of the training window", {
  x <- random_walk(8)
  o <- run_optim(x, training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  smoothed <- NNS.reg(1:80, x[1:80], plot = FALSE, smooth = TRUE)$Fitted.xy$y.hat
  p <- NNS.ARMA(smoothed, training.set = 80, h = 20, seasonal.factor = o$periods,
                method = o$method, weights = o$weights, shrink = TRUE, plot = FALSE)
  expect_true(o$nns.regress)
  expect_equal(o$results, p + o$bias.shift)
})

test_that("nns.regress must beat the accepted shrink score", {
  # Smoothed 1.670 beats raw 3.809 but loses to raw with shrink, 1.447.
  o <- run_optim(random_walk(5), training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  expect_true(o$shrink)
  expect_false(o$nns.regress)
})
