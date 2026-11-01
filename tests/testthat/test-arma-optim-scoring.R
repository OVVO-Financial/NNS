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
