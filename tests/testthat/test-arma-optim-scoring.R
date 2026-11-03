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

validation_forecast <- function(o, actual, negative.values = FALSE) {
  p <- o$errors + actual + o$bias.shift
  if (negative.values) p else pmax(0, p)
}

test_that("the lin model forecasts with negative.values", {
  t <- 1:160
  x <- pmax(0.1, 20 - 0.15 * t + sin(2 * pi * t / 12))
  o <- run_optim(x, training.set = 128, seasonal.factor = 12, lin.only = TRUE,
             negative.values = TRUE, obj.fn = mse)
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = 12, method = "lin",
                negative.values = TRUE, plot = FALSE)
  expect_equal(o$errors, p - tail(x, 32))
  expect_equal(o$obj.fn, mean((validation_forecast(o, tail(x, 32), TRUE) - tail(x, 32))^2))
})

test_that("linear.approximation = FALSE evaluates the nonlin method", {
  x <- AirPassengers[1:96]
  o <- run_optim(x, training.set = 84, seasonal.factor = c(6, 12), linear.approximation = FALSE)
  p <- NNS.ARMA(x, training.set = 84, h = 12, seasonal.factor = 12, method = "nonlin", plot = FALSE)
  expect_equal(o$method, "nonlin")
  expect_equal(as.numeric(o$periods), 12)
  expect_equal(o$errors, p - tail(x, 12))
  expect_equal(o$obj.fn, default_score(validation_forecast(o, tail(x, 12)), tail(x, 12)))
})

test_that("the both stage is scored with NNS.ARMA(method = \"both\")", {
  x <- random_walk(5)
  o <- run_optim(x, training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  p <- NNS.ARMA(x, training.set = 80, h = 20, seasonal.factor = 2, method = "both",
                shrink = o$shrink, plot = FALSE)
  expect_equal(o$method, "both")
  expect_equal(o$errors, p - tail(x, 20))
  expect_equal(o$obj.fn, mean((validation_forecast(o, tail(x, 20)) - tail(x, 20))^2))
})

test_that("equal weights are selected when they score better", {
  o <- run_optim(two_cycles(), h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
  expect_equal(as.numeric(o$periods), c(7, 12))
  expect_equal(o$weights, c(0.5, 0.5))
})


# A deterministic branch fixture tests smoothing data flow independently of
# changes to gravity or the empirically selected model on one random seed.
smoothing_fixture <- function(x, h = NULL) {
  calls <- integer()
  testthat::local_mocked_bindings(
    NNS.reg = function(x, y, ...) {
      calls <<- c(calls, length(y))
      list(Fitted.xy = data.frame(y.hat = rep(mean(y), length(y))))
    },
    NNS.ARMA = function(variable, h, training.set = NULL, ...) {
      if (!is.null(training.set)) variable <- head(variable, training.set)
      rep(tail(variable, 1), h)
    },
    .package = "NNS"
  )
  result <- run_optim(x, h = h, training.set = 80, seasonal.factor = 1,
                      obj.fn = mse, bias.correction = "none")
  list(result = result, calls = calls)
}

test_that("the out-of-sample forecast applies nns.regress to the full series", {
  x <- c(rep(10, 79), 100, rep(12, 20))
  fit <- smoothing_fixture(x, h = 20)
  expect_true(fit$result$nns.regress)
  expect_true(fit$result$shrink)
  expect_equal(fit$calls, c(80L, 100L))
  expect_equal(fit$result$results, rep(mean(x), 20))
})

test_that("obj.fn is the equal-weights score when equal weights win (min)", {
  x <- two_cycles()
  o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = o$method,
                weights = c(0.5, 0.5), shrink = o$shrink, plot = FALSE)
  expect_equal(o$errors, p - tail(x, 32))
  expect_equal(o$obj.fn, mean((validation_forecast(o, tail(x, 32)) - tail(x, 32))^2))
})

test_that("obj.fn is the equal-weights score when equal weights win (max)", {
  x <- two_cycles()
  o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), objective = "max",
             obj.fn = expression(-mean((predicted - actual)^2)))
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = o$method,
                weights = c(0.5, 0.5), shrink = o$shrink, plot = FALSE)
  expect_equal(o$errors, p - tail(x, 32))
  expect_equal(o$obj.fn, -mean((validation_forecast(o, tail(x, 32)) - tail(x, 32))^2))
})

test_that("nns.regress is not selected by a smooth that saw the validation window", {
  x <- c(rep(0, 79), 14.5, rep(8, 20))
  # A full-series mean would beat the last-value forecast; the prefix mean loses.
  expect_lt((mean(x) - 8)^2, (14.5 - 8)^2)
  expect_gt((mean(x[1:80]) - 8)^2, (14.5 - 8)^2)
  fit <- smoothing_fixture(x)
  expect_false(fit$result$nns.regress)
  expect_equal(fit$calls, 80L)
  expect_equal(fit$result$results, rep(14.5, 20))
})

test_that("validation results with nns.regress use a smooth of the training window", {
  x <- c(rep(10, 79), 100, rep(12, 20))
  fit <- smoothing_fixture(x)
  expect_true(fit$result$nns.regress)
  expect_equal(fit$calls, 80L)
  expect_equal(fit$result$results, rep(mean(x[1:80]), 20))
})

test_that("nns.regress must beat the accepted shrink score", {
  # Smoothed 1.670 beats raw 3.809 but loses to raw with shrink, 1.447.
  o <- run_optim(random_walk(5), training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
  expect_true(o$shrink)
  expect_false(o$nns.regress)
})

trend_cycles <- function() {
  set.seed(1)
  t <- 1:160
  50 + 0.05 * t + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
}

shrink_vs_equal_weights <- function(objective, obj.fn) {
  x <- trend_cycles()
  o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), lin.only = TRUE, objective = objective, obj.fn = obj.fn)
  p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = "lin",
                weights = c(0.5, 0.5), plot = FALSE)
  expect_equal(o$weights, c(0.5, 0.5))
  expect_false(o$shrink)
  expect_equal(o$errors, p - tail(x, 32))
}

# Equal weights score 12.04 and no weights 12.67, so shrink at 12.40 must lose to the equal-weights score.
test_that("shrink is tested against the equal-weights score (min)", {
  shrink_vs_equal_weights("min", mse)
})

test_that("shrink is tested against the equal-weights score (max)", {
  shrink_vs_equal_weights("max", expression(-mean((predicted - actual)^2)))
})
