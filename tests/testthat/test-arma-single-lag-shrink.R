shrunk_lin <- function(h, ...) {
  NNS.ARMA(AirPassengers, h = h, method = "lin", shrink = TRUE, plot = FALSE, ...)
}

shrunk_recursion <- function(x, lag, h) {
  for (j in seq_len(h)) {
    y <- rev(x[seq(length(x) + 1 - lag, 1, -lag)])
    i <- seq_along(y)
    fit <- coef(lm(y ~ i))
    x <- c(x, max(0, (fit[[1]] + fit[[2]] * (length(y) + 1) + mean(y)) / 2))
  }
  tail(x, h)
}

test_that("a single shrunk lag refits on its own forecasts past the lag", {
  single <- shrunk_lin(30, seasonal.factor = 12)
  expect_equal(single, shrunk_lin(30, seasonal.factor = c(12, 24), weights = c(1, 0)))
  expect_equal(single, shrunk_recursion(as.numeric(AirPassengers), 12, 30))
})

test_that("a single shrunk lag is unchanged up to the lag", {
  expect_equal(shrunk_lin(1, seasonal.factor = 12), 332.0568182, tolerance = 1e-6)
  expect_equal(shrunk_lin(12, seasonal.factor = 12), shrunk_recursion(as.numeric(AirPassengers), 12, 12))
})

test_that("a single shrunk lag longer than the history forecasts NaN, not 0", {
  expect_equal(NNS.ARMA(rep(100, 11), seasonal.factor = 12, h = 3, method = "lin", shrink = TRUE, plot = FALSE),
               c(NaN, 100, 100))
})
