ap <- as.numeric(AirPassengers)

arma <- function(x, ...) NNS.ARMA(x, plot = FALSE, ...)

one_step_path <- function(x, h, ...) {
  out <- numeric(h)
  for (i in seq_len(h)) {
    out[i] <- arma(x, h = 1, ...)
    x <- c(x, out[i])
  }
  out
}

test_that("method = 'means' forecasts the seasonal subset means", {
  monthly <- as.numeric(tapply(ap, rep(1:12, 12), mean))
  expect_equal(arma(ap, h = 14, seasonal.factor = 12, method = "means"),
               rep(monthly, length.out = 14))
})

test_that("method = 'means' with several lags weights each lag's subset mean", {
  expected <- 0.25 * mean(ap[seq(1, 144, 12)]) + 0.75 * mean(ap[seq(1, 144, 24)])
  expect_equal(arma(ap, h = 1, seasonal.factor = c(12, 24), weights = c(1, 3),
                    method = "means"),
               expected)
})

test_that("nonlinear forecasts keep negative values when the series has them", {
  set.seed(1)
  s <- 3 * sin(2 * pi * (1:240) / 12) + rnorm(240, 0, 0.5)
  f <- arma(s, h = 12, seasonal.factor = 12)
  expect_lt(max(f[7:11]), -0.5)
})

test_that("forecasts do not depend on the scale of the weights", {
  for (m in c("lin", "nonlin", "both")) {
    expect_equal(arma(ap, h = 3, seasonal.factor = c(12, 24), weights = c(1, 1), method = m),
                 arma(ap, h = 3, seasonal.factor = c(12, 24), weights = c(0.5, 0.5), method = m))
  }
})

test_that("a single-lag weight is relative", {
  expect_equal(arma(ap, h = 3, seasonal.factor = 12, weights = 0.5),
               arma(ap, h = 3, seasonal.factor = 12))
})

test_that("weights need one value per seasonal lag", {
  expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = 1),
               "one value per seasonal lag")
})

test_that("weights must be finite and must not sum to 0", {
  expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = c(0, 0)), "finite")
  expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = c(NA, 1)), "finite")
})

test_that("weights = 'equal' works with a numeric seasonal.factor", {
  expect_equal(arma(ap, h = 3, seasonal.factor = c(12, 24), weights = "equal"),
               arma(ap, h = 3, seasonal.factor = c(12, 24), weights = c(1, 1)))
})

test_that("seasonal.factor = TRUE with best.periods = NULL uses the best period", {
  best <- NNS.seas(ap, plot = FALSE)$best.period
  expect_equal(arma(ap, h = 3, seasonal.factor = TRUE, best.periods = NULL),
               arma(ap, h = 3, seasonal.factor = best))
})

test_that("a constant series forecasts the constant with automatic seasonality", {
  expect_equal(arma(rep(5, 60), h = 3), c(5, 5, 5))
})

test_that("all-zero seasonal weights fall back to equal weights", {
  set.seed(1)
  x <- rnorm(100, 10)
  expect_equal(arma(x, h = 2, modulo = 7),
               arma(x, h = 2, modulo = 7, weights = "equal"))
})

test_that("a constant series forecasts the constant with a numeric seasonal.factor", {
  for (m in c("lin", "nonlin", "both", "means")) {
    expect_equal(arma(rep(5, 60), h = 3, seasonal.factor = 12, method = m), c(5, 5, 5))
    expect_equal(arma(rep(5, 60), h = 3, seasonal.factor = c(12, 24), method = m), c(5, 5, 5))
    expect_equal(arma(rep(-2, 30), h = 2, seasonal.factor = 7, method = m), c(-2, -2))
  }
})

test_that("dynamic = TRUE matches repeated one-step forecasts with best.periods", {
  expect_equal(arma(ap, h = 3, seasonal.factor = FALSE, best.periods = 2, dynamic = TRUE),
               one_step_path(ap, 3, seasonal.factor = FALSE, best.periods = 2))
})

test_that("dynamic = TRUE honours modulo", {
  expect_equal(arma(ap, h = 3, seasonal.factor = FALSE, best.periods = 2, modulo = 12,
                    dynamic = TRUE),
               one_step_path(ap, 3, seasonal.factor = FALSE, best.periods = 2, modulo = 12))
})

test_that("dynamic = TRUE uses periods that appear later in the forecast path", {
  set.seed(1)
  x <- rnorm(20, 10)
  expect_equal(arma(x, h = 6, seasonal.factor = FALSE, best.periods = 4, method = "lin",
                    dynamic = TRUE),
               one_step_path(x, 6, seasonal.factor = FALSE, best.periods = 4, method = "lin"))
})

test_that("dynamic = TRUE rejects numeric weights", {
  expect_error(arma(ap, h = 2, seasonal.factor = FALSE, best.periods = 2, dynamic = TRUE,
                    weights = c(1, 0)),
               "dynamic")
})

test_that("dynamic = TRUE applies weights = 'equal' at every step", {
  expect_equal(arma(ap, h = 3, seasonal.factor = FALSE, best.periods = 2, dynamic = TRUE,
                    weights = "equal"),
               one_step_path(ap, 3, seasonal.factor = FALSE, best.periods = 2, weights = "equal"))
})

test_that("h must be a positive integer", {
  for (h in list(0, 2.5, -1, NA, c(1, 2))) {
    expect_error(arma(ap, h = h, seasonal.factor = 12), "positive integer")
  }
  expect_length(arma(ap, h = 3L, seasonal.factor = 12), 3)
})

test_that("method = 'means' returns prediction intervals", {
  set.seed(123)
  r <- arma(ap, h = 4, seasonal.factor = 12, method = "means", pred.int = 0.95)
  expect_equal(dim(r), c(4L, 3L))
  expect_true(all(is.finite(as.matrix(r))))
  expect_true(all(r[[2]] <= r[[1]] & r[[1]] <= r[[3]]))
})

test_that("the warn option is restored when NNS.ARMA fails", {
  old <- getOption("warn")
  on.exit(options(warn = old), add = TRUE)
  options(warn = 1)
  expect_error(arma(ap[1:3], h = 2, seasonal.factor = 12))
  expect_equal(getOption("warn"), 1)
  expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = 1))
  expect_equal(getOption("warn"), 1)
})
