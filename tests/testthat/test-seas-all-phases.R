cv <- function(v) sd(v) / mean(v)

test_that("NNS.seas ranks period 12 first on AirPassengers", {
  r <- NNS.seas(as.numeric(AirPassengers), plot = FALSE)
  expect_equal(r$best.period, 12L)
  expect_equal(r$periods, c(12L, 13L, 6L))
  expect_equal(r$all.periods$Coefficient.of.Variation,
               c(0.412018727459241, 0.426207016328105, 0.427938478152998), tolerance = 1e-12)
})

test_that("NNS.seas reproduces the ch25 lag-2 toy example", {
  r <- NNS.seas(c(10, 18, 11, 21, 12, 20, 13, 23), plot = FALSE)
  expect_equal(r$periods, 2L)
  expect_equal(r$all.periods$Coefficient.of.Variation,
               (4 * cv(c(10, 11, 12, 13)) + 4 * cv(c(18, 21, 20, 23))) / 8)
})

test_that("each lag m scores all m component series, weighted by length", {
  x <- as.numeric(lynx)
  n <- length(x)
  ref <- sapply(1:(n %/% 2 - 1), function(m) {
    sum(sapply(1:m, function(j) { v <- x[seq(j, n, m)]; length(v) * cv(v) })) / n
  })
  pass <- which(ref < cv(x))
  r <- NNS.seas(x, plot = FALSE)$all.periods
  expect_equal(sort(r$Period), pass)
  expect_equal(r$Coefficient.of.Variation[order(r$Period)], ref[pass])
})

test_that("a series with mean zero up to rounding compares standard deviations", {
  x <- sin(2 * pi * (1:96) / 12)
  x <- x - mean(x)
  r <- NNS.seas(x, plot = FALSE)$all.periods
  expect_equal(r$Variable.Coefficient.of.Variation[1], sd(x))
  expect_equal(sort(r$Period[1:3]), c(12L, 24L, 36L))
  expect_lt(r$Coefficient.of.Variation[3], 1e-12)
})

test_that("no lag is reported unless its CV is strictly below CV(X)", {
  expect_equal(NNS.seas(rep(5, 20), plot = FALSE)$periods, 1L)
})

test_that("the no-seasonality row gets weight 1, not 0", {
  M <- data.frame(Period = 1L, Coefficient.of.Variation = 0.5, Variable.Coefficient.of.Variation = 0.5)
  expect_equal(NNS:::ARMA.seas.weighting(FALSE, M)$Weights, 1)
})

test_that("NNS.ARMA forecasts a nonseasonal trend from lag 1", {
  f <- NNS.ARMA(1:100, h = 2, plot = FALSE, seasonal.plot = FALSE)
  expect_equal(f, NNS.ARMA(1:100, h = 2, seasonal.factor = 1, plot = FALSE))
  expect_equal(f, c(101, 102), tolerance = 1e-3)
})

test_that("a constant zero phase scores 0 and keeps its lag", {
  x <- rep(c(0, 10), 10)
  r <- NNS.seas(x, plot = FALSE)$all.periods
  expect_equal(sort(r$Period), c(2L, 4L, 6L, 8L))
  expect_equal(r$Coefficient.of.Variation, rep(0, 4))
  expect_equal(NNS.ARMA(x, h = 2, plot = FALSE, seasonal.plot = FALSE), c(0, 10))
})
