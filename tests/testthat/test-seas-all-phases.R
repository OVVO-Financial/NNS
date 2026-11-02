# Each component sd is divided by c4(k), and lags are scanned only up to n / 5
# so every component has at least 5 points (#69).
c4 <- function(k) sqrt(2 / (k - 1)) * exp(lgamma(k / 2) - lgamma((k - 1) / 2))
cv <- function(v) sd(v) / c4(length(v)) / mean(v)

test_that("NNS.seas ranks period 12 first on AirPassengers", {
  r <- NNS.seas(as.numeric(AirPassengers), plot = FALSE)
  expect_equal(r$best.period, 12L)
  expect_equal(r$periods, 12L)
  expect_equal(r$all.periods$Coefficient.of.Variation, 0.421476943243904, tolerance = 1e-12)
})

test_that("NNS.seas reproduces the ch25 lag-2 toy example once components have 5 points", {
  # The 8-point book series has 4 points per lag-2 component: below the minimum.
  expect_equal(NNS.seas(c(10, 18, 11, 21, 12, 20, 13, 23), plot = FALSE)$periods, 1L)
  r <- NNS.seas(c(10, 18, 11, 21, 12, 20, 13, 23, 14, 22), plot = FALSE)
  expect_equal(r$periods, 2L)
  expect_equal(r$all.periods$Coefficient.of.Variation,
               (5 * cv(c(10, 11, 12, 13, 14)) + 5 * cv(c(18, 21, 20, 23, 22))) / 10)
})

test_that("NNS.seas scans only lags with at least 5 points per component", {
  set.seed(69)
  for (n in c(60, 200, 1000)) {
    r <- NNS.seas(rnorm(n, 10), plot = FALSE)
    expect_true(all(r$all.periods$Period <= n %/% 5))
  }
  # Pure noise no longer flags nearly every lag (493 of 500 before).
  set.seed(1)
  flagged <- replicate(5, sum(NNS.seas(rnorm(1000, 10), plot = FALSE)$all.periods$Period > 1))
  expect_lt(mean(flagged), 120)
})

test_that("NNS.seas finds the true period of a noisy seasonal series", {
  # The best period is the true period or one of its multiples (which share
  # its seasonal structure); before #69 it was a long non-multiple lag.
  for (p in c(7, 12, 24)) for (s in 1:5) {
    set.seed(100 * p + s); x <- 10 + 2 * sin(2 * pi * (1:240) / p) + rnorm(240)
    r <- NNS.seas(x, plot = FALSE)
    expect_equal(r$best.period %% p, 0)
    expect_true(p %in% r$periods)
  }
})

test_that("each lag m scores all m component series, weighted by length", {
  x <- as.numeric(lynx)
  n <- length(x)
  ref <- sapply(1:(n %/% 5), function(m) {
    sum(sapply(1:m, function(j) { v <- x[seq(j, n, m)]; length(v) * cv(v) })) / n
  })
  pass <- which(ref < sd(x) / c4(n) / mean(x))
  r <- NNS.seas(x, plot = FALSE)$all.periods
  expect_equal(sort(r$Period), pass)
  expect_equal(r$Coefficient.of.Variation[order(r$Period)], ref[pass])
})

test_that("a series with mean zero up to rounding compares standard deviations", {
  x <- sin(2 * pi * (1:96) / 12)
  x <- x - mean(x)
  r <- NNS.seas(x, plot = FALSE)$all.periods
  expect_equal(r$Variable.Coefficient.of.Variation[1], sd(x) / c4(96))
  # 24 and 36 have fewer than 5 points per component at n = 96.
  expect_equal(r$Period[1], 12L)
  expect_lt(r$Coefficient.of.Variation[1], 1e-12)
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
  expect_equal(sort(r$Period), c(2L, 4L))
  expect_equal(r$Coefficient.of.Variation, rep(0, 2))
  expect_equal(NNS.ARMA(x, h = 2, plot = FALSE, seasonal.plot = FALSE), c(0, 10))
})

test_that("a constant decimal series has no seasonal period", {
  expect_equal(NNS.seas(rep(.1, 20), plot = FALSE)$best.period, 1)
  expect_equal(NNS.seas(rep(c(.1, 1), 10), plot = FALSE)$best.period, 2)
})
