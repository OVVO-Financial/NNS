# Reference for the scan (#97): the linear trend is removed (keeping the mean),
# each component's dispersion is the Gini mean difference, the coefficient of
# variation is used when the series stays on one side of zero, and lags are
# scanned only up to n / 5 so every component has at least 5 points (#69).
gmd <- function(v) {
  k <- length(v)
  if (k < 2 || min(v) == max(v)) return(0)
  s <- sort(v)
  2 * sum((2 * seq_len(k) - k - 1) * s) / (k * (k - 1))
}
detrend <- function(x) { t <- seq_along(x); x - unname(coef(lm(x ~ t))[2]) * (t - mean(t)) }
score <- function(x, m) {
  n <- length(x); use_cv <- all(x > 0) || all(x < 0)
  sum(sapply(1:m, function(j) {
    v <- x[seq(j, n, m)]
    length(v) * (if (use_cv) gmd(v) / abs(mean(v)) else gmd(v))
  })) / n
}

test_that("NNS.seas ranks period 12 first on AirPassengers", {
  r <- NNS.seas(as.numeric(AirPassengers), plot = FALSE)
  expect_equal(r$best.period, 12L)
  x <- detrend(as.numeric(AirPassengers))
  expect_equal(r$all.periods$Coefficient.of.Variation[1], score(x, 12), tolerance = 1e-10)
  expect_equal(r$all.periods$Variable.Coefficient.of.Variation[1], score(x, 1), tolerance = 1e-10)
})

test_that("NNS.seas reproduces the ch25 lag-2 toy example once components have 5 points", {
  # The 8-point book series has 4 points per lag-2 component: below the minimum.
  expect_equal(NNS.seas(c(10, 18, 11, 21, 12, 20, 13, 23), plot = FALSE)$periods, 1L)
  x <- c(10, 18, 11, 21, 12, 20, 13, 23, 14, 22)
  r <- NNS.seas(x, plot = FALSE)
  expect_equal(r$periods, 2L)
  expect_equal(r$all.periods$Coefficient.of.Variation, score(detrend(x), 2), tolerance = 1e-10)
  # Without detrending the components are the book's (10, ..., 14) and (18, ..., 22).
  r0 <- NNS.seas(x, plot = FALSE, detrend = FALSE)
  expect_equal(r0$all.periods$Coefficient.of.Variation, score(x, 2), tolerance = 1e-12)
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
  x <- detrend(as.numeric(lynx))
  n <- length(x)
  ref <- sapply(1:(n %/% 5), function(m) score(x, m))
  pass <- which(ref < ref[1])
  r <- NNS.seas(as.numeric(lynx), plot = FALSE)$all.periods
  expect_equal(sort(r$Period), pass)
  expect_equal(r$Coefficient.of.Variation[order(r$Period)], ref[pass], tolerance = 1e-10)
})

test_that("a series that crosses zero compares dispersions, not CVs", {
  x <- sin(2 * pi * (1:96) / 12)
  x <- x - mean(x)
  r <- NNS.seas(x, plot = FALSE, detrend = FALSE)$all.periods
  expect_equal(r$Variable.Coefficient.of.Variation[1], gmd(x))
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
  r <- NNS.seas(x, plot = FALSE, detrend = FALSE)$all.periods
  expect_equal(sort(r$Period), c(2L, 4L))
  expect_equal(r$Coefficient.of.Variation, rep(0, 2))
  expect_equal(NNS.ARMA(x, h = 2, plot = FALSE, seasonal.plot = FALSE), c(0, 10))
})

test_that("a constant decimal series has no seasonal period", {
  expect_equal(NNS.seas(rep(.1, 20), plot = FALSE)$best.period, 1)
  expect_equal(NNS.seas(rep(c(.1, 1), 10), plot = FALSE)$best.period, 2)
})

test_that("a linear trend does not hide the seasonal period (#97)", {
  bp <- sapply(1201:1210, function(s) {
    set.seed(s); t <- 1:240
    NNS.seas(10 + .07 * t + 2 * sin(2 * pi * t / 12) + rnorm(240), plot = FALSE)$best.period
  })
  expect_true(all(bp %% 12 == 0))                 # 0 of 10 before
  expect_equal(NNS.seas(1:50 + 0.5, plot = FALSE)$best.period, 1)
})

test_that("an exactly linear series detrends to a constant, not a period", {
  for (a in c(0, 10, -50, 1e6, 1e9)) for (b in c(.3, -2, 1e-3)) for (n in c(12L, 50L, 144L)) {
    r <- NNS.seas(a + b * seq_len(n), plot = FALSE)
    expect_equal(r$best.period, 1, info = paste(a, b, n))
    expect_true(all(r$all.periods$Coefficient.of.Variation >= 0), info = paste(a, b, n))
  }
  r <- NNS.seas(10 + .3 * (1:50), plot = FALSE)$all.periods
  expect_equal(c(r$Coefficient.of.Variation, r$Variable.Coefficient.of.Variation), c(0, 0))
  t <- 1:120
  expect_equal(NNS.seas(1e6 + .3 * t + sin(2 * pi * t / 12), plot = FALSE)$best.period, 12)
  # Small or slow periods on a high, long series are not rounding noise.
  t <- seq_len(12000)
  expect_equal(NNS.seas(1e9 + .3 * t + .01 * sin(2 * pi * t / 12), plot = FALSE, max.period = 12)$best.period, 12)
  x <- 1e9 + .3 * t + .1 * cos(2 * pi * (t - .5) / 1200)
  expect_equal(NNS.seas(x, plot = FALSE, max.period = 1200)$best.period, 1200)
  expect_equal(NNS.seas(1e9 + cumsum(rep(.3, 1e5)), plot = FALSE, max.period = 24)$best.period, 1)
})

test_that("noise flags do not depend on the level or the tail shape (#97)", {
  flagged <- function(gen) mean(sapply(1:5, function(s) {
    set.seed(s); sum(NNS.seas(gen(), plot = FALSE)$periods > 1)
  }))
  # Before: 170 of 200 at mean 0, 0 at mean 1, 198 for t3 noise at mean 10.
  for (v in c(flagged(function() rnorm(1000)), flagged(function() 1 + rnorm(1000)),
              flagged(function() 10 + rt(1000, 3))))
    expect_true(v > 60 && v < 140)
})

test_that("p.value keeps periods that beat the best lag of shuffled copies (#97)", {
  fp <- sapply(1:10, function(s) {
    set.seed(s); NNS.seas(10 + rnorm(300), plot = FALSE, p.value = TRUE, nperm = 40)$best.period
  })
  expect_lte(sum(fp > 1), 2)
  set.seed(3); t <- 1:240
  r <- NNS.seas(10 + 2 * sin(2 * pi * t / 12) + rnorm(240), plot = FALSE, p.value = TRUE, nperm = 40)
  expect_equal(r$best.period %% 12, 0)
  expect_true(all(r$all.periods$p.value[match(r$periods, r$all.periods$Period)] <= 0.05))
})

test_that("max.period caps the scan", {
  set.seed(4); x <- 10 + rnorm(500)
  expect_true(all(NNS.seas(x, plot = FALSE, max.period = 7)$all.periods$Period <= 7))
  expect_error(NNS.seas(x, plot = FALSE, max.period = 0), "max.period")
})

test_that("p.value tests modulo periods beyond the scan against a null that scores them", {
  # Each shuffled copy scores the same modulo-added periods as the observed
  # series; with a null that only scanned up to n / 5, noise passed in 27 of
  # 40 series here.
  hits <- sapply(1:20, function(s) {
    set.seed(s)
    a <- NNS.seas(10 + rnorm(96), modulo = 20:47, p.value = TRUE, nperm = 40, plot = FALSE)
    any(a$periods > 1)
  })
  expect_lte(sum(hits), 3)
  set.seed(3)
  x <- 10 + sin(2 * pi * (1:96) / 24) + rnorm(96, sd = .3)
  a <- NNS.seas(x, modulo = 24, mod.only = FALSE, p.value = TRUE, nperm = 40, plot = FALSE)
  expect_equal(a$best.period, 24L)
})
