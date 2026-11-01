direct_lpm <- function(d, t, x) vapply(t, function(tt) LPM(d, tt, x), numeric(1))
direct_upm <- function(d, t, x) vapply(t, function(tt) UPM(d, tt, x), numeric(1))

test_that("vectorised LPM/UPM (> 32 targets) match the direct scan at every degree", {
  set.seed(7)
  price <- 400 * cumprod(1 + rnorm(1500, 4e-4, 0.012))
  ret <- rt(1000, df = 3) * 0.01
  for (x in list(price, ret)) {
    t <- c(seq(min(x) - 0.1 * diff(range(x)), max(x) + 0.1 * diff(range(x)), length.out = 80),
           sort(x)[seq(1, length(x), by = 25)])
    for (d in c(0, 1, 2, 3, 4, 8, 12)) {
      l <- LPM(d, t, x); u <- UPM(d, t, x)
      expect_true(all(l >= 0) && all(u >= 0))
      expect_equal(l, direct_lpm(d, t, x), tolerance = 1e-10)
      expect_equal(u, direct_upm(d, t, x), tolerance = 1e-10)
      r <- LPM.ratio(d, t, x)
      expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
    }
  }
  # Degree 2 just below the maximum used to come out negative.
  x <- rnorm(1e4, 100, 1)
  t <- sort(x)[(1e4 - 99):1e4] - 1e-9
  expect_true(all(UPM(2, t, x) >= 0))
})

test_that("negative and NaN degrees are rejected", {
  expect_error(LPM(-1, 0, c(-2, 1)), "non-negative")
  expect_error(UPM(NaN, 0, c(-2, 1)), "non-negative")
  expect_error(LPM.ratio(-0.5, 0, c(-2, 1)), "non-negative")
  expect_error(Co.LPM(-1, c(-2, 1), c(-2, 1), 0, 0), "non-negative")
  expect_error(D.UPM(1, -2, c(-2, 1), c(-2, 1), 0, 0), "non-negative")
  expect_error(PM.matrix(-1, 1, NULL, matrix(rnorm(20), 10), TRUE), "non-negative")
  expect_error(Co.LPM_nD(matrix(rnorm(20), 10), c(0, 0), -1), "non-negative")
  expect_equal(LPM(0.5, 0, c(-4, 1)), 1)
})

test_that("discrete NNS.mode handles values beyond the int range", {
  expect_equal(NNS.mode(c(3e9, 3e9, 3e9, 1, 2), discrete = TRUE), 3e9)
  expect_equal(NNS.mode(c(-5e9, -5e9, -5e9, 1, 2), discrete = TRUE), -5e9)
  expect_equal(NNS.mode(c(1, 2, 2, 3, 3, 4, 4, 5), discrete = TRUE, multi = TRUE), c(2, 3, 4))
})

test_that("NNS.gravity / NNS.mode stay fast with a tight IQR and an outlier", {
  set.seed(1)
  x <- c(5 + 1e-7 * rnorm(1000), 6)
  elapsed <- system.time({ g <- NNS.gravity(x); m <- NNS.mode(x, multi = FALSE) })[["elapsed"]]
  expect_lt(elapsed, 2)
  expect_equal(g, 5, tolerance = 1e-3)
  expect_equal(m, 5, tolerance = 1e-5)
  # Tiny histograms (<= 3 bins) are well defined.
  expect_equal(NNS.mode(c(0, 0, 1, 1), multi = FALSE), 0.5)
  expect_equal(NNS.gravity(c(0, 0, 1, 1)), 0.375)
})

test_that("vector generators do not read out of bounds or divide by zero", {
  g <- NNS:::generate.vectors(c(1, 2, 3), 7L)
  expect_length(g$Component.series[[1]], 0)
  expect_equal(NNS:::generate.vectors(as.numeric(1:10), 3L)$Component.series[[1]], c(2, 5, 8))
  expect_equal(NNS:::generate.vectors(as.numeric(1:10), 3L, 7L)$Component.series[[1]], c(2, 5))
  expect_error(NNS:::generate.vectors(c(1, 2, 3), 2L, 5L), "len")
  expect_error(NNS:::generate.lin.vectors(1:4, 0L), "positive")
  lv <- NNS:::generate.lin.vectors(c(1, 2, 3), 7L, 2L)
  expect_length(lv$Component.series[[1]], 0)
})
