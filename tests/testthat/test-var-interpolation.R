# NNS.VAR's interpolated_and_extrapolated values for mixed-frequency panels.

test_that("interior NAs are interpolated linearly when trailing NAs follow", {
  # A quarterly series in a monthly panel whose last two months are not yet
  # observed. Before the fix every NA, interior ones included, came from an
  # NNS.stack fit on the index, which lagged a trending series by several
  # quarters.
  set.seed(11)
  n <- 120
  t <- seq_len(n)
  q <- 100 + 0.5 * t + 2 * sin(2 * pi * t / 12) + rnorm(n, sd = 0.3)
  m <- 50 + cumsum(rnorm(n, 0.1, 0.5))
  keep <- seq(1, n - 2, by = 3)
  x <- rep(NA_real_, n)
  x[keep] <- q[keep]
  out <- NNS.VAR(data.frame(q = x, m = m), h = 0, tau = 3, status = FALSE, ncores = 1)

  last <- max(keep)
  inner <- setdiff(seq_len(last), keep)
  expect_identical(out$q[keep], x[keep])
  expect_equal(out$q[inner], approx(keep, x[keep], inner)$y, tolerance = 1e-12)
  expect_identical(out$m, m)

  # Trailing NAs are NNS.ARMA forecasts from the completed series.
  prefix <- approx(keep, x[keep], seq_len(last))$y
  expect_equal(out$q[(last + 1):n],
               NNS.ARMA(prefix, h = n - last, plot = FALSE, seasonal.plot = FALSE),
               tolerance = 1e-12)
  expect_true(all(abs(out$q[(last + 1):n] - x[last]) < 5))
})

test_that("a series too short to forecast repeats its last observation", {
  set.seed(12)
  m <- 50 + cumsum(rnorm(40))
  z <- c(3, 4, rep(NA, 38))
  out <- NNS.VAR(data.frame(z = z, m = m), h = 0, status = FALSE, ncores = 1)
  expect_identical(out$z, c(3, 4, rep(4, 38)))
})
