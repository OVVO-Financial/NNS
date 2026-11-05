# NNS.VAR's interpolated_and_extrapolated values for mixed-frequency panels.

test_that("the index fit is NNS 13.0's one-neighbour fit on the duplicated index", {
  # Regression points 1, 7 and 13 give three cells. The newest observation
  # (13) is a cell of its own; the others are summarised by the gravity (here
  # the median) of their index and values: (2.5, 11) and (8.5, 15.5).
  x <- c(1, 4, 7, 10, 13)
  y <- c(10, 12, 15, 16, 20)
  fit <- .nns_var_index_fit(x, y, c(0, 2, 5, 8, 11, 14, 16), points = c(1, 7, 13))
  # Inside the data: the value of the nearest regression point.
  expect_equal(fit[2:5], c(11, 11, 15.5, 20))
  # After the data: from 20, the 3:2:1 average of the slopes to the centre of
  # the regression points (8.5) and to 10.75 and 11.875, whose nearest points
  # are 8.5 (the first of two tied) and 13: (3 * 1 + 2 * 2 + 0) / 6 per period.
  expect_equal(fit[6:7], 20 + c(1, 3) * 7 / 6, tolerance = 1e-12)
  # Before the data: from the first cell's 11, the same rule gives -0.3 per
  # period.
  expect_equal(fit[1], 10.7, tolerance = 1e-12)
})

test_that("with trailing NAs, every NA comes from the index fit", {
  # A quarterly series in a monthly panel whose last two months are not yet
  # observed.
  set.seed(11)
  n <- 120
  t <- seq_len(n)
  q <- 100 + 0.5 * t + 2 * sin(2 * pi * t / 12) + rnorm(n, sd = 0.3)
  m <- 50 + cumsum(rnorm(n, 0.1, 0.5))
  keep <- seq(1, n - 2, by = 3)
  x <- rep(NA_real_, n)
  x[keep] <- q[keep]
  out <- NNS.VAR(data.frame(q = x, m = m), h = 0, tau = 3, status = FALSE, ncores = 1)

  miss <- setdiff(t, keep)
  expect_identical(out$q[keep], x[keep])
  expect_identical(out$m, m)
  expect_equal(out$q[miss], .nns_var_index_fit(keep, x[keep], miss), tolerance = 1e-12)
  # The newest observation keeps its own regression point, so the trailing
  # values lie on a line through it.
  last <- max(keep)
  expect_equal(2 * out$q[last + 1] - out$q[last + 2], x[last], tolerance = 1e-10)
})

test_that("without trailing NAs, interior NAs are interpolated linearly", {
  set.seed(13)
  m <- 50 + cumsum(rnorm(30))
  z <- 10 + 0.3 * seq_len(30) + rnorm(30)
  z[c(5, 6, 17)] <- NA
  out <- NNS.VAR(data.frame(z = z, m = m), h = 0, status = FALSE, ncores = 1)
  obs <- which(!is.na(z))
  expect_equal(out$z[c(5, 6, 17)], approx(obs, z[obs], c(5, 6, 17))$y, tolerance = 1e-12)
  expect_identical(out$z[obs], z[obs])
})

test_that("a series with one or two observations is extended from them", {
  set.seed(12)
  m <- 50 + cumsum(rnorm(40))
  out <- NNS.VAR(data.frame(z = c(3, 4, rep(NA, 38)), m = m), h = 0, status = FALSE, ncores = 1)
  expect_equal(out$z, c(3, 4, 5:42))
  out <- NNS.VAR(data.frame(z = c(3, rep(NA, 39)), m = m), h = 0, status = FALSE, ncores = 1)
  expect_identical(out$z, rep(3, 40))
})
