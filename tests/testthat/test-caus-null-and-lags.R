# NNS.caus after #96: aligned dim.red weights, a null centred for any margins,
# one lag per permutation test, and circular shifts for time series.

test_that("dim.red with ts.test weights the aligned pairs the model fits", {
  set.seed(9)
  X <- cbind(signal = rnorm(200), noise = rnorm(200))
  y <- X[, 1]^2 + rnorm(200, sd = .1)
  fit <- NNS.stack(X[1:160, ], y[1:160], IVs.test = X[161:200, ], ts.test = 32,
                   folds = 1, method = 2, stack = FALSE, dim.red.method = "NNS.caus",
                   order = 5, optimize.threshold = FALSE, ncores = 1, status = FALSE)
  w <- fit$NNS.dim.red.coefficients
  expect_gt(w[1], 4 * w[2])                               # 0 / 0.0186 before
  expect_lt(sqrt(mean((fit$dim.red - y[161:200])^2)), 0.5)  # 1.618 before
  # cs and ts give the same weight: both score aligned pairs.
  expect_identical(NNS:::Uni.caus(y, X[, 1], "ts"), NNS:::Uni.caus(y, X[, 1], "cs"))
})

test_that("independent series with different tie structures centre at 0.5", {
  set.seed(71)
  v <- replicate(20, NNS.caus(rbinom(1000, 1, .05), rnorm(1000))[[1]])
  expect_lt(abs(mean(v) - 0.5), 0.05)                     # 0.408 before
  set.seed(72)
  w <- replicate(20, NNS.caus(rexp(1000), sample(1:5, 1000, TRUE))[[1]])
  expect_lt(abs(mean(w) - 0.5), 0.05)
})

test_that("P is unbiased under independence for a binary effect", {
  # The degree-2 deviation with the m / (m - 1) factor is a U-statistic, so
  # the within-cell average is unbiased whatever the effect's distribution.
  set.seed(73)
  p <- replicate(200, {
    cause <- NNS:::.nns_caus_u(rnorm(500)); effect <- NNS:::.nns_caus_u(rbinom(500, 1, .1))
    n <- length(cause); k <- floor(sqrt(n)); cell <- pmax(1L, ceiling(cause * k))
    d <- vapply(split(effect, cell), NNS:::.nns_caus_dev, numeric(1)); w <- lengths(split(effect, cell))
    1 - sum(w * d) / sum(w) / NNS:::.nns_caus_dev(effect)
  })
  expect_lt(abs(mean(p)), 0.01)
  v <- c(3, 1, 4, 1, 5)
  expect_equal(NNS:::.nns_caus_dev(v), var(v))
  # About its centre, an unresolved sign is not resolved: knowing y = x^2
  # leaves x's deviation in place.
  set.seed(74); x <- rnorm(1000)
  r <- NNS:::.nns_caus_raw(x, x^2, 0L)
  expect_gt(r$P[1], 0.8)
  expect_lt(r$P[2], 0.1)
})

test_that("tau = \"ts\" with p.value resolves the lag once", {
  set.seed(1); x <- as.numeric(arima.sim(list(ar = .9), 50)); y <- rnorm(50)
  r <- NNS.caus(x, y, tau = "ts", p.value = TRUE, nperm = 20, seed = 1)
  expect_equal(r$causation, NNS.caus(x, y, tau = "ts"))
  expect_identical(r$p.value$permute, "shift")
  expect_length(r$p.value$null_distribution, 20)
})

test_that("time series use a circular-shift null, cross-sections a shuffle", {
  set.seed(2); x <- rnorm(100); y <- c(0, x[-100]) + rnorm(100)
  expect_identical(NNS.caus(x, y, tau = 1, p.value = TRUE, nperm = 5, seed = 1)$p.value$permute, "shift")
  expect_identical(NNS.caus(x, y, p.value = TRUE, nperm = 5, seed = 1)$p.value$permute, "y")
  # A shift is a rotation at least max(tau + 1, n / 10) steps from alignment.
  set.seed(3); v <- 1:100
  for (i in 1:50) {
    s <- NNS:::.nns_caus_shift(v, 2L)
    k <- which(s == 1) - 1L                              # s starts at v[k + 1]
    expect_identical(s, c(v[(101 - k):100], v[1:(100 - k)]))
    expect_true(min(k, 100 - k) >= 10)
  }
})

test_that("the matrix p-value path keeps factor columns numeric", {
  r <- NNS.caus(iris, p.value = TRUE, nperm = 2, seed = 1)
  expect_equal(dim(r$causality), c(5, 5))
  expect_true(all(is.finite(r$p.value_matrix)))
})
