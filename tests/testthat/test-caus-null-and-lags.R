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
  expect_identical(NNS:::.nns_caus_shift_offsets(100, 2L), 10:90)
  expect_identical(NNS:::.nns_caus_shift(1:100, 30L), c(31:100, 1:30))
})

test_that("the shift null uses each admissible offset at most once", {
  # n = 100, tau = 49 admits only offset 50. Drawn 19 times with replacement it
  # gave p = 1 / 20 whenever the observed statistic beat it.
  p <- vapply(1:10, function(s) {
    set.seed(s); x <- rnorm(100); y <- rnorm(100)
    NNS.caus(x, y, tau = 49, p.value = TRUE, nperm = 19, seed = s)$p.value$two.sided
  }, numeric(1))
  expect_true(all(p >= 0.5))
  set.seed(1); M <- matrix(rnorm(120), 60, 2)
  r <- NNS.caus(M, tau = 29, p.value = TRUE, nperm = 19, seed = 1)
  expect_gte(r$p.value_matrix[1, 2], 0.5)
  set.seed(2); d <- NNS:::.nns_caus_shift_draw(100, 2L, 19)
  expect_length(d, 19)
  expect_false(anyDuplicated(d) > 0)
  expect_true(all(d %in% 10:90))
})

test_that("the matrix p-value path keeps factor columns numeric", {
  r <- NNS.caus(iris, p.value = TRUE, nperm = 2, seed = 1)
  expect_equal(dim(r$causality), c(5, 5))
  expect_true(all(is.finite(r$p.value_matrix)))
})

test_that("the matrix time-series null shifts one series per pair", {
  # Shifting every column independently can leave a pair aligned when the two
  # shifts nearly coincide; each pair now keeps its first column and shifts
  # the second, so the shift helper runs once per pair and permutation.
  calls <- 0L
  testthat::local_mocked_bindings(
    .nns_caus_shift = function(v, tau = 0L) { calls <<- calls + 1L; n <- length(v); c(v[31:n], v[1:30]) },
    .package = "NNS")
  set.seed(5); M <- matrix(rnorm(240), 60, 4)
  r <- NNS.caus(M, tau = 1, p.value = TRUE, nperm = 2, seed = 1)
  expect_identical(calls, 2L * ncol(utils::combn(4, 2)))   # 12; 8 with column shifts
  expect_true(all(is.finite(r$p.value_matrix)))
})

test_that("the matrix p-value honours permute like the vector path", {
  set.seed(107); M <- matrix(rnorm(120), 60, 2)
  for (tau in c(0, 1)) for (m in c("x", "y", "both", "shift")) {
    mat <- NNS.caus(M, tau = tau, p.value = TRUE, nperm = 19, permute = m, seed = 7)$p.value_matrix[1, 2]
    vec <- NNS.caus(M[, 1], M[, 2], tau = tau, p.value = TRUE, nperm = 19, permute = m, seed = 7)$p.value$two.sided
    expect_equal(mat, vec, info = paste(tau, m))
  }
})

test_that("the matrix p-value uses the same signed statistic as the vector path on tied data", {
  set.seed(17); invisible(rnorm(120))
  M <- cbind(sample(0:1, 60, TRUE), sample(0:1, 60, TRUE))
  mat <- NNS.caus(M, p.value = TRUE, nperm = 19, permute = "both", seed = 7)$p.value_matrix[1, 2]
  vec <- NNS.caus(M[, 1], M[, 2], p.value = TRUE, nperm = 19, permute = "both", seed = 7)$p.value$two.sided
  expect_equal(mat, vec)
})
