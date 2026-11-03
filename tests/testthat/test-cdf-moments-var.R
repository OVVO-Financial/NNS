strict_joint_survival <- function(A, t) mean(rowSums(sweep(A, 2, t, ">")) == ncol(A))

test_that("NNS.moments sample skewness is the adjusted Fisher-Pearson G1", {
  expect_equal(NNS.moments(c(1, 2, 4, 8, 16), population = FALSE)$skewness, 1.325315, tolerance = 1e-6)
})

test_that("degree-0 survival for d = 3 matches strict joint counts with ties", {
  set.seed(7)
  A <- matrix(round(rnorm(90), 1), ncol = 3)
  S <- NNS.CDF(A, 0, type = "survival", plot = FALSE)$Function[[4]]
  expect_equal(S, apply(A, 1, strict_joint_survival, A = A), tolerance = 1e-12)
  t <- colMeans(A)
  expect_equal(NNS.CDF(A, 0, target = t, type = "survival", plot = FALSE)$target.value,
               strict_joint_survival(A, t), tolerance = 1e-12)
})

test_that("unsupported multivariate hazard never returns a CDF divided by survival", {
  set.seed(7)
  A <- matrix(rnorm(80), ncol = 2)
  expect_error(NNS.CDF(A, 0, target = A[5, ], type = "hazard", plot = FALSE),
               "joint density or coordinate-hazard definition")
  expect_error(NNS.CDF(A, 1, type = "hazard", plot = FALSE), "not implemented")
})

test_that("bivariate cumulative hazard target.value equals the row value at that target", {
  set.seed(7)
  A <- matrix(rnorm(80), ncol = 2)
  out <- NNS.CDF(A, 0, target = A[5, ], type = "cumulative hazard", plot = FALSE)
  expect_equal(out$target.value, out$Function[[3]][5], tolerance = 1e-12)
})

test_that("multivariate target is bounds-checked on every coordinate", {
  A <- cbind(c(0, 1, 2), c(0, 1, 2), c(0, 1, 2))
  expect_error(NNS.CDF(A, 0, target = c(0, 0, 1e6), type = "survival", plot = FALSE),
               "target out of bounds")
})

test_that("multivariate target must have one value per column", {
  A <- cbind(c(0, 1, 2), c(0, 1, 2), c(0, 1, 2))
  expect_error(NNS.CDF(A, 0, target = c(1, 1), plot = FALSE), "one value per column")
})

test_that("degree > 0 multivariate survival is the normalised joint UPM", {
  set.seed(7)
  A <- matrix(rnorm(90), ncol = 3)
  S <- NNS.CDF(A, 1, type = "survival", plot = FALSE)$Function[[4]]
  expect_equal(S, apply(A, 1, function(t) Co.UPM_nD(A, t, degree = 1)), tolerance = 1e-12)
})

test_that("degree > 0 survival is 0 with a coordinate at its maximum", {
  B <- cbind(c(0, 1, 2, 4), c(1, 0, 3, 2))
  expect_equal(NNS.CDF(B, 1, target = c(4, 1), type = "survival", plot = FALSE)$target.value, 0)
})

test_that("hazard is finite with tie runs at the minimum, interior and maximum", {
  x <- c(rep(0, 12), 1, 2, rep(3, 12), 4, 5, rep(6, 12))
  expect_true(all(is.finite(NNS.CDF(x, type = "hazard", plot = FALSE)$Function[[2]])))
})

test_that("hazard of a constant variable stops with a clear message", {
  expect_error(NNS.CDF(rep(2, 5), type = "hazard", plot = FALSE), "at least two distinct values")
})

test_that("LPM.ratio and UPM.ratio use integer targets", {
  x <- c(0, 1, 2, 3, 10)
  expect_equal(LPM.ratio(0, 1L, x), 0.4)
  expect_equal(UPM.ratio(1, c(1L, 2L), x), UPM.ratio(1, c(1, 2), x))
})

test_that("NNS.CDF of integer data matches the same data as double", {
  x <- rep(0:2, each = 20)
  expect_equal(NNS.CDF(x, plot = FALSE)$Function[[2]], rep(c(1, 2, 3) / 3, each = 20))
  expect_equal(NNS.CDF(x, type = "hazard", plot = FALSE)$Function[[2]],
               NNS.CDF(as.double(x), type = "hazard", plot = FALSE)$Function[[2]])
})

test_that("fractional-degree LPM.VaR inverts LPM.ratio", {
  set.seed(1)
  x <- rnorm(500) / 50
  p <- c(0.01, 0.05, 0.5, 0.95)
  for (degree in c(0.5, 1.5)) {
    expect_lt(max(abs(LPM.ratio(degree, LPM.VaR(p, degree, x), x) - p)), 1e-7)
  }
})

test_that("fractional-degree UPM.VaR inverts UPM.ratio", {
  set.seed(1)
  x <- rnorm(500) / 50
  p <- c(0.01, 0.05, 0.5, 0.95)
  for (degree in c(0.5, 1.5)) {
    expect_lt(max(abs(UPM.ratio(degree, UPM.VaR(p, degree, x), x) - p)), 1e-7)
  }
})

test_that("fractional-degree VaR returns the sample range at p = 0 and 1", {
  set.seed(1)
  x <- rnorm(500) / 50
  for (degree in c(0.5, 1.5)) {
    expect_identical(LPM.VaR(c(0, 1), degree, x), range(x))
    expect_identical(UPM.VaR(c(1, 0), degree, x), range(x))
  }
})

test_that("degree-0 joint survival with many columns is the joint upper share", {
  set.seed(9)
  x <- matrix(sample(0:4, 400, TRUE), 20, 20)
  s <- sapply(1:20, function(i) mean(colSums(t(x) > x[i, ]) == 20))
  h <- NNS.CDF(x, degree = 0, type = "cumulative hazard", plot = FALSE)$Function$CDF
  expect_equal(h, pmax(-log(pmax(s, 1e-10)), 0))
})

test_that("multivariate cumulative hazard needs degree 0 (#100)", {
  X <- cbind(c(0, 1, 2, 3), c(0, 2, 1, 3))
  expect_error(NNS.CDF(X, degree = 1, type = "hazard", plot = FALSE), "not implemented")
  expect_error(NNS.CDF(X, degree = 1, type = "cumulative hazard", plot = FALSE), "needs degree = 0")
  expect_true(is.list(NNS.CDF(X, degree = 0, type = "cumulative hazard", plot = FALSE)))
  # One dimension at degree 1 is a proper distribution function and keeps its hazards.
  set.seed(1)
  expect_true(all(is.finite(NNS.CDF(rnorm(50), degree = 1, type = "hazard", plot = FALSE)$Function[[2]])))
  # The degree-0 joint CDF gives every rectangle nonnegative mass.
  A <- matrix(rnorm(200), 100, 2)
  f <- function(t) NNS.CDF(A, degree = 0, target = t, plot = FALSE)$target.value
  expect_gte(f(c(.5, .5)) - f(c(-.5, .5)) - f(c(.5, -.5)) + f(c(-.5, -.5)), 0)
})
