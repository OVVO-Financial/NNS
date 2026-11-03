test_that("NNS.dep is 1 for identical binary data at any group size", {
  b <- rep(0:1, c(51, 49))
  expect_equal(unname(unlist(NNS.dep(b, b))), c(1, 1))
})

test_that("NNS.dep keeps the sign of a fully discordant cell", {
  x <- c(-2, -1, 1, 2)
  expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1))
  expect_equal(unname(unlist(NNS.dep(x, -x, asym = TRUE))), c(-1, 1))
  expect_equal(NNS.dep(cbind(x, -x))$Correlation[1, 2], -1)
})

test_that("NNS.copula does not return 1 for an interior pair sum of 1", {
  X <- cbind(c(5, 3, 4, 5, 3, 3, 2, 5, 2, 2),
             c(2, 5, 1, 1, 5, 5, 5, 2, 2, 1),
             c(1, 1, 4, 2, 4, 5, 2, 4, 4, 5))
  expect_equal(NNS.copula(X), 0.577972677175, tolerance = 1e-10)
})

test_that("NNS.copula is 1 for perfect dependence with mixed signs", {
  x <- 1:100
  expect_equal(NNS.copula(cbind(x, x, -x)), 1)
  expect_equal(NNS.copula(cbind(x, x, -x), continuous = FALSE), 1)
  expect_equal(NNS.copula(cbind(x, x, x, -x)), 1)
})

test_that("degree-0 nD co-moments partition the sample when points tie the target", {
  X <- cbind(c(1, 2, 3, 3, 4, 5), c(3, 3, 1, 5, 3, 2))
  t <- c(3, 3)
  expect_equal(c(Co.LPM_nD(X, t), Co.UPM_nD(X, t), DPM_nD(X, t)), c(0.5, 0, 0.5))
  G <- as.matrix(expand.grid(-1:1, -1:1, -1:1))
  t <- c(0, 0, 0)
  expect_equal(c(Co.LPM_nD(G, t), Co.UPM_nD(G, t), DPM_nD(G, t)), c(8, 1, 18) / 27)
})

test_that("NNS.dep degree-0 discordant count excludes ties with the cell mean", {
  # The middle point ties both cell means, so it counts as co-lower and the
  # cell is not all-discordant. Three points then carry no evidence beyond
  # the within-cell permutation null (#70), so the result is 0, not -1.
  x <- -1:1
  expect_equal(unname(unlist(NNS.dep(x, -x))), c(0, 0))
  # With five points the countermonotone signal clears the null: the debiased
  # degree-0 and degree-1 terms average to 2/3 (no square root since #95).
  x <- -2:2
  expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1) * 2 / 3)
})

test_that("NNS.dep reads an exactly independent tied grid as 0", {
  # The discrete blend used to add 1/6 here (#95).
  G <- as.matrix(expand.grid(-1:1, -1:1))[rep(1:9, each = 100), ]
  expect_equal(unname(unlist(NNS.dep(G[, 1], G[, 2]))), c(0, 0))
  expect_equal(unname(unlist(NNS.dep(G[, 1], G[, 2], asym = TRUE))), c(0, 0))
})

test_that("symmetric NNS.dep on tied data ignores argument and row order", {
  # Zero covariances and values equal to a cell mean rounded either way
  # depending on the order the cell was summed in, more so at a high level.
  set.seed(487); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
  expect_equal(unlist(NNS.dep(y, x)), unlist(NNS.dep(x, y)), tolerance = 1e-12)
  for (lev in c(0, 1e9)) for (div in c(1, 10)) {
    d <- vapply(1:300, function(s) {
      set.seed(s); n <- sample(10:80, 1); x <- lev + sample(0:5, n, TRUE) / div; y <- lev + sample(0:5, n, TRUE) / div
      o <- sample(n)
      a <- unlist(NNS.dep(x, y))
      max(abs(a - unlist(NNS.dep(y, x))), abs(a - unlist(NNS.dep(x[o], y[o]))))
    }, numeric(1))
    expect_lt(max(d), 1e-12)
  }
  # Decimal ties give the result of their exact integer image.
  set.seed(208); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
  expect_equal(NNS.dep(x / 10, y / 10), NNS.dep(x, y), tolerance = 1e-12)   # 0.167 vs -0.189 before
  set.seed(78); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
  # 1e9 + 0.1 is stored to 1.2e-7, so the shifted data matches its image to about 1e-6.
  expect_equal(NNS.dep(1e9 + x / 10, 1e9 + y / 10), NNS.dep(x, y), tolerance = 1e-6)   # 0.110 vs 0.191
  # Distinct values at a high level are not ties.
  ts <- 1.7e9 + seq_len(100) / 1e6
  expect_equal(unname(unlist(NNS.dep(ts, -ts))), c(-1, 1))
  expect_equal(NNS.dep(iris[, 1], iris[, 2]), NNS.dep(iris[, 2], iris[, 1]), tolerance = 1e-12)
})

test_that("the local dependence NNS.reg's default order reads is unchanged on tied data", {
  set.seed(388); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
  # Pinned again after #117 corrected the gravity splits (0.054, order 1 before).
  expect_equal(NNS:::.nns_dep_local(x, y), 0.281843356364165, tolerance = 1e-12)
  expect_identical(NNS:::.nns_reg_default_order(x, y), 2L)
  expect_equal(NNS.reg(x, y, point.est = 4, plot = FALSE, residual.plot = FALSE)$Point.est,
               2.4009920634920636, tolerance = 1e-12)
})
