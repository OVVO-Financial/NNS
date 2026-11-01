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
  x <- -1:1
  expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1) * sqrt(2 / 3))
})

test_that("NNS.dep discrete blend on an independent tied grid uses the same split", {
  G <- as.matrix(expand.grid(-1:1, -1:1))[rep(1:9, each = 100), ]
  expect_equal(unname(unlist(NNS.dep(G[, 1], G[, 2]))), c(0, 1 / 6))
  expect_equal(unname(unlist(NNS.dep(G[, 1], G[, 2], asym = TRUE))), c(0, 1 / 6))
})
