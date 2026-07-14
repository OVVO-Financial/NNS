test_that("Racine-Hastie R2 is the squared observed-fitted correlation", {
  actual <- c(-2, -1, 0, 1, 2)
  predicted <- c(8, 4, 0, -4, -8)

  expected <- stats::cor(actual, predicted)^2
  expect_equal(NNS:::.nns_reg_r2(actual, predicted), expected, tolerance = 1e-14)
  expect_equal(expected, 1)
})


test_that("Racine-Hastie R2 remains bounded when predictive R2 is negative", {
  actual <- c(-2, -1, 0, 1, 2)
  predicted <- c(20, 10, 0, -10, -20)

  predictive_r2 <- 1 - sum((actual - predicted)^2) /
    sum((actual - mean(actual))^2)
  nns_r2 <- NNS:::.nns_reg_r2(actual, predicted)

  expect_lt(predictive_r2, 0)
  expect_gte(nns_r2, 0)
  expect_lte(nns_r2, 1)
  expect_equal(nns_r2, stats::cor(actual, predicted)^2, tolerance = 1e-14)
})


test_that("Racine-Hastie R2 handles constant-series degeneracy explicitly", {
  expect_equal(NNS:::.nns_reg_r2(rep(3, 5), rep(3, 5)), 1)
  expect_equal(NNS:::.nns_reg_r2(rep(3, 5), rep(4, 5)), 0)
  expect_equal(NNS:::.nns_reg_r2(1:5, rep(3, 5)), 0)
})


test_that("reported NNS.reg and NNS.M.reg R2 values are bounded", {
  x <- seq(-2, 2, length.out = 40)
  y <- sin(4 * x) + x^2
  univariate <- NNS.reg(x, y, order = 1, plot = FALSE, residual.plot = FALSE)

  X <- cbind(x, x^2)
  multivariate <- NNS.M.reg(
    X, y, order = 1, n.best = 1,
    plot = FALSE, residual.plot = FALSE
  )

  expect_gte(univariate$R2, 0)
  expect_lte(univariate$R2, 1)
  expect_gte(multivariate$R2, 0)
  expect_lte(multivariate$R2, 1)
})
