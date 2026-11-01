# A product grid is exactly independent, so the copula is 0.
grid_X <- function() {
  as.matrix(expand.grid(c(0.1, 0.2, 0.3, 0.5, 1.7, 4.1),
                        c(1, 1.5, 2, 9),
                        c(0.05, 0.4, 0.45, 3.3, 7)))
}

test_that("NNS.copula is 0 on independent skewed margins at the default target", {
  X <- grid_X()
  expect_equal(NNS.copula(X), 0, tolerance = 1e-6)
  expect_equal(NNS.copula(X, continuous = FALSE), 0, tolerance = 1e-6)
  expect_equal(NNS.copula(X[, 1:2]), 0, tolerance = 1e-6)
  x <- X[, 1]
  expect_equal(NNS.copula(cbind(x, x, -x)), 1)
})

test_that("NNS.copula is 0 on independent margins at a non-mean target", {
  X <- grid_X()
  t <- c(0.25, 1.75, 0.42)
  expect_equal(NNS.copula(X, target = t), 0, tolerance = 1e-6)
  expect_equal(NNS.copula(X, target = t, continuous = FALSE), 0, tolerance = 1e-6)
  Y <- rbind(c(-10, 0.1, 0.1), c(0.1, -10, 0.1), c(0.1, 0.1, -10))
  expect_lte(NNS.copula(Y, target = rep(0, 3)), 1)
})

test_that("NNS.copula drops constant margins", {
  expect_equal(NNS.copula(cbind(-2:2, 1)), 0)
  expect_equal(NNS.copula(cbind(-2:2, 1), continuous = FALSE), 0)
  x <- 1:5
  expect_equal(NNS.copula(cbind(x, x, 3)), 1)
})

test_that("A constant margin leaves the copula of the varying margins unchanged", {
  # p = 0.5 in both margins, where the old fixed anchors were exact; values match the old cbind(x, y)
  x <- c(1, 2, 3, 4, 5, 6)
  y <- c(3, 1, 6, 2, 5, 4)
  expect_equal(NNS.copula(cbind(x, y, 5)), 0.690065559342, tolerance = 1e-10)
  expect_equal(NNS.copula(cbind(x, y, 5), continuous = FALSE), sqrt(1/3), tolerance = 1e-10)
  Z <- iris[, 1:3]
  expect_equal(NNS.copula(cbind(Z, k = 1)), NNS.copula(Z))
  expect_equal(NNS.copula(cbind(Z, k = 1), target = c(5, 3, 4, 0)), NNS.copula(Z, target = c(5, 3, 4)))
})
