arma_pi <- function(x, ...) {
  set.seed(123)
  NNS.ARMA(x, plot = FALSE, pred.int = 0.95, h = 6, ...)
}

test_that("lower prediction bound sits below the forecast on a negative-valued series", {
  set.seed(1)
  x <- sin(1:200 / 2) * 5 - 10 + rnorm(200, sd = 0.5)
  r <- arma_pi(x, seasonal.factor = 13)
  expect_true(all(r[[1]] < 0))
  expect_true(all(r[[2]] < r[[1]]))
  expect_true(all(r[[1]] < r[[3]]))
})

test_that("lower prediction bound sits below the forecast on a mean-zero series", {
  set.seed(2)
  x <- rnorm(250, 0, 0.01)
  r <- arma_pi(x, seasonal.factor = FALSE)
  expect_true(all(r[[2]] < r[[1]]))
  expect_true(all(r[[1]] < r[[3]]))
})
