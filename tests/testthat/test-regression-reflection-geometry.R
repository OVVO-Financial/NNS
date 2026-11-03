test_that("response reflection preserves default mean-regression geometry", {
  set.seed(128)
  x <- data.frame(a = rnorm(140), b = rnorm(140))
  y <- x$a^2 + 2 * x$b + rnorm(140, sd = .3)
  fit <- function(response) NNS.reg(x, response, point.est = x[1:10, ],
    noise.reduction = "mean", plot = FALSE, ncores = 1)
  a <- fit(y); b <- fit(-y)
  expect_equal(a$n.best, b$n.best)
  expect_equal(a$rhs.partitions, b$rhs.partitions)
  expect_equal(a$Point.est, -b$Point.est, tolerance = 1e-10)
  expect_equal(a$R2, b$R2, tolerance = 1e-10)
})

test_that("Racine-Hastie fit is squared correlation rather than one minus relative SSE", {
  y <- seq(-3, 4, length.out = 30)
  prediction <- 2 * y + 5
  expect_equal(NNS:::.nns_reg_r2(y, prediction), 1)
  expect_lt(1 - sum((y - prediction)^2) / sum((y - mean(y))^2), 0)
})
