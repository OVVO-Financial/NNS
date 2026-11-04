# The default neighbour count of multivariate NNS.reg counts regression
# points, so it scales with their number rather than with n (issue #153).
# The point count stops growing with n at a fixed order; a sqrt(n) count then
# averaged over nearly every point and pulled estimates toward the mean.

g2 <- function(n) {
  x <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  list(x = x, p = plogis(1.5 * x$x1 - x$x2))
}

test_that("default class probabilities do not flatten as n grows", {
  set.seed(99)
  te <- g2(2000)
  for (n in c(2000, 5000)) {
    set.seed(n)
    tr <- g2(n)
    y <- rbinom(n, 1, tr$p)
    r <- NNS.reg(tr$x, y, point.est = te$x, plot = FALSE, ncores = 1)
    cells <- nrow(r$RPM)
    expect_lte(r$n.best, floor(sqrt(cells)))
    # sqrt(n) neighbours gave 0.17 here; k = 3 gives 0.04.
    expect_lt(mean(abs(r$Point.probabilities[, "1"] - te$p)), 0.07)
  }
})

test_that("default multivariate regression keeps local structure as n grows", {
  set.seed(99)
  te <- g2(2000)
  set.seed(2000)
  tr <- g2(2000)
  y <- tr$p + rnorm(2000, sd = 0.2)
  r <- NNS.reg(tr$x, y, point.est = te$x, plot = FALSE, ncores = 1, type = "REG")
  expect_lte(r$n.best, floor(sqrt(nrow(r$RPM))))
  # sqrt(n) neighbours gave 0.10 here.
  expect_lt(mean(abs(r$Point.est - te$p)), 0.06)
})
