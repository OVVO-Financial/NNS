test_that("NNS.SS gives an interval when one group is constant", {
  x <- rep(1, 30)
  y <- rep(c(0, 1), 15)
  set.seed(1)
  r <- NNS.SS(x, y, confidence.interval = TRUE, reps = 99)
  expect_equal(r$p_star, 0.75)
  expect_true(all(is.finite(r$boot_vals)))
  expect_true(r$lower <= 0.75 && 0.75 <= r$upper)

  set.seed(1)
  r <- NNS.SS(y, x, confidence.interval = TRUE, reps = 99)
  expect_equal(r$p_star, 0.25)
  expect_true(r$lower <= 0.25 && 0.25 <= r$upper)
})

test_that("NNS.SS interval for two constant groups is the point estimate", {
  r <- NNS.SS(rep(2, 10), rep(1, 15), confidence.interval = TRUE, reps = 9)
  expect_equal(c(r$lower, r$upper), c(1, 1))
  r <- NNS.SS(rep(1, 10), rep(1, 15), confidence.interval = TRUE, reps = 9)
  expect_equal(c(r$lower, r$upper), c(0.5, 0.5))
})

test_that("NNS.SS bootstraps non-constant groups with NNS.meboot as before", {
  set.seed(42)
  x <- rnorm(40, 0.3)
  y <- rnorm(30)
  set.seed(1)
  r <- NNS.SS(x, y, confidence.interval = TRUE, reps = 49)
  set.seed(1)
  xb <- NNS.meboot(x, reps = 49, rho = 1)["replicates", ]$replicates
  yb <- NNS.meboot(y, reps = 49, rho = 1)["replicates", ]$replicates
  vals <- vapply(1:49, function(i) NNS.SS(xb[, i], yb[, i])$p_star, numeric(1))
  expect_identical(r$boot_vals, vals)
})
