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

test_that("NNS.SS rejects factors and malformed numbers in either sample", {
  expect_error(NNS.SS(factor("100"), 10), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(10, factor("100")), "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(c("10", "N/A"), 10), "`x` must contain only numeric values.", fixed = TRUE)
  expect_error(NNS.SS(10, c("10", "1,000"), confidence.interval = TRUE, reps = 9),
               "`y` must contain only numeric values.", fixed = TRUE)
})

test_that("NNS.SS accepts numeric strings and omits missing observations", {
  ref <- NNS.SS(c(10, 20), c(5, 10))
  expect_identical(NNS.SS(c("10", "20"), c("5", "10")), ref)
  expect_identical(NNS.SS(c(10, NA_real_, 20), c(NaN, 5, 10)), ref)
})

test_that("NNS.SS rejects factor columns and list elements before comparing samples", {
  fund <- data.frame(amount = factor("100"))
  returns <- c(5, 10, 20)
  expect_error(NNS.SS(fund, returns), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(returns, fund), "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(fund, returns, confidence.interval = TRUE, reps = 5),
               "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(returns, fund, confidence.interval = TRUE, reps = 5),
               "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(list(10, factor("100")), returns), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_identical(NNS.SS(data.frame(amount = 100), returns), NNS.SS(100, returns))
  expect_identical(NNS.SS(list(10, 100), returns), NNS.SS(c(10, 100), returns))
})
