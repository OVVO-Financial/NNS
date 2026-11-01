anova <- function(..., seed = 1) {
  set.seed(seed)
  tryCatch(NNS.ANOVA(..., plot = FALSE), error = conditionMessage)
}

x_na <- c(-2, NA, 0, NaN, 1, 3)
y_na <- c(NA, -1, 1, 2)
x_ok <- c(-2, 0, 1, 3)
y_ok <- c(-1, 1, 2)

test_that("two-sample NNS.ANOVA omits NA and NaN", {
  expect_equal(anova(x_na, y_na), anova(x_ok, y_ok))
})

test_that("two-sample NNS.ANOVA omits NA and NaN with medians = TRUE", {
  expect_equal(anova(x_na, y_na, medians = TRUE), anova(x_ok, y_ok, medians = TRUE))
})

test_that("two-sample NNS.ANOVA omits NA and NaN with robust = TRUE", {
  expect_equal(anova(x_na, y_na, robust = TRUE), anova(x_ok, y_ok, robust = TRUE))
})

x <- c(-2, -1, 0, 2, 6, 8, 10)
y <- c(-3, 1, 4, 5)

test_that("robust mode bootstraps each group independently from its full sample", {
  set.seed(314)
  cx <- replicate(100, x[sample.int(length(x), replace = TRUE)])
  ty <- replicate(100, y[sample.int(length(y), replace = TRUE)])
  cs <- vapply(seq_len(100), function(i) {
    NNS:::NNS.ANOVA.bin(cx[, i], ty[, i], plot = FALSE)$Certainty
  }, numeric(1))
  expected <- c(NNS.gravity(cs), unname(quantile(cs, c(0.025, 0.975))))

  r <- anova(x, y, robust = TRUE, seed = 314)
  expect_type(r, "list")
  expect_equal(
    unname(unlist(r[c("Robust Certainty Estimate", "Lower 95% CI", "Upper 95% CI")])),
    expected
  )
})

test_that("robust mode resamples past the first min(n1, n2) values of the longer group", {
  set.seed(1)
  long <- c(rnorm(10), rnorm(90, 5))
  short <- rnorm(10)
  r <- anova(long, short, robust = TRUE, seed = 2)
  expect_equal(r$Certainty, 0)
  expect_equal(r[["Robust Certainty Estimate"]], 0)
  expect_equal(r[["Upper 95% CI"]], 0)
})

boot_stats <- function(seed) {
  set.seed(seed)
  cb <- matrix(sample(x, 1000 * length(x), replace = TRUE), nrow = length(x))
  tb <- matrix(sample(y, 1000 * length(y), replace = TRUE), nrow = length(y))
  list(c = colMeans(cb), t = colMeans(tb))
}

bounds <- function(r) {
  if (is.list(r)) unname(c(r$Effect_Size_LB, r$Effect_Size_UB)) else r
}

test_that("tails = 'left' gives a finite lower effect bound and no upper bound", {
  b <- boot_stats(32)
  expect_equal(
    bounds(anova(x, y, tails = "left", seed = 32)),
    c(unname(quantile(b$t, 0.05) - quantile(b$c, 0.95)), Inf)
  )
})

test_that("tails = 'right' gives a finite upper effect bound and no lower bound", {
  b <- boot_stats(32)
  expect_equal(
    bounds(anova(x, y, tails = "right", seed = 32)),
    c(-Inf, unname(quantile(b$t, 0.95) - quantile(b$c, 0.05)))
  )
})

test_that("one-sided tails also work with robust = TRUE", {
  left <- bounds(anova(x, y, tails = "left", robust = TRUE))
  right <- bounds(anova(x, y, tails = "right", robust = TRUE))
  expect_true(is.numeric(left) && is.finite(left[1]) && left[2] == Inf)
  expect_true(is.numeric(right) && right[1] == -Inf && is.finite(right[2]))
})
