test_that("NNS.caus p.value statistic is oriented x -> y", {
  set.seed(1); x <- rnorm(100); y <- x^2 + rnorm(100, sd = 0.1)
  xy <- NNS.caus(x, y, p.value = TRUE, nperm = 10, seed = 2)
  yx <- NNS.caus(y, x, p.value = TRUE, nperm = 10, seed = 2)
  expect_identical(names(yx$causation), c("C(x--->y)", "C(y--->x)"))
  # The statistic is C(x--->y) - C(y--->x): y = x^2 is driven by x.
  expect_gt(xy$p.value$observed_signed, 0.5)
  expect_equal(yx$p.value$observed_signed, -xy$p.value$observed_signed)
})

test_that("NNS.dep p.value without print.map opens no graphics device", {
  graphics.off()
  set.seed(11); x <- rnorm(50); y <- rnorm(50)
  expect_silent(NNS.dep(x, y, p.value = TRUE))
  expect_null(dev.list())
})

test_that("NNS.dep Dependence p.value is the upper-tail permutation p-value", {
  set.seed(5); x <- rnorm(50); y <- x + rnorm(50)
  obs <- NNS.dep(x, y)$Dependence
  set.seed(6); perms <- replicate(100L, sample.int(50))
  null <- apply(perms, 2, function(i) NNS.dep(x, y[i])$Dependence)
  set.seed(6); p <- NNS.dep(x, y, p.value = TRUE)
  expect_equal(p$`Dependence p.value`, (1 + sum(null >= obs)) / 101)
  # 7 / 101 before the per-cell debiasing (#70) removed the independence floor.
  expect_equal(p$`Dependence p.value`, 1 / 101)
})

test_that("NNS.dep Dependence p.value is at least 1 / 101", {
  set.seed(7); x <- rnorm(50)
  expect_equal(NNS.dep(x, x, p.value = TRUE)$`Dependence p.value`, 1 / 101)
})

test_that("NNS.dep Correlation p.value is the two-sided permutation p-value", {
  set.seed(5); x <- rnorm(50); y <- x + rnorm(50)
  obs <- NNS.dep(x, y)$Correlation
  set.seed(6); perms <- replicate(100L, sample.int(50))
  null <- apply(perms, 2, function(i) NNS.dep(x, y[i])$Correlation)
  set.seed(6); p <- NNS.dep(x, y, p.value = TRUE)
  expect_equal(p$`Correlation p.value`,
               min(1, 2 * min(1 + sum(null <= obs), 1 + sum(null >= obs)) / 101))
  expect_equal(p$`Correlation p.value`, 2 / 101)
})

test_that("NNS.dep Correlation p.value counts tied null correlations", {
  # 99 of the 100 positions of y's event give Correlation 0; it was 0 before.
  x <- c(1, rep(0, 99)); y <- c(0, 1, rep(0, 98))
  set.seed(871); p <- NNS.dep(x, y, p.value = TRUE)
  expect_equal(p$Correlation, 0)
  expect_equal(p$`Correlation p.value`, 1)
})

test_that("NNS.dep Correlation p.value holds its size under independence", {
  # 28 / 200 rejections at 0.05 when it was the smaller one-sided tail.
  set.seed(5)
  p <- replicate(200L, NNS.dep(rnorm(30), rnorm(30), p.value = TRUE)$`Correlation p.value`)
  expect_lte(mean(p <= 0.05), 0.10)
})

test_that("NNS.dep p.value = TRUE honours asym", {
  set.seed(3); x <- rnorm(100); y <- x^2 + rnorm(100, sd = 0.1)
  direct <- NNS.dep(y, x, asym = TRUE)
  set.seed(4); p <- NNS.dep(y, x, p.value = TRUE, asym = TRUE)
  expect_equal(p$Correlation, direct$Correlation)
  expect_equal(p$Dependence, direct$Dependence)
  # x is not a function of y = x^2: 0.6005389 before the per-cell debiasing
  # (#70), 0.1027830 before the multi-level measure (#95).
  expect_equal(p$Dependence, 0.1156500, tolerance = 1e-6)
})

test_that("NNS.dep rejects infinite values", {
  set.seed(8); x <- rnorm(30); y <- x + rnorm(30)
  expect_error(NNS.dep(replace(x, 3, Inf), y), "x must be finite")
  expect_error(NNS.dep(replace(x, 3, -Inf), y), "x must be finite")
  expect_error(NNS.dep(x, replace(y, 5, Inf)), "y must be finite")
  expect_error(NNS.dep(cbind(replace(x, 3, Inf), y)), "x must be finite")
  expect_error(NNS.dep(replace(x, 3, NaN), y), "missing values")
})

test_that("NNS.caus treats single-column data frames like vectors", {
  x <- as.double(1:20); y <- x^2
  v <- NNS.caus(x, y)
  expect_equal(NNS.caus(data.frame(x), y), v)
  expect_equal(NNS.caus(x, data.frame(y)), v)
  expect_equal(NNS.caus(data.frame(x), data.frame(y)), v)
  expect_error(NNS.caus(data.frame(x = replace(x, 2, Inf)), y), "must be finite")
})

test_that("NNS.caus keeps a factor in a single-column data frame", {
  f <- factor(rep(c("a", "b"), 30))
  y <- as.double(1:60)^2
  v <- NNS.caus(f, y, factor.2.dummy = TRUE)
  expect_equal(unname(v[["C(x--->y)"]] + v[["C(y--->x)"]]), 1)
  expect_equal(NNS.caus(data.frame(f), y, factor.2.dummy = TRUE), v)
})

test_that("NNS.caus rejects multi-column x or y when y is supplied", {
  set.seed(1); x <- rnorm(60); y <- x^2 + rnorm(60, sd = 0.1)
  expect_error(NNS.caus(data.frame(x, x), y), "single variable")
  expect_error(NNS.caus(x, cbind(y, y)), "single variable")
  expect_error(NNS.caus(factor(rep(1:3, 20)), y, factor.2.dummy = TRUE), "single variable")
  expect_length(NNS.caus(factor(rep(1:2, 30)), y, factor.2.dummy = TRUE), 2)
})

test_that("NNS.caus rejects unequal lengths", {
  set.seed(1); x <- rnorm(60); y <- rnorm(60)
  expect_error(NNS.caus(x, y[-1]), "same length")
})

test_that("NNS.caus rejects infinite values", {
  set.seed(1); x <- rnorm(60); y <- rnorm(60)
  expect_error(NNS.caus(replace(x, 4, Inf), y), "must be finite")
  expect_error(NNS.caus(x, replace(y, 4, -Inf)), "must be finite")
})

test_that("NNS.caus rejects constant inputs", {
  set.seed(1); x <- rnorm(60)
  expect_error(NNS.caus(x, rep(1, 60)), "must not be constant")
  expect_error(NNS.caus(rep(1, 60), x, tau = 3), "must not be constant")
})

test_that("NNS.caus rejects tau that is not \"cs\", \"ts\" or a non-negative integer", {
  set.seed(1); x <- rnorm(60); y <- x^2 + rnorm(60, sd = 0.1)
  for (tau in list(1.5, -1, "3", "foo", c(1, 2), NA, TRUE))
    expect_error(NNS.caus(x, y, tau = tau), "tau must be")
  expect_equal(NNS.caus(x, y, tau = 3L), NNS.caus(x, y, tau = 3))
  expect_equal(NNS.caus(x, y, tau = "cs"), NNS.caus(x, y, tau = 0))
})

test_that("NNS.caus rejects tau that leaves fewer than three lagged rows", {
  set.seed(1); x <- rnorm(20); y <- rnorm(20)
  expect_length(NNS.caus(x, y, tau = 17), 2)
  expect_error(NNS.caus(x, y, tau = 18), "tau must be")
  expect_error(NNS.caus(c(1, 3), c(2, 1)), "tau must be")
})

test_that("NNS.caus tau = \"ts\" stops when no seasonal period qualifies", {
  ar1 <- function(seed) {
    set.seed(seed); e <- rnorm(50); a <- numeric(50)
    for (i in 2:50) a[i] <- 0.9 * a[i - 1] + e[i]
    list(a = a, b = rnorm(50))
  }
  # NNS.seas reports only lags with >= 5 points per component (#69), so an
  # AR(1) usually returns Period 1 (no seasonality), which qualifies as
  # tau = 1. For seed 6 every period NNS.seas reports for x exceeds sqrt(50).
  s <- ar1(6)
  expect_true(all(NNS.seas(s$a, plot = FALSE)$periods > sqrt(50)))
  expect_error(NNS.caus(s$a, s$b, tau = "ts"), "no seasonal period")
  s <- ar1(8)
  expect_length(NNS.caus(s$a, s$b, tau = "ts"), 2)
})
