# NNS.caus follows book chapter 14 (#68): explicitly lagged pairs on the copula
# scale, lagged co-partial-moment dependence, and a [0, 1] index with 0.5 neutral.

test_that("the index is in [0, 1], sums to 1 and is neutral for x = y", {
  set.seed(1); x <- rnorm(300); y <- x^2 + rnorm(300, sd = 0.2)
  v <- NNS.caus(x, y)
  expect_named(v, c("C(x--->y)", "C(y--->x)"))
  expect_true(all(v >= 0 & v <= 1))
  expect_equal(unname(sum(v)), 1)
  expect_equal(unname(NNS.caus(x, x)), c(0.5, 0.5))
  expect_equal(unname(NNS.caus(x, x, tau = 2)), c(0.5, 0.5))
})

test_that("swapping the arguments swaps the directions", {
  set.seed(2); x <- rnorm(300); y <- sin(2 * x) + rnorm(300, sd = 0.2)
  expect_equal(unname(NNS.caus(y, x)), rev(unname(NNS.caus(x, y))))
  expect_equal(unname(NNS.caus(y, x, tau = 1)), rev(unname(NNS.caus(x, y, tau = 1))))
})

test_that("a non-invertible function is attributed to its argument", {
  set.seed(3); x <- rnorm(1000)
  expect_gt(NNS.caus(x, x^2)[["C(x--->y)"]], 0.9)
  expect_lt(NNS.caus(x^2, x)[["C(x--->y)"]], 0.1)
})

test_that("an invertible relationship has no cross-sectional direction", {
  set.seed(4); x <- rnorm(500)
  expect_equal(unname(NNS.caus(x, 2 * x + 1)), c(0.5, 0.5))
})

test_that("the result is invariant to shifts, rescaling and increasing transforms", {
  set.seed(5); x <- rnorm(500, 1); y <- x - 1.5 + rnorm(500, sd = 0.3)
  v <- NNS.caus(x, y)
  expect_equal(NNS.caus(x, y + 10), v)
  expect_equal(NNS.caus(3 * x - 7, exp(y)), v)
  expect_equal(NNS.caus(x, y, tau = 2), NNS.caus(x + 100, 0.01 * y, tau = 2))
})

test_that("tau lags x against y and identifies lagged causation", {
  set.seed(6)
  lagged <- replicate(10, {
    x <- rnorm(400); y <- c(0, 0.8 * x[-400]) + 0.3 * rnorm(400)
    NNS.caus(x, y, tau = 1)[["C(x--->y)"]]
  })
  expect_gt(min(lagged), 0.7)
  # A contemporaneous comparison cannot see a pure lag-1 effect.
  set.seed(6); x <- rnorm(400); y <- c(0, 0.8 * x[-400]) + 0.3 * rnorm(400)
  expect_lt(abs(NNS.caus(x, y, tau = 0)[["C(x--->y)"]] - 0.5), 0.2)
})

test_that("independent series stay near the neutral 0.5", {
  set.seed(7)
  v <- replicate(20, NNS.caus(rnorm(500), rnorm(500), tau = 1)[["C(x--->y)"]])
  expect_lt(abs(mean(v) - 0.5), 0.1)
  expect_lt(max(abs(v - 0.5)), 0.4)
})

test_that("the causation matrix holds C(row ---> column) with a 0 diagonal", {
  set.seed(8); a <- rnorm(300); b <- a^2 + rnorm(300, sd = 0.1); cc <- rnorm(300)
  m <- NNS.caus(cbind(a, b, cc))
  expect_equal(unname(diag(m)), c(0, 0, 0))
  expect_equal(m["a", "b"], NNS.caus(a, b)[["C(x--->y)"]])
  expect_equal(m["b", "a"], NNS.caus(a, b)[["C(y--->x)"]])
  expect_equal(m["a", "b"] + m["b", "a"], 1)
  expect_gt(m["a", "b"], 0.8)
})

test_that("NNS.caus weights still drive dimension reduction", {
  set.seed(9); X <- cbind(x1 = rnorm(200), x2 = rnorm(200)); y <- X[, 1]^2 + rnorm(200, sd = 0.1)
  fit <- NNS.reg(X, y, dim.red.method = "NNS.caus", plot = FALSE, ncores = 1)
  expect_true(all(is.finite(fit$Fitted.xy$y.hat)))
  w <- fit$equation$Coefficient
  expect_gt(abs(w[1]), abs(w[2]))
})

test_that("p.value permutes the rows of a single-column data frame", {
  set.seed(1); x <- rnorm(80); y <- x^2 + rnorm(80, sd = 0.2)
  v <- NNS.caus(x, y, p.value = TRUE, nperm = 20, seed = 3)
  d <- NNS.caus(x, data.frame(y = y), p.value = TRUE, nperm = 20, seed = 3)
  expect_equal(d$p.value$null_distribution, v$p.value$null_distribution)
  expect_gt(length(unique(d$p.value$null_distribution)), 1)
  expect_equal(NNS.caus(data.frame(x = x), y, p.value = TRUE, nperm = 20, seed = 3, permute = "x")$p.value$two.sided,
               NNS.caus(x, y, p.value = TRUE, nperm = 20, seed = 3, permute = "x")$p.value$two.sided)
})
