# NNS.dep measures each cell's concordance against the within-cell
# permutation null, so independence no longer leaves a floor (#70).

test_that("NNS.dep has no independence floor", {
  set.seed(70)
  dep <- replicate(20, NNS.dep(rnorm(1000), rnorm(1000))$Dependence)
  # Before the debiasing the mean was about 0.32 at n = 1000.
  expect_lt(mean(dep), 0.06)
  set.seed(71)
  dep100 <- replicate(40, NNS.dep(rnorm(100), rnorm(100))$Dependence)
  # About 0.60 before.
  expect_lt(mean(dep100), 0.15)
})

test_that("NNS.dep keeps perfect dependence at 1", {
  set.seed(72); x <- rnorm(500)
  expect_equal(unname(unlist(NNS.dep(x, x))), c(1, 1))
  expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1))
  b <- c(rep(0, 51), rep(1, 49))
  expect_equal(unname(unlist(NNS.dep(b, b))), c(1, 1))
})

test_that("NNS.dep Dependence stays in [0, 1] and bounds |Correlation|", {
  set.seed(73)
  for (i in 1:20) {
    n <- sample(c(20, 60, 200), 1)
    x <- rnorm(n); y <- 0.5 * x + rnorm(n)
    for (asym in c(FALSE, TRUE)) {
      r <- NNS.dep(x, y, asym = asym)
      expect_gte(r$Dependence, 0)
      expect_lte(r$Dependence, 1)
      expect_lte(abs(r$Correlation), r$Dependence + 1e-12)
    }
  }
})

test_that("NNS.dep separates weak dependence from independence", {
  set.seed(74)
  pairs <- replicate(10, {
    x <- rnorm(1000)
    c(NNS.dep(x, x + rnorm(1000))$Dependence, NNS.dep(x, rnorm(1000))$Dependence)
  })
  expect_gt(mean(pairs[1, ]), mean(pairs[2, ]) + 0.1)
  # Nonlinear dependence is still detected.
  x <- rnorm(1000)
  expect_gt(NNS.dep(x, x^2 + rnorm(1000, sd = 0.1))$Dependence, 0.6)
})

test_that("a cell with constant x counts only when y is determined", {
  set.seed(75)
  xb <- rbinom(500, 1, 0.5); yb <- rbinom(500, 1, 0.5)
  expect_lt(NNS.dep(xb, yb)$Dependence, 0.1)
  expect_equal(NNS.dep(xb, xb)$Dependence, 1)
})
