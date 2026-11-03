# NNS.dep reads the debiased cell dependence at every level of the recursive
# mean split, on x, on y and on the joint quadrants, and reports the strongest
# level (#95).

test_that("a strong but noisy predictor is not scored 0", {
  # Before #95 the iris signal scored 0 and the noise column 0.065.
  set.seed(3)
  X <- cbind(signal = iris$Sepal.Length, noise = rnorm(150))
  y <- iris$Petal.Width
  expect_gt(NNS.dep(X[, 1], y, asym = TRUE)$Dependence, 0.75)
  expect_lt(NNS.dep(X[, 2], y, asym = TRUE)$Dependence, 0.3)
  f <- NNS.reg(X, y, dim.red.method = "NNS.dep", plot = FALSE, residual.plot = FALSE)
  expect_gt(f$equation$Coefficient[1], 4 * f$equation$Coefficient[2])
  expect_gt(f$R2, 0.6)

  # y = x + e (r = 0.71) ranged from 0.03 to 0.32 across seeds.
  dep <- vapply(1:5, function(s) {
    set.seed(s); x <- rnorm(1000)
    NNS.dep(x, x + rnorm(1000))$Dependence
  }, numeric(1))
  expect_true(all(dep > 0.55 & dep < 0.8))
})

test_that("relations without a function are detected by the joint quadrants", {
  set.seed(123)
  df <- data.frame(x = runif(10000, -1, 1), y = runif(10000, -1, 1))
  df <- subset(df, x^2 + y^2 <= 1 & x^2 + y^2 >= 0.95)
  circle <- NNS.dep(df$x, df$y)
  expect_gt(circle$Dependence, 0.8)
  expect_lt(abs(circle$Correlation), 0.2)

  set.seed(5); x <- rnorm(1000)
  expect_gt(NNS.dep(x, x * rnorm(1000))$Dependence, 0.4)
})

test_that("asym = TRUE measures the dependence of y on x only", {
  set.seed(6); x <- rnorm(500); y <- x^2 + rnorm(500, sd = 0.1)
  expect_gt(NNS.dep(x, y, asym = TRUE)$Dependence, 0.8)
  expect_lt(NNS.dep(y, x, asym = TRUE)$Dependence, 0.3)
  # The symmetric measure includes both conditional partitions.
  sym <- NNS.dep(x, y)$Dependence
  expect_gte(sym, NNS.dep(x, y, asym = TRUE)$Dependence)
  expect_gte(sym, NNS.dep(y, x, asym = TRUE)$Dependence)
})

test_that("discrete margins: no blend floor, and |Correlation| <= Dependence", {
  g <- expand.grid(x = c(rep(0, 9), 1), y = c(rep(0, 9), 1))   # exactly independent
  expect_equal(unname(unlist(NNS.dep(g$x, g$y))), c(0, 0))
  set.seed(4); x <- sample(0:4, 100, TRUE); y <- c(0, 3, 0, 4, 3)[x + 1]
  for (asym in c(FALSE, TRUE)) {
    r <- NNS.dep(x, y, asym = asym)
    expect_equal(r$Dependence, 1)
    expect_lte(abs(r$Correlation), r$Dependence)
  }
  b <- rep(0:1, each = 50)
  expect_equal(unname(unlist(NNS.dep(b, 1 - b))), c(-1, 1))
})

test_that("small tied cells use the exact permutation null", {
  # Seven points form a single cell. x has six tied values, so the null of the
  # degree-1 term is exact over the 7 distinct arrangements (the normal
  # approximation used for larger cells gives 0.659 here instead of 0.483).
  x <- c(rep(0, 6), 1); y <- c(1, 1, 1, 2, 2, 3, 9)
  dx <- x - mean(x); dy <- y - mean(y)
  t1 <- function(v) abs(sum(v * dy)) / sum(abs(v * dy))
  e1 <- mean(vapply(1:7, function(k) t1(replace(rep(dx[1], 7), k, dx[7])), numeric(1)))
  a <- sum(dx <= 0); b <- sum(dy <= 0); conc <- sum((dx <= 0) == (dy <= 0))
  k <- max(0, a + b - 7):min(a, b)
  e0 <- sum(dhyper(k, a, 7 - a, b) * abs(2 * (7 - a - b + 2 * k) / 7 - 1))
  t0 <- abs(2 * conc / 7 - 1)
  expected <- max(0, ((t0 - e0) / (1 - e0) + (t1(dx) - e1) / (1 - e1)) / 2)
  expect_equal(e1, 0.48295, tolerance = 1e-5)
  expect_equal(NNS.dep(x, y)$Dependence, expected, tolerance = 1e-12)
  expect_equal(NNS.dep(x, y, asym = TRUE)$Dependence, expected, tolerance = 1e-12)
})

test_that("the native and partition-map paths agree", {
  set.seed(9); x <- rnorm(300); y <- x^2 + rnorm(300)
  native <- c(unlist(NNS.dep(x, y)), unlist(NNS.dep(x, y, asym = TRUE)))
  old <- options(NNS.native = FALSE); on.exit(options(old))
  expect_identical(c(unlist(NNS.dep(x, y)), unlist(NNS.dep(x, y, asym = TRUE))), native)
})

test_that("NNS.reg's order heuristic keeps the local (finest-cell) dependence", {
  # The default order reads dependence as a fine-scale signal-to-noise
  # measure. The multi-level NNS.dep also counts dependence between coarse
  # cells (0.64 for y = x + e below) and would over-partition noisy data, so
  # the heuristic keeps the finest-cell measure. Its copula component now
  # averages both response orientations, preserving the sign-invariance
  # expected of a smoothing-complexity choice.
  set.seed(1); x <- rnorm(500); y <- x + rnorm(500)
  expect_equal(NNS:::.nns_dep_local(x, y), 0)
  expect_gt(NNS.dep(x, y, asym = TRUE)$Dependence, 0.55)
  expect_equal(NNS:::.nns_reg_default_order(x, y), NNS:::.nns_reg_default_order(x, -y))
  set.seed(2); x <- runif(400); y <- sin(6 * x) + rnorm(400, sd = .2)
  expect_equal(NNS:::.nns_dep_local(x, y), 0.394552637766, tolerance = 1e-10)
  expect_equal(NNS:::.nns_reg_default_order(x, y), NNS:::.nns_reg_default_order(x, -y))
  set.seed(3); u <- sample(1:4, 200, TRUE); v <- u^2 + rnorm(200, sd = .5)
  expect_equal(NNS:::.nns_dep_local(u, v, asym = FALSE), 0.0905745380689, tolerance = 1e-10)
})
