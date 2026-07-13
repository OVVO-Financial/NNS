# REVISION: 2026-07-13-NNS-REGRESSION-AUDIT-TESTS

test_that("NNS.reg preserves warning settings on success and failure", {
  old <- getOption("warn")
  on.exit(options(warn = old), add = TRUE)
  options(warn = 1)

  expect_silent(NNS.reg(1:20, (1:20)^2, plot = FALSE))
  expect_equal(getOption("warn"), 1)

  expect_error(NNS.reg(1:20, 1:10, plot = FALSE), "rows")
  expect_equal(getOption("warn"), 1)
})

test_that("training and response dimensions are checked before fitting", {
  expect_error(
    NNS.reg(matrix(1:20, ncol = 2), 1:5, plot = FALSE),
    "rows"
  )
})

test_that("named prediction columns are reordered and schema checked", {
  set.seed(1)
  x <- data.frame(a = rnorm(40), b = rnorm(40))
  y <- x$a - 2 * x$b
  p1 <- x[1:5, c("a", "b")]
  p2 <- x[1:5, c("b", "a")]

  f1 <- NNS.reg(x, y, point.est = p1, plot = FALSE,
                residual.plot = FALSE, ncores = 1)
  f2 <- NNS.reg(x, y, point.est = p2, plot = FALSE,
                residual.plot = FALSE, ncores = 1)
  expect_equal(f1$Point.est, f2$Point.est, tolerance = 1e-12)

  bad <- p1
  names(bad) <- c("a", "c")
  expect_error(
    NNS.reg(x, y, point.est = bad, plot = FALSE,
            residual.plot = FALSE),
    "exactly match"
  )
})

test_that("dimension reduction is invariant to unrelated prediction rows", {
  set.seed(2)
  x <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
  y <- x$a + x$b^2 - x$c
  base <- x[1:3, ]
  extra <- rbind(base, data.frame(a = 1e6, b = -1e6, c = 1e6))

  f1 <- NNS.reg(x, y, point.est = base, dim.red.method = "cor",
                plot = FALSE, residual.plot = FALSE)
  f2 <- NNS.reg(x, y, point.est = extra, dim.red.method = "cor",
                plot = FALSE, residual.plot = FALSE)
  expect_equal(f1$Point.est, f2$Point.est[1:3], tolerance = 1e-12)
})

test_that("factor encoding is training-fitted and unseen levels fail", {
  x <- data.frame(group = factor(rep(c("a", "b"), each = 15)),
                  z = seq_len(30))
  y <- rep(c(1, 2), each = 15)
  ok <- data.frame(group = factor(c("a", "b"), levels = c("a", "b")),
                   z = c(4, 20))
  expect_silent(
    NNS.reg(x, y, point.est = ok, type = "CLASS",
            plot = FALSE, residual.plot = FALSE)
  )

  bad <- data.frame(group = "c", z = 10)
  expect_error(
    NNS.reg(x, y, point.est = bad, type = "CLASS",
            plot = FALSE, residual.plot = FALSE),
    "unseen level"
  )
})

test_that("classification returns numeric observed class values", {
  train <- iris[-(141:150), ]
  test <- iris[141:150, 1:4]
  fit <- NNS.reg(train[, 1:4], train[, 5], point.est = test,
                 type = "CLASS", plot = FALSE,
                 residual.plot = FALSE, ncores = 1)
  expect_true(is.numeric(fit$Point.est))
  expect_true(all(fit$Point.est %in% 1:3))
  expect_identical(fit$class.levels, levels(iris$Species))
})

test_that("count responses remain regression unless CLASS is explicit", {
  x <- seq_len(60)
  y <- rep(1:3, each = 20)
  regression <- NNS.reg(x, y, plot = FALSE)
  expect_null(regression$Prediction.Accuracy)

  classification <- NNS.reg(x, y, type = "CLASS", plot = FALSE)
  expect_true(is.numeric(classification$Prediction.Accuracy))
})

test_that("XONLY does not activate classification", {
  x <- seq_len(40)
  y <- sin(x / 4)
  fit <- NNS.reg(x, y, type = "XONLY", plot = FALSE)
  expect_null(fit$Prediction.Accuracy)
})

test_that("numeric dimension-reduction equation reproduces X star", {
  set.seed(3)
  x <- cbind(a = rnorm(40), b = rnorm(40), c = rnorm(40))
  y <- rowSums(x)
  coef <- c(2, -1, 0.5)
  fit <- NNS.reg(x, y, dim.red.method = coef, plot = FALSE)

  mins <- apply(x, 2, min)
  ranges <- apply(x, 2, max) - mins
  norm <- matrix(0.5, nrow(x), ncol(x))
  active <- ranges > 0
  norm[, active] <- sweep(sweep(x[, active, drop = FALSE], 2,
                                mins[active], "-"), 2, ranges[active], "/")
  expected <- as.numeric(norm %*% coef / sum(abs(coef) > 0))
  expect_equal(as.numeric(fit$x.star$x), expected, tolerance = 1e-12)
  expect_equal(tail(fit$equation$Coefficient, 1), sum(abs(coef) > 0))
})

test_that("all-zero dimension weights fall back to a defined model", {
  x <- cbind(a = rep(1, 30), b = rep(2, 30))
  y <- seq_len(30)
  fit <- NNS.reg(x, y, dim.red.method = c(0, 0), plot = FALSE)
  expect_true(all(is.finite(fit$x.star$x)))
  expect_gt(tail(fit$equation$Coefficient, 1), 0)
})

test_that("order max fitted values use the same prediction function", {
  x <- c(1, 1, 2, 3, 4, 5)
  y <- c(1, 3, 4, 9, 16, 25)
  fit <- NNS.reg(x, y, order = "max", point.est = x, plot = FALSE)
  expect_equal(fit$Fitted.xy$y.hat, fit$Point.est, tolerance = 1e-12)
})

test_that("one-row matrix and vector predictions are identical", {
  set.seed(4)
  x <- cbind(a = rnorm(50), b = rnorm(50))
  y <- x[, 1] - x[, 2]
  v <- unname(x[5, ])
  m <- matrix(v, nrow = 1)

  fv <- NNS.reg(x, y, point.est = v, plot = FALSE,
                residual.plot = FALSE, ncores = 1)
  fm <- NNS.reg(x, y, point.est = m, plot = FALSE,
                residual.plot = FALSE, ncores = 1)
  expect_equal(fv$Point.est, fm$Point.est, tolerance = 1e-12)
})

test_that("n.best one fitted values equal predictions on training rows", {
  set.seed(5)
  x <- cbind(a = rnorm(40), b = rnorm(40))
  y <- x[, 1]^2 + x[, 2]
  fit <- NNS.reg(x, y, point.est = x, n.best = 1,
                 plot = FALSE, residual.plot = FALSE, ncores = 1)
  expect_equal(fit$Fitted.xy$y.hat, fit$Point.est, tolerance = 1e-12)
})

test_that("distance modes are validated and recorded", {
  set.seed(6)
  x <- cbind(a = rnorm(30), b = rnorm(30))
  y <- x[, 1] + x[, 2]
  for (d in c("L1", "L2", "FACTOR")) {
    fit <- NNS.reg(x, y, point.est = x[1:2, ], dist = d,
                   plot = FALSE, residual.plot = FALSE)
    expect_identical(fit$dist, d)
  }
  expect_error(NNS.reg(x, y, dist = "DTW", plot = FALSE), "dist")
})

test_that("prediction intervals preserve prediction row count", {
  x <- seq_len(40)
  y <- x + sin(x)
  points <- c(-10, 1, 20, 50)
  fit <- NNS.reg(x, y, point.est = points,
                 confidence.interval = 0.9, plot = FALSE)
  expect_equal(nrow(fit$pred.int), length(points))
  expect_true(all(fit$pred.int$pred.int.neg <= fit$pred.int$pred.int.pos))
})

test_that("missing and nonfinite prediction values fail", {
  x <- cbind(a = 1:20, b = 21:40)
  y <- 1:20
  expect_error(
    NNS.reg(x, y, point.est = matrix(c(NA, 2), nrow = 1),
            plot = FALSE, residual.plot = FALSE),
    "finite"
  )
})

test_that("NNS.part order max never passes a character value to C++", {
  fit <- NNS.part(c(1, 1, 2, 3), c(1, 3, 4, 9),
                  order = "max", type = "XONLY")
  expect_equal(nrow(fit$regression.points), 3)
  expect_identical(fit$order, 3L)
})
