# The continuous NNS.mode returns data values on discrete-valued data, finds a
# mode at the edge of the data, and backs the regression "mode" reducer (#114).

test_that("discrete-valued data return their own values, end values included", {
  expect_identical(NNS.mode(c(1, 1, 1, 5, 5, 5, 9, 9, 9)), c(1, 5, 9))   # 1, 3.33, 8 before
  expect_identical(NNS.mode(rep(0:4, each = 50)), as.numeric(0:4))       # 0.885 1.897 2.909
  expect_identical(NNS.mode(rep(0:19, each = 400)), as.numeric(0:19))    # 18 centres, no 0 or 19
  # Tied modes: the one nearest the gravity, so never a value in a gap.
  expect_identical(NNS.mode(c(1, 1, 1, 5, 5, 5, 9, 9, 9), multi = FALSE), 5)
  expect_identical(NNS.mode(c(0, 0, 1, 1), multi = FALSE), 0)
  # A mode is a value more frequent than its neighbours: two clusters, two modes.
  x <- rep(c(1, 2, 3, 7, 8, 9), c(10, 30, 10, 10, 25, 10))
  expect_identical(NNS.mode(x), c(2, 8))
  expect_identical(NNS.mode(x, multi = FALSE), 2)
  # Values need not be integers.
  expect_identical(NNS.mode(rep(c(0.1, 0.25, 0.4), c(20, 50, 20))), 0.25)
})

test_that("a density highest at its minimum has its mode there", {
  e <- vapply(1:20, function(s) { set.seed(s); NNS.mode(rexp(1000), multi = FALSE) }, numeric(1))
  expect_lt(mean(e), 0.3)                                                # 0.64 before
  # Isolated extremes of a symmetric sample are not modes.
  k <- vapply(1:20, function(s) { set.seed(s); x <- rnorm(100); m <- NNS.mode(x); sum(m == min(x) | m == max(x)) }, numeric(1))
  expect_identical(sum(k), 0)
})

test_that("peaks are centred on the observations, not on bin edges", {
  # A sharp spike at 0.5 on a uniform background: the left-edge bin names put
  # the mode up to a bin (about 0.11) below it (0.418 to 0.493 before).
  m <- vapply(1:10, function(s) {
    set.seed(s); NNS.mode(c(0.5 + rnorm(1000, sd = 0.001), runif(1000, 0, 10)), multi = FALSE)
  }, numeric(1))
  expect_lt(max(abs(m - 0.5)), 0.01)
})

test_that("the regression mode reducer works on any scale", {
  # It used the discrete mode, which rounds to integers: on x in [0, 1] the
  # fit had 3 regression points and R2 0.72.
  set.seed(1); x <- runif(300); y <- sin(6 * x) + rnorm(300, sd = .1)
  f <- NNS.reg(x, y, noise.reduction = "mode", plot = FALSE, residual.plot = FALSE)
  expect_gt(nrow(f$regression.points), 10)
  expect_gt(f$R2, 0.9)
  g <- NNS.reg(100 * x, y, noise.reduction = "mode", plot = FALSE, residual.plot = FALSE)
  expect_equal(g$Fitted.xy$y.hat, f$Fitted.xy$y.hat, tolerance = 1e-8)
  p <- NNS.part(x, y, order = 2, noise.reduction = "mode")$regression.points
  expect_false(all(p$x == round(p$x)))
  # Classification keeps the discrete class mode.
  expect_identical(NNS:::mode_class(c(1, 1, 3, 3)), 1)
})
