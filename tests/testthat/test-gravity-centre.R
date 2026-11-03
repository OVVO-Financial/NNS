# NNS.gravity averages the median, a histogram mode, the mean and the midpoint
# of the quartiles. It used floor-position order statistics (the lower middle
# value for even n) and named histogram bins by their left edges, so it sat
# below the centre of symmetric data and did not reflect with the data (#117).

test_that("gravity is centred on symmetric data and reflects with it", {
  expect_identical(NNS.gravity(1:4), 2.5)                 # 2.25 before
  expect_identical(-NNS.gravity(-(1:4)), 2.5)             # 2.75 before
  g <- vapply(1:500, function(s) { set.seed(s); x <- rnorm(10); c(NNS.gravity(x), -NNS.gravity(-x)) }, numeric(2))
  expect_lt(abs(mean(g[1, ])), 0.05)                      # -0.176 before
  expect_lt(mean(abs(g[1, ] - g[2, ])), 0.1)              # 0.35 before
})

test_that("one outlier moves gravity a quarter as far as the mean", {
  expect_equal(NNS.gravity(c(rep(0, 99), 1e6)), 2500)
  expect_identical(NNS.gravity(c(1, 5, 9)), 5)            # up to three values: the median
})

test_that("default NNS.reg fits are not pulled below the data on pure noise", {
  # The default reducer is gravity: on noise its fits sat 0.047 below the
  # noise.reduction = "mean" fits.
  d <- vapply(1:40, function(s) {
    set.seed(s); x <- rnorm(200); y <- rnorm(200)
    mean(NNS.reg(x, y, plot = FALSE, residual.plot = FALSE)$Fitted.xy$y.hat -
           NNS.reg(x, y, noise.reduction = "mean", plot = FALSE, residual.plot = FALSE)$Fitted.xy$y.hat)
  }, numeric(1))
  expect_lt(abs(mean(d)), 0.02)
})
