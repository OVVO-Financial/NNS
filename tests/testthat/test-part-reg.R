test_that("NNS.part keeps the x centroid of integer-valued x", {
  p <- NNS.part(1:20, 2 * (1:20), type = "XONLY", order = 2,
                noise.reduction = "mean", obs.req = 0)
  expect_equal(p$regression.points$x, c(5.5, 15.5))
  expect_equal(p$regression.points$y, c(11, 31))
})

test_that("NNS.reg fits an integer-valued line exactly", {
  x <- 1:10
  fit <- NNS.reg(x, 3 * x + 1, plot = FALSE, residual.plot = FALSE, ncores = 1)
  expect_equal(fit$R2, 1)
  expect_equal(fit$Fitted.xy$y.hat, 3 * x + 1)
})

test_that("native and pure-R multivariate mode reducers agree", {
  set.seed(4)
  X <- matrix(round(rnorm(120), 1), ncol = 2)
  y <- round(X[, 1] + X[, 2]^2 + rnorm(60, sd = 0.3), 1)
  fit <- function(native) {
    op <- options(NNS.native.mreg = native)
    on.exit(options(op))
    NNS.reg(X, y, noise.reduction = "mode", point.est = X[1:4, ],
            plot = FALSE, residual.plot = FALSE, ncores = 1)
  }
  a <- fit(TRUE)
  b <- fit(FALSE)
  expect_equal(a$RPM, b$RPM)
  expect_equal(a$Fitted.xy$y.hat, b$Fitted.xy$y.hat)
  expect_equal(a$Point.est, b$Point.est)
})
