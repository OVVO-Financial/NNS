test_that("univariate fast helper matches point estimates from public path", {
  set.seed(13); x <- rnorm(20); y <- sin(x); test <- x[1:5]
  fast <- .nns_reg_univariate_fast(x, y, test)
  pub <- NNS.reg(x, y, point.est = test, plot = FALSE, residual.plot = FALSE, point.only = TRUE)
  expect_equal(fast$prediction, as.numeric(pub$Point.est), tolerance = 1e-12)
})
