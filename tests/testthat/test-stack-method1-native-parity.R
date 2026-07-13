test_that("method 1 native and reference select identical n.best on regression smoke case", {
  set.seed(11); x <- data.frame(a = rnorm(24), b = rnorm(24)); y <- x$a - x$b^2
  ref <- stack_with_native(FALSE, NNS.stack(x, y, IVs.test = x[1:3,], method = 1, folds = 2, ncores = 1, status = FALSE))
  nat <- stack_with_native(TRUE, NNS.stack(x, y, IVs.test = x[1:3,], method = 1, folds = 2, ncores = 1, status = FALSE))
  expect_identical(nat$NNS.reg.n.best, ref$NNS.reg.n.best)
  expect_numeric_close(nat$reg, ref$reg)
})
