test_that("method 2 native and reference select identical dimension count on smoke case", {
  set.seed(12); x <- data.frame(a = rnorm(24), b = rnorm(24), c = rnorm(24)); y <- x$a + x$b
  ref <- stack_with_native(FALSE, NNS.stack(x, y, IVs.test = x[1:3,], method = 2, folds = 2, ncores = 1, status = FALSE))
  nat <- stack_with_native(TRUE, NNS.stack(x, y, IVs.test = x[1:3,], method = 2, folds = 2, ncores = 1, status = FALSE))
  expect_identical(nat$NNS.dim.red.threshold, ref$NNS.dim.red.threshold)
  expect_numeric_close(nat$dim.red, ref$dim.red)
})
