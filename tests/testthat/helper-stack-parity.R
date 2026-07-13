stack_with_native <- function(native, expr) {
  old <- options(NNS.native.stack = native, NNS.native.mreg = native,
                 NNS.native.univariate = native)
  on.exit(options(old), add = TRUE)
  force(expr)
}
expect_numeric_close <- function(x, y, tol = 1e-12) {
  expect_equal(as.numeric(x), as.numeric(y), tolerance = tol, scale = 1)
}
