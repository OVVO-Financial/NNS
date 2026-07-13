stack_with_native <- function(native, expr) {
  old <- options(NNS.native.stack = native, NNS.native.mreg = native,
                 NNS.native.univariate = native)
  on.exit(options(old), add = TRUE)
  force(expr)
}
expect_numeric_close <- function(actual, expected, tolerance = 1e-12) {
  testthat::expect_identical(dim(actual), dim(expected))
  testthat::expect_identical(names(actual), names(expected))
  a <- as.numeric(actual)
  b <- as.numeric(expected)
  testthat::expect_identical(is.na(a), is.na(b))
  keep <- is.finite(a) & is.finite(b)
  if (any(keep)) {
    testthat::expect_lte(max(abs(a[keep] - b[keep])), tolerance)
  }
}
