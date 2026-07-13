test_that("native mreg setup returns RPM and duplicate map smoke outputs", {
  X <- cbind(a = c(1,1,2,2), b = c(1,1,2,2)); y <- c(1,2,3,4)
  prep <- .nns_mreg_prepare_model(X, y, order = "max", use.native = TRUE)
  expect_true(is.data.frame(prep$RPM)); expect_identical(as.integer(NNS_duplicate_column_map_cpp(X)), c(1L,1L))
})

test_that("native mreg setup preserves one-based findInterval IDs and RPM order", {
  X <- matrix(seq_len(12), ncol = 1, dimnames = list(NULL, "x"))
  y <- seq_len(12)
  boundaries <- list(seq_len(12))
  native <- NNS_mreg_setup_cpp(X, y, boundaries, 3L, FALSE)
  reference_ids <- as.character(findInterval(X[, 1], boundaries[[1]],
                                             left.open = FALSE,
                                             rightmost.closed = TRUE))
  expect_identical(as.character(native$ids), reference_ids)

  native_rpm <- as.data.frame(native$RPM)
  names(native_rpm) <- c("x", "y.hat")
  reference_rpm <- .nns_mreg_build_rpm(as.data.frame(X), y, reference_ids,
                                       "off", FALSE)
  expect_identical(native_rpm$x, reference_rpm$x)
  expect_identical(native_rpm$y.hat, reference_rpm$y.hat)
})
