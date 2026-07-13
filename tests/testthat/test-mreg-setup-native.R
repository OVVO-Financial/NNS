test_that("native mreg setup returns RPM and duplicate map smoke outputs", {
  X <- cbind(a = c(1,1,2,2), b = c(1,1,2,2)); y <- c(1,2,3,4)
  prep <- .nns_mreg_prepare_model(X, y, order = "max", use.native = TRUE)
  expect_true(is.data.frame(prep$RPM)); expect_identical(as.integer(NNS_duplicate_column_map_cpp(X)), c(1L,1L))
})
