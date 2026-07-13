test_that("native xstar path matches active coefficient projections including zero fallback", {
  train <- matrix(1:12, 4); test <- matrix(13:18, 2); coef <- c(2, 0, -1)
  ord <- order(abs(coef), decreasing = TRUE, na.last = NA, method = "radix")
  out <- NNS_xstar_path_cpp(train, test, coef, ord, 1L)
  expect_equal(dim(out$train), c(4L, 3L)); expect_equal(out$column_order, ord)
  z <- NNS_xstar_path_cpp(train, test, c(0,0,0), seq_len(3), 1L)
  expect_equal(z$train[,2], rowMeans(train[,1:2]))
})

test_that("native xstar path is deterministic with multiple threads", {
  train <- matrix(seq_len(60), 20, 3); test <- matrix(seq_len(15), 5, 3)
  coef <- c(3, 0, -2); ord <- order(abs(coef), decreasing = TRUE, na.last = NA, method = "radix")
  one <- NNS_xstar_path_cpp(train, test, coef, ord, 1L)
  many <- NNS_xstar_path_cpp(train, test, coef, ord, 2L)
  expect_equal(many$train, one$train, tolerance = 1e-14)
  expect_equal(many$test, one$test, tolerance = 1e-14)
  expect_identical(many$representative, one$representative)
})
