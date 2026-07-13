test_that("native mreg setup returns RPM and duplicate map smoke outputs", {
  X <- cbind(a = c(1,1,2,2), b = c(1,1,2,2)); y <- c(1,2,3,4)
  prep <- .nns_mreg_prepare_model(X, y, order = "max", use.native = TRUE)
  expect_true(is.data.frame(prep$RPM)); expect_identical(as.integer(NNS_duplicate_column_map_cpp(X)), c(1L,1L))
})

test_that("native mreg setup preserves one-based rightmost-closed findInterval IDs", {
  X <- matrix(as.numeric(1:12), ncol = 1, dimnames = list(NULL, "x"))
  y <- as.numeric(1:12)
  boundaries <- list(as.numeric(1:12))
  native <- NNS_mreg_setup_cpp(X, y, boundaries, 3L, FALSE)
  reference_ids <- as.character(findInterval(X[, 1], boundaries[[1]],
                                             left.open = FALSE,
                                             rightmost.closed = TRUE))
  expect_identical(as.character(native$ids), reference_ids)
  expect_identical(tail(as.character(native$ids), 2), c("11", "11"))

  native_rpm <- as.data.frame(native$RPM)
  names(native_rpm) <- c("x", "y.hat")
  reference_rpm <- .nns_mreg_build_rpm(as.data.frame(X), y, reference_ids,
                                       "off", FALSE)
  expect_identical(native_rpm$x, reference_rpm$x)
  expect_identical(native_rpm$y.hat, reference_rpm$y.hat)

  reference_groups <- split(seq_len(nrow(X)), reference_ids)
  expect_identical(
    names(reference_groups),
    c("1", "10", "11", "2", "3", "4", "5", "6", "7", "8", "9")
  )
  merged_group_position <- match("11", names(reference_groups))
  expect_false(is.na(merged_group_position))
  expect_equal(native_rpm$x[merged_group_position], 11.5)
  expect_equal(native_rpm$y.hat[merged_group_position], 11.5)
  expect_equal(tail(native_rpm$x, 1), 9)
  expect_equal(tail(native_rpm$y.hat, 1), 9)
})

test_that("native interval IDs match R for below first, exact boundaries, and above final", {
  boundaries <- list(as.numeric(c(1, 2, 3, 4)))
  X <- matrix(as.numeric(c(0, 1, 1.5, 2, 3.5, 4, 5)), ncol = 1)
  native <- NNS_mreg_setup_cpp(X, seq_len(nrow(X)), boundaries, 3L, FALSE)
  reference_ids <- as.character(findInterval(X[, 1], boundaries[[1]],
                                             left.open = FALSE,
                                             rightmost.closed = TRUE))
  expect_identical(as.character(native$ids), reference_ids)
  expect_identical(as.character(native$ids), as.character(c(0, 1, 1, 2, 3, 3, 4)))
})

test_that("native interval IDs match R across multiple predictors and boundary counts", {
  X <- cbind(
    a = as.numeric(c(0, 1, 1.5, 2, 3.5, 4, 5, 4, 4)),
    b = as.numeric(c(-1, 0, 5, 10, 20, 20, 25, 20, 20)),
    c = rep(7, 9)
  )
  boundaries <- list(
    as.numeric(c(1, 2, 3, 4)),
    as.numeric(c(0, 10, 20)),
    as.numeric(7)
  )
  native <- NNS_mreg_setup_cpp(X, seq_len(nrow(X)), boundaries, 3L, FALSE)
  id_parts <- lapply(seq_along(boundaries), function(j) {
    findInterval(X[, j], boundaries[[j]], left.open = FALSE,
                 rightmost.closed = TRUE)
  })
  reference_ids <- do.call(paste, c(as.data.frame(id_parts), sep = "."))
  expect_identical(as.character(native$ids), reference_ids)

  native_rpm <- as.data.frame(native$RPM)
  names(native_rpm) <- c(colnames(X), "y.hat")
  reference_rpm <- .nns_mreg_build_rpm(as.data.frame(X), seq_len(nrow(X)),
                                       reference_ids, "off", FALSE)
  expect_identical(native_rpm$a, reference_rpm$a)
  expect_identical(native_rpm$b, reference_rpm$b)
  expect_identical(native_rpm$c, reference_rpm$c)
  expect_identical(native_rpm$y.hat, reference_rpm$y.hat)
})
