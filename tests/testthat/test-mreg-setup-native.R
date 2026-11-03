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

test_that("binary regression coordinates preserve both terminal atoms", {
  X <- data.frame(a = rep(c(0, 1), 20), constant = 0)
  y <- -5 + 10 * X$a
  for (depth in list(NULL, 1L, 3L, "max")) {
    native <- .nns_mreg_prepare_model(X, y, order = depth, use.native = TRUE)
    reference <- .nns_mreg_prepare_model(X, y, order = depth, use.native = FALSE)
    expect_equal(native$RPM, reference$RPM)
    expect_equal(native$ids, reference$ids)
    expect_equal(native$row.ids, reference$row.ids)
    expect_equal(native$RPM$a, c(0, 1))
    expect_equal(native$RPM$y.hat, c(-5, 5))
  }
  fit <- function(native) withr::with_options(list(NNS.native.mreg = native),
    NNS.reg(X, y, type = "REG", point.est = X[1:2, ], n.best = 1,
             plot = FALSE, ncores = 1))
  expect_equal(fit(TRUE)$Point.est, c(-5, 5))
  expect_equal(fit(FALSE)$Point.est, c(-5, 5))
})

test_that("mixed factor and numeric regression does not collapse dummy columns", {
  X <- data.frame(group = factor(rep(c("left", "right"), 20)),
                   value = rep(1:20, each = 2))
  y <- ifelse(X$group == "left", -5, 5)
  fit <- function(native) withr::with_options(list(NNS.native.mreg = native),
    NNS.reg(X, y, type = "REG", point.est = X[c(1, 2, 39, 40), ],
             n.best = 1, plot = FALSE, ncores = 1))
  a <- fit(TRUE); b <- fit(FALSE)
  expect_equal(a$RPM, b$RPM)
  expect_true(all(a$RPM$group_left %in% c(0, 1)))
  expect_true(all(a$RPM$group_right %in% c(0, 1)))
  expect_equal(a$RPM$group_left + a$RPM$group_right, rep(1, nrow(a$RPM)))
  expect_equal(a$Point.est, c(-5, 5, -5, 5))
  expect_equal(a$Point.est, b$Point.est)
})

test_that("native exact support guard leaves ordinary two-boundary intervals intact", {
  # An interior observed value makes these two boundaries an interval, not
  # a binary support. Different cutpoints and constant columns also retain
  # their existing rightmost-closed interval IDs.
  for (case in list(
      list(x = c(0, .5, 1), cut = c(0, 1)),
      list(x = c(0, 1, 0, 1), cut = c(-.5, 1.5)),
      list(x = rep(1, 3), cut = c(0, 1)))) {
    X <- matrix(case$x, ncol = 1)
    native <- NNS_mreg_setup_cpp(X, seq_len(nrow(X)), list(case$cut), 0L, FALSE)
    expected <- as.character(findInterval(case$x, case$cut, rightmost.closed = TRUE))
    expect_equal(as.character(native$ids), expected)
  }
  X <- matrix(rep(c(0, 1), 2), ncol = 1)
  native <- NNS_mreg_setup_cpp(X, 1:4, list(c(0, 1)), 0L, FALSE)
  expect_equal(as.character(native$ids), rep(c("1", "2"), 2))
  expect_equal(unname(native$RPM), cbind(c(0, 1), c(2, 3)))
})
