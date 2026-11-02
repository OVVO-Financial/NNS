test_that("1-D fast path matches generic distances with unequal ranges", {
  x <- c(-2, -1, 0, 0, 1, 2)
  rpm_x <- cbind(x, x)
  y_reg <- c(4, 1, 2, 6, 1, 4)
  y_cls <- c(1, 1, 2, 2, 3, 3)
  xt_scalar <- c(-0.5, 0, 0.5, 3)
  Xtest <- cbind(xt_scalar, xt_scalar)

  # Deliberately inconsistent supplied ranges for duplicated coordinates.
  # The 1-D shortcut must still reproduce the generic per-column scaling.
  mins <- c(-2, -4)
  maxs <- c(2, 8)

  rpm_reg <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_reg)
  rpm_cls <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_cls)
  dist_names <- c("NNS", "L2", "L1")

  for (dist_code in 0:2) {
    dist_name <- dist_names[dist_code + 1L]

    for (k in c(1L, 4L, nrow(rpm_x), nrow(rpm_x) + 2L)) {
      kk <- min(k, nrow(rpm_x))

      reference <- .nns_mreg_predict_reference(
        Xtest, rpm_reg, kk, dist_name, mins, maxs, FALSE
      )
      legacy <- NNS_mreg_predict_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE
      )
      single_one <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
      )
      single_many <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
      )
      path_one <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
      )
      path_many <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
      )

      expect_equal(legacy, reference, tolerance = 1e-12)
      expect_equal(single_one, reference, tolerance = 1e-12)
      expect_equal(single_many, single_one, tolerance = 1e-12)
      expect_equal(path_one[, kk], reference, tolerance = 1e-12)
      expect_equal(path_many, path_one, tolerance = 1e-12)

      reference_cls <- .nns_mreg_predict_reference(
        Xtest, rpm_cls, kk, dist_name, mins, maxs, TRUE
      )
      legacy_cls <- NNS_mreg_predict_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE
      )
      single_cls_one <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
      )
      single_cls_many <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
      )
      path_cls_one <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
      )
      path_cls_many <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
      )

      expect_equal(legacy_cls, reference_cls, tolerance = 1e-12)
      expect_equal(single_cls_one, reference_cls, tolerance = 1e-12)
      expect_equal(single_cls_many, single_cls_one, tolerance = 1e-12)
      expect_equal(path_cls_one[, kk], reference_cls, tolerance = 1e-12)
      expect_equal(path_cls_many, path_cls_one, tolerance = 1e-12)
    }
  }
})

test_that("1-D fast path does not depend on the other queries in the batch", {
  # A non-duplicated query turns the 1-D shortcut off for the whole batch.
  # After range scaling (.07, .07) is equidistant from .01 and .13 under NNS.
  rpm_x <- cbind(c(.01, .13, 1), c(.01, .13, 1))
  y <- c(0, 100, 999)
  mins <- c(.01, .01)
  maxs <- c(1, 1)
  one <- matrix(.07, 1, 2)
  batch <- rbind(c(.07, .07), c(.2, .3))
  rpm <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y)
  expect_equal(.nns_mreg_predict_reference(one, rpm, 1L, "NNS", mins, maxs, FALSE), 50)

  set.seed(7)
  g <- round(runif(60) * 50) / 50 + .01
  gx <- cbind(g, g)
  gy <- round(sin(9 * g), 3)
  gc <- (seq_along(g) %% 3) + 1
  q <- round(runif(100) * 100) / 100
  gmaxs <- c(1, 1.3)
  q_one <- cbind(c(.07, q), c(.07, q))
  q_batch <- rbind(q_one, c(.2, .3))
  m <- nrow(q_one)

  for (dist_code in 0:2) {
    for (k in 1:3) {
      single <- NNS_mreg_predict_v2_cpp(rpm_x, y, one, k, dist_code, mins, maxs, FALSE, 1L)
      batched <- NNS_mreg_predict_v2_cpp(rpm_x, y, batch, k, dist_code, mins, maxs, FALSE, 1L)
      expect_identical(single, batched[1])
    }
    for (cls in c(FALSE, TRUE)) {
      yy <- if (cls) gc else gy
      p1 <- NNS_mreg_predict_path_v2_cpp(gx, yy, q_one, 8L, dist_code, mins, gmaxs, cls, 1L)
      pb <- NNS_mreg_predict_path_v2_cpp(gx, yy, q_batch, 8L, dist_code, mins, gmaxs, cls, 1L)
      expect_identical(p1, pb[seq_len(m), , drop = FALSE])
      for (k in c(1L, 4L, 8L)) {
        s1 <- NNS_mreg_predict_v2_cpp(gx, yy, q_one, k, dist_code, mins, gmaxs, cls, 1L)
        sb <- NNS_mreg_predict_v2_cpp(gx, yy, q_batch, k, dist_code, mins, gmaxs, cls, 1L)
        expect_identical(s1, sb[seq_len(m)])
        expect_identical(s1, p1[, k])
      }
    }
  }
  expect_identical(NNS_mreg_predict_v2_cpp(rpm_x, y, one, 1L, 0L, mins, maxs, FALSE, 1L), 50)
})
