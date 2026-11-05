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
