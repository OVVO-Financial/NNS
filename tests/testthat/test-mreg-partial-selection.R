test_that("partial-selection and single-k kernels match references", {
  rpm_x <- rbind(
    c(0, 0), c(0, 2), c(2, 0), c(2, 2),
    c(1, 1), c(1, 1), c(3, 1), c(-1, 1)
  )
  y_reg <- c(0, 2, 4, 6, 3, 5, 7, -1)
  y_cls <- c(1, 1, 2, 2, 3, 3, 4, 4)
  Xtest <- rbind(c(1, 0), c(1, 1), c(1, 2), c(0.5, 1.5))
  mins <- apply(rpm_x, 2, min)
  maxs <- apply(rpm_x, 2, max)

  rpm_reg <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_reg)
  rpm_cls <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_cls)
  dist_names <- c("NNS", "L2", "L1")

  for (dist_code in 0:2) {
    dist_name <- dist_names[dist_code + 1L]

    for (k in c(1L, 2L, 4L, nrow(rpm_x), nrow(rpm_x) + 3L)) {
      kk <- min(k, nrow(rpm_x))

      ref_reg <- .nns_mreg_predict_reference(
        Xtest, rpm_reg, kk, dist_name, mins, maxs, FALSE
      )
      legacy_reg <- NNS_mreg_predict_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE
      )
      single_reg_1 <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
      )
      single_reg_2 <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
      )
      path_reg_1 <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
      )
      path_reg_2 <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
      )

      expect_equal(legacy_reg, ref_reg, tolerance = 1e-12)
      expect_equal(single_reg_1, ref_reg, tolerance = 1e-12)
      expect_equal(single_reg_2, single_reg_1, tolerance = 1e-12)
      expect_equal(path_reg_1[, kk], ref_reg, tolerance = 1e-12)
      expect_equal(path_reg_2, path_reg_1, tolerance = 1e-12)

      ref_cls <- .nns_mreg_predict_reference(
        Xtest, rpm_cls, kk, dist_name, mins, maxs, TRUE
      )
      legacy_cls <- NNS_mreg_predict_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE
      )
      single_cls_1 <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
      )
      single_cls_2 <- NNS_mreg_predict_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
      )
      path_cls_1 <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
      )
      path_cls_2 <- NNS_mreg_predict_path_v2_cpp(
        rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
      )

      expect_equal(legacy_cls, ref_cls, tolerance = 1e-12)
      expect_equal(single_cls_1, ref_cls, tolerance = 1e-12)
      expect_equal(single_cls_2, single_cls_1, tolerance = 1e-12)
      expect_equal(path_cls_1[, kk], ref_cls, tolerance = 1e-12)
      expect_equal(path_cls_2, path_cls_1, tolerance = 1e-12)
    }
  }
})

test_that("NNS.reg keeps both nearest rows in a batch", {
  x <- data.frame(a = c(-2, 2, seq(1e6, 1e8, length.out = 64)),
                  b = c(rep(0, 65), 1))
  y <- c(0, 10, rep(3, 64))
  q <- data.frame(a = 0, b = 0)
  batch <- q[rep(1, 4), , drop = FALSE]
  fit <- function(query, dist = "NNS") {
    NNS.reg(x, y, type = "REG", order = "max", n.best = 1, ncores = 1,
            point.est = query, dist = dist, plot = FALSE, point.only = TRUE)
  }

  one <- fit(q)
  ref <- .nns_mreg_predict_reference(as.matrix(batch), one$RPM, 1L, "NNS",
                                     vapply(x, min, 0), vapply(x, max, 0), FALSE)
  expect_equal(ref, rep(5, 4), tolerance = 1e-12)
  expect_equal(rep(one$Point.est, 4), ref, tolerance = 1e-12)
  expect_equal(fit(batch)$Point.est, ref, tolerance = 1e-12)
  expect_equal(fit(batch, "L1")$Point.est, ref, tolerance = 1e-12)
  expect_equal(fit(batch, "L2")$Point.est, ref, tolerance = 1e-12)
})
