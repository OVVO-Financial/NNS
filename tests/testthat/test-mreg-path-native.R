test_that("v2 all-k path matches legacy path on smoke data", {
  rpm <- cbind(c(0,1,2), c(0,1,4)); x <- cbind(c(.5,1.5)); mins <- 0; maxs <- 2
  old <- NNS_mreg_predict_path_cpp(rpm[,1,drop=FALSE], rpm[,2], x, 3L, 0L, mins, maxs, FALSE)
  nat <- NNS_mreg_predict_path_v2_cpp(rpm[,1,drop=FALSE], rpm[,2], x, 3L, 0L, mins, maxs, FALSE, 1L)
  expect_equal(nat, old, tolerance = 1e-12)
})

test_that("v2 all-k path is deterministic with multiple native threads", {
  rpm_x <- matrix(seq(0, 1, length.out = 40), ncol = 2)
  y <- seq_len(nrow(rpm_x)) / 10
  xt <- rpm_x[c(1, 5, 10), , drop = FALSE]
  one <- NNS_mreg_predict_path_v2_cpp(rpm_x, y, xt, 6L, 0L,
                                      apply(rpm_x, 2, min), apply(rpm_x, 2, max),
                                      FALSE, 1L)
  many <- NNS_mreg_predict_path_v2_cpp(rpm_x, y, xt, 6L, 0L,
                                       apply(rpm_x, 2, min), apply(rpm_x, 2, max),
                                       FALSE, 2L)
  expect_equal(many, one, tolerance = 1e-12)
})
