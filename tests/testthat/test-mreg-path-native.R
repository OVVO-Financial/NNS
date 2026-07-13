test_that("v2 all-k path matches legacy path on smoke data", {
  rpm <- cbind(c(0,1,2), c(0,1,4)); x <- cbind(c(.5,1.5)); mins <- 0; maxs <- 2
  old <- NNS_mreg_predict_path_cpp(rpm[,1,drop=FALSE], rpm[,2], x, 3L, 0L, mins, maxs, FALSE)
  nat <- NNS_mreg_predict_path_v2_cpp(rpm[,1,drop=FALSE], rpm[,2], x, 3L, 0L, mins, maxs, FALSE, 1L)
  expect_equal(nat, old, tolerance = 1e-12)
})
