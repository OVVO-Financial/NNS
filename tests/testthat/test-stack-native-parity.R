test_that("native backend options are accepted", {
  old <- options(NNS.native.stack = TRUE, NNS.native.mreg = TRUE, NNS.native.univariate = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(isTRUE(getOption("NNS.native.stack")))
})
