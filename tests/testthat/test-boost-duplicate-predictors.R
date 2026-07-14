test_that("NNS.boost accepts duplicate predictor columns like NNS.reg (cbind(x, x))", {
  set.seed(123)
  x <- rnorm(60)
  y <- as.numeric(x > 0)
  train <- cbind(x, x)

  # Same cbind(x, x) dimension trick as NNS.reg; NNS.boost must not error with
  # "[IVs.train] predictor names must be unique." any more, and duplicate-named
  # test columns must align rather than error.
  expect_error(
    NNS.boost(IVs.train = train, DV.train = y, IVs.test = train[1:5, ],
              learner.trials = 5, epochs = 5, CV.size = 0.25,
              status = FALSE),
    NA
  )
})
