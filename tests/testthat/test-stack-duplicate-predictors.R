test_that("NNS.stack accepts duplicate predictor columns like NNS.reg (cbind(x, x))", {
  set.seed(123)
  x <- rnorm(100); y <- rnorm(100)

  # The cbind(x, x) dimension trick works for NNS.reg; NNS.stack must not
  # error with "[IVs.train] predictor names must be unique." any more.
  expect_error(NNS.reg(x = cbind(x, x), y = y, plot = FALSE), NA)
  res <- NNS.stack(IVs.train = cbind(x, x), DV.train = y,
                   method = 1, ncores = 1, folds = 2, status = FALSE)
  expect_type(res, "list")
  expect_true("reg" %in% names(res))
})

test_that("NNS.stack aligns duplicate-named test columns (cbind(x, x)) instead of erroring", {
  set.seed(123)
  x <- rnorm(60); y <- rnorm(60)
  train <- cbind(x, x)

  expect_error(
    NNS.stack(IVs.train = train, DV.train = y, IVs.test = train[1:5, ],
              method = 1, ncores = 1, folds = 2, status = FALSE),
    NA
  )
})
