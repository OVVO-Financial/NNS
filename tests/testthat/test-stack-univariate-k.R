capture_frame <- function(fun, call) {
  captured <- new.env()
  suppressMessages(trace(
    fun, where = asNamespace("NNS"), print = FALSE,
    exit = bquote(assign("frame", environment(), envir = .(captured)))
  ))
  on.exit(suppressMessages(untrace(fun, where = asNamespace("NNS"))))
  force(call)
  captured$frame
}

wt <- mtcars[, "wt", drop = FALSE]

test_that("NNS.stack reports no n.best for a single predictor", {
  s <- NNS.stack(wt, mtcars$mpg, method = 1, seed = 1, ncores = 1, status = FALSE)
  expect_identical(s$NNS.reg.n.best, NA_integer_)
})

test_that("NNS.stack OBJfn.reg for a single predictor scores the deployed NNS.reg fit", {
  y <- mtcars$mpg
  frame <- capture_frame("NNS.stack", s <- NNS::NNS.stack(
    wt, y, method = 1, seed = 1, ncores = 1, status = FALSE
  ))
  univariate_fit <- function(train_x, train_y, test_x) {
    NNS.reg(train_x, train_y, point.est = test_x, order = frame$order, type = "REG",
            factor.2.dummy = FALSE, plot = FALSE, residual.plot = FALSE,
            ncores = 1)$Point.est
  }
  oof <- rep(NA_real_, length(y))
  for (sp in frame$splits) {
    design <- frame$.method1_design_for_split(sp$train, sp$validation)
    oof[sp$validation] <- univariate_fit(design$train[[1L]], frame$y_fit_all[sp$train],
                                         design$test[[1L]])
  }
  expect_equal(s$OBJfn.reg, sum((oof - y)^2), tolerance = 1e-12)
  expect_equal(s$reg, univariate_fit(frame$full_design$train[[1L]], frame$y_fit_all,
                                     frame$full_design$test[[1L]]), tolerance = 1e-12)
})
