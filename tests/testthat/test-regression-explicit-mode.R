test_that("REG keeps numeric binary responses on the regression surface", {
  x <- c(0, 1, 2, 3)
  y <- c(0, 0, 1, 1)
  fit <- NNS.reg(x, y, type = "REG", order = "max",
                 point.est = 1.5, plot = FALSE)
  expect_equal(fit$Point.est, .5)
  expect_null(fit$Prediction.Accuracy)
  expect_null(fit$class.levels)
  expect_equal(NNS.reg(x, y, type = "reg", order = "max",
                       point.est = 1.5, plot = FALSE)$Point.est, .5)
  expect_type(NNS.reg(x, y, plot = FALSE)$Prediction.Accuracy, "double")
  expect_error(NNS.reg(x, factor(y), type = "REG", plot = FALSE),
               "requires a numeric response")
})

test_that("REG reaches full multivariate and dimension reduction fits", {
  t <- seq_len(40)
  x <- data.frame(a = t, b = sin(t / 5))
  y <- rep(c(0, 1), 20)
  for (method in list(NULL, "cor", "NNS.dep")) {
    fit <- NNS.reg(x, y, type = "REG", dim.red.method = method,
                   point.est = x[1:3, ], plot = FALSE, ncores = 1)
    expect_null(fit$Prediction.Accuracy)
    expect_null(fit$class.levels)
    expect_true(all(is.finite(fit$Point.est)))
    if (!is.null(method)) {
      expect_s3_class(fit$equation, "data.frame")
      expect_length(fit$equation$Coefficient, ncol(x) + 1L)
      expect_true(all(is.finite(fit$equation$Coefficient)))
    }
  }
})

test_that("Stack REG preserves numeric response scale in both components", {
  t <- seq_len(40)
  x <- data.frame(a = t, b = sin(t / 5))
  run <- function(y) NNS.stack(x, y, IVs.test = x[1:3, ], type = "REG",
    method = c(1, 2), folds = 2, ncores = 1, seed = 17, status = FALSE,
    pred.int = .8)
  binary <- run(rep(c(0, 1), 20))
  signed <- run(rep(c(-5, 5), 20))
  expect_null(binary$class.levels)
  expect_null(signed$class.levels)
  for (component in c("reg", "dim.red", "stack"))
    expect_equal(signed[[component]], 10 * binary[[component]] - 5,
                 tolerance = 1e-10)
  for (component in c("reg.pred.int", "dim.red.pred.int", "pred.int"))
    expect_equal(as.matrix(signed[[component]]),
                 10 * as.matrix(binary[[component]]) - 5, tolerance = 1e-10)
})

test_that("numeric derivative responses use regression even when binary", {
  y <- rep(c(0, 1), 15)
  expect_equal(dy.dx(seq_along(y), y, eval.point = "overall"), 1 / 29)
  expect_equal(dy.dx(seq_along(y), 10 * y - 5, eval.point = "overall"), 10 / 29)
})
