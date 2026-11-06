# NNS.stack reports the dimension-reduction cutoff that keeps its tuned count,
# and NNS.reg accepts what NNS.stack reports, NA included.

kept <- function(equation) {
  vapply(equation, function(eq) {
    sum(abs(eq$Coefficient[eq$Variable != "DENOMINATOR"]) > 0)
  }, integer(1))
}

test_that("a two-class stack reports one cutoff that keeps its tuned count", {
  # As in the logistic regression example: a 0/1 response, stacked as 1/2,
  # whose cutoff is then passed to NNS.reg.
  set.seed(21)
  n <- 400
  x <- data.frame(matrix(rnorm(n * 6), n, dimnames = list(NULL, paste0("v", 1:6))))
  y <- rbinom(n, 1, plogis(1.5 * x$v1 - x$v2 + .5 * x$v3))
  s <- NNS.stack(IVs.train = x, DV.train = y + 1, type = "CLASS", status = FALSE,
                 ncores = 1, folds = 3)
  expect_length(s$NNS.dim.red.threshold, 1)
  expect_true(is.finite(s$NNS.dim.red.threshold))
  expect_lt(s$NNS.dim.red.count, 6)
  r <- NNS.reg(x, y, dim.red.method = "cor", threshold = s$NNS.dim.red.threshold,
               plot = FALSE)
  expect_equal(unname(kept(r$equation)), rep(s$NNS.dim.red.count, 2))
})

test_that("with more classes, the stack reports one cutoff per class", {
  set.seed(9)
  n <- 600
  x <- data.frame(matrix(rnorm(n * 6), n))
  cls <- factor(max.col(cbind(x$X1, x$X2, -x$X1 - x$X2 + rnorm(n, sd = .5))))
  s <- NNS.stack(IVs.train = x, DV.train = cls, type = "CLASS", status = FALSE,
                 ncores = 1, folds = 5)
  expect_named(s$NNS.dim.red.threshold, levels(cls))
  expect_lt(s$NNS.dim.red.count, 6)
  r <- NNS.reg(x, cls, dim.red.method = "cor", threshold = s$NNS.dim.red.threshold,
               plot = FALSE)
  expect_equal(unname(kept(r$equation)), rep(s$NNS.dim.red.count, 3))
})

test_that("NNS.reg applies a per-class threshold by class name or in class order", {
  # Each cutoff keeps two different predictors for its class; 0.8 applied to
  # versicolor would keep none.
  th <- c(setosa = .8, versicolor = .1, virginica = .7)
  fit <- function(threshold) {
    NNS.reg(iris[, 1:4], iris[, 5], dim.red.method = "cor", threshold = threshold,
            plot = FALSE)$equation
  }
  in.order <- fit(unname(th))
  expect_equal(fit(th[c(3, 1, 2)]), in.order)
  for (k in names(th)) expect_equal(in.order[[k]], fit(th[[k]])[[k]])
  expect_equal(unname(kept(in.order)), c(2L, 2L, 2L))
  expect_error(fit(c(.1, .2)), "one cutoff per class")
  expect_error(
    NNS.reg(iris[, 2:4], iris[, 1], dim.red.method = "cor", threshold = c(.1, .2),
            plot = FALSE),
    "single number for regression"
  )
})

test_that("NA thresholds and neighbour counts select the defaults", {
  x <- iris[, 2:4]
  y <- iris[, 1]
  expect_identical(
    NNS.reg(x, y, dim.red.method = "cor", threshold = NA, plot = FALSE)$equation,
    NNS.reg(x, y, dim.red.method = "cor", threshold = 0, plot = FALSE)$equation
  )
  expect_identical(
    NNS.reg(iris[, 1:4], iris[, 5], dim.red.method = "cor", threshold = NULL,
            plot = FALSE)$equation,
    NNS.reg(iris[, 1:4], iris[, 5], dim.red.method = "cor", threshold = 0,
            plot = FALSE)$equation
  )
  expect_identical(
    NNS.reg(x, y, n.best = NA, point.est = x[1:5, ], plot = FALSE, ncores = 1)$Point.est,
    NNS.reg(x, y, point.est = x[1:5, ], plot = FALSE, ncores = 1)$Point.est
  )
  # Without a dimension-reduction component the stack reports NA, which
  # NNS.reg accepts.
  s <- NNS.stack(IVs.train = x, DV.train = y, method = 1, status = FALSE,
                 ncores = 1, folds = 2, seed = 1)
  expect_true(is.na(s$NNS.dim.red.threshold))
  expect_error(
    NNS.reg(x, y, dim.red.method = "cor", threshold = s$NNS.dim.red.threshold,
            plot = FALSE),
    NA
  )
})

test_that("a classification stack reports the cell count when all cells win", {
  # A single predictor's only neighbour candidate is all cells.
  s <- suppressWarnings(NNS.stack(IVs.train = iris[, 3], DV.train = iris[, 5],
                                  type = "CLASS", status = FALSE, ncores = 1,
                                  folds = 3))
  expect_true(is.finite(s$NNS.reg.n.best))
  expect_gte(s$NNS.reg.n.best, 1)
})

test_that("the logistic regression example runs with a renamed test column", {
  # read.csv of an Excel "CSV UTF-8" test file names its first column X...Male.
  set.seed(23)
  n <- 300
  x <- data.frame(Male = rbinom(n, 1, .4), Pregnancy.Test = rbinom(n, 1, .2),
                  Folic.Acid = rbinom(n, 1, .2), Wine = rbinom(n, 1, .3),
                  Ginger.Ale = rbinom(n, 1, .2))
  y <- rbinom(n, 1, plogis(-1 + 2 * x$Pregnancy.Test + 1.5 * x$Folic.Acid - x$Wine))
  test <- setNames(x[1:40, ], c("X...Male", names(x)[-1]))
  s <- NNS.stack(IVs.train = x, DV.train = y + 1, IVs.test = test, balance = TRUE,
                 order = "max", type = "CLASS", status = FALSE, ncores = 1, folds = 3)
  expect_length(s$stack, 40)
  expect_true(is.finite(s$NNS.dim.red.threshold))
  r <- NNS.reg(x = x, y = y, point.est = test, plot = FALSE, dim.red.method = "cor",
               threshold = s$NNS.dim.red.threshold)
  expect_length(r$Point.est, 40)
  expect_true(all(kept(r$equation) >= s$NNS.dim.red.count))
})
