test_that("regression uses observed class labels for codes and probabilities", {
  x <- rep(1:2, each = 8)
  y <- rep(c("b", "a"), each = 8)
  fit <- function(y) NNS.reg(x, y, point.est = 1:2, plot = FALSE)
  expected <- fit(factor(y, levels = c("b", "a")))
  observed <- fit(factor(y, levels = c("b", "missing", "a")))
  expect_identical(observed$class.levels, expected$class.levels)
  expect_equal(observed$Point.probabilities, expected$Point.probabilities)
  expect_identical(observed$Point.est, expected$Point.est)
  ordered <- fit(ordered(y, levels = c("b", "missing", "a")))
  expect_identical(ordered$class.levels, c("b", "a"))
  expect_identical(ordered$class.levels[ordered$Point.est], y[c(1, 16)])
})

test_that("classification interpolates class frequencies, not arbitrary class codes", {
  x <- rep(0:2, each = 4)
  labels <- rep(c("A", "C", "B"), each = 4)
  run <- function(levels) NNS.reg(x, factor(labels, levels = levels),
    order = "max", point.est = c(.25, .5, .75, 1.25, 1.5), plot = FALSE)
  a <- run(c("A", "B", "C")); b <- run(c("C", "A", "B"))
  expect_equal(a$class.levels[a$Point.est], b$class.levels[b$Point.est])
  expect_equal(unclass(a$Point.probabilities), b$Point.probabilities[, c("A", "B", "C")])
  expect_equal(a$Point.probabilities[2, ], c(A = .5, B = 0, C = .5))
  expect_equal(a$class.levels[a$Point.est[2]], "A")
  expect_equal(rowSums(a$Point.probabilities), rep(1, 5))
})

test_that("class cells retain counts and numeric labels keep their identity", {
  x <- factor(rep(c("left", "right"), each = 4))
  y <- c(10, 10, 10, -7, -7, -7, -7, 10)
  fit <- NNS.reg(x, y, type = "CLASS", point.est = c("left", "right"), plot = FALSE)
  expect_equal(fit$Point.est, c(10, -7))
  expect_equal(unname(unclass(fit$Point.probabilities)), rbind(c(.25, .75), c(.75, .25)))
  relabel <- ifelse(y == 10, -100, 8)
  other <- NNS.reg(x, relabel, type = "CLASS", point.est = c("left", "right"), plot = FALSE)
  expect_equal(other$Point.est, c(-100, 8))
})

test_that("stack probabilities and decoded decisions are invariant to nominal label order", {
  x <- iris[, 1:3]
  run <- function(y) NNS.stack(x, y, IVs.test = x[c(1, 52, 103), ],
    folds = 2, seed = 17, ncores = 1, status = FALSE)
  a <- run(iris$Species)
  b <- run(factor(iris$Species, levels = rev(levels(iris$Species))))
  expect_equal(a$class.levels[a$stack], b$class.levels[b$stack])
  expect_equal(unclass(a$probabilities), b$probabilities[, colnames(a$probabilities)])
  expect_equal(a$weights, b$weights)
  expect_equal(rowSums(a$probabilities), rep(1, 3))
  expect_true(all(is.finite(a$NNS.dim.red.threshold)))
  expect_true(all(is.finite(b$NNS.dim.red.threshold)))
  expect_true(a$NNS.dim.red.count %in% 1:3)
})

test_that("numeric stack outputs and custom objectives see original class values", {
  x <- iris[1:100, 1:2]
  y <- rep(c(-5, 8), each = 50)
  observed <- numeric()
  record <- function(actual) observed <<- c(observed, actual)
  loss <- as.expression(bquote({ .(record)(actual); mean(predicted == actual) }))
  fit <- NNS.stack(x, y, x[c(1, 100), ], method = 1, obj.fn = loss,
    objective = "max", folds = 2, ncores = 1, status = FALSE, seed = 9)
  expect_setequal(unique(observed), c(-5, 8))
  expect_equal(fit$stack, c(-5, 8))
})

test_that("nominal predictor level order does not alter categorical group regression", {
  y <- rep(c(2.1, 5.3, -1.8), each = 6) + rep(c(-.2, 0, .2), 6)
  labels <- rep(c("a", "b", "c"), each = 6)
  run <- function(levels) NNS.reg(factor(labels, levels = levels), y,
    point.est = c("a", "b", "c"), noise.reduction = "mean", plot = FALSE)
  a <- run(c("a", "b", "c")); b <- run(c("c", "a", "b"))
  expect_equal(a$Point.est, b$Point.est)
  expect_equal(a$Point.est, c(2.1, 5.3, -1.8))
})

test_that("both stack components canonicalize unordered predictor levels", {
  set.seed(128)
  x <- data.frame(a = rnorm(90), group = factor(rep(c("z", "a", "m"), 30)))
  reordered <- x; reordered$group <- factor(x$group, levels = c("m", "z", "a"))
  response <- x$a + rep(c(1, 2, -1), 30) + rnorm(90, sd = .1)
  run <- function(data, y) NNS.stack(data, y, IVs.test = data[1:6, ],
    method = c(1, 2), folds = 2, seed = 19, ncores = 1, status = FALSE)
  expect_equal(run(x, response), run(reordered, response))
  target <- factor(rep(c("B", "A", "C"), 30))
  expect_equal(run(x, target), run(reordered, target))
})

test_that("two-valued numeric coordinates retain distinct class cells", {
  x <- rep(c(0, 1), 20)
  y <- factor(ifelse(x == 0, "zero", "one"))
  fit <- NNS.reg(x, y, point.est = c(0, .5, 1), plot = FALSE)
  expect_equal(nrow(fit$RPM), 2L)
  expect_equal(fit$Prediction.Accuracy, 1)
  expect_equal(unname(unclass(fit$Point.probabilities)),
               rbind(c(0, 1), c(.5, .5), c(1, 0)))
  mixed <- data.frame(group = factor(x), continuous = rep(1:20, each = 2))
  cells <- NNS.reg(mixed, y, point.est = mixed[1:2, ], n.best = 1,
                   plot = FALSE, ncores = 1)
  expect_equal(cells$class.levels[cells$Point.est], c("zero", "one"))
  expect_true(all(apply(cells$RPM.probabilities, 1, max) == 1))
})

test_that("named dimension weights align in regression and nominal classification", {
  set.seed(12)
  x <- data.frame(a = rep(c(-1, 1), 20), b = rnorm(40))
  labels <- factor(ifelse(x$a > 0, "A", "B"))
  for (response in list(labels, seq_len(nrow(x)) / 10)) {
    fit <- function(w) NNS.reg(x, response, point.est = x[1:4, ],
      dim.red.method = w, plot = FALSE)
    a <- fit(c(a = 1, b = 0)); b <- fit(c(b = 0, a = 1))
    expect_equal(a$equation, b$equation)
    expect_equal(a$Point.est, b$Point.est)
    expect_equal(a$Point.probabilities, b$Point.probabilities)
    expect_error(fit(c(a = 1, unknown = 0)), "exactly match")
  }
  fit <- function(w) NNS.stack(x, labels, x[1:4, ], dim.red.method = w,
    method = c(1, 2), folds = 2, ncores = 1, seed = 1, status = FALSE)
  expect_equal(fit(c(a = 1, b = 0)), fit(c(b = 0, a = 1)))
  expect_error(fit(c(a = 1, unknown = 0)), "exactly match")
})

test_that("nominal Stack shares the factor coefficient names and fold columns", {
  x <- data.frame(group = factor(rep(c("a", "b"), 20)),
                  value = rep(1:20, each = 2))
  y <- factor(ifelse(x$group == "a", "L", "R"))
  fit <- function(w) NNS.stack(x, y, x[1:4, ], dim.red.method = w,
    method = 2, folds = 2, ncores = 1, seed = 9, status = FALSE)
  weights <- c(group__a = 1, group__b = -1, value = 0)
  a <- fit(weights); b <- fit(rev(weights))
  expect_equal(a, b)
  expect_equal(a$class.levels[a$stack], as.character(y[1:4]))
  expect_equal(a$NNS.dim.red.coefficients[[1]]$Variable,
               c(names(weights), "DENOMINATOR"))
  # A level present only in one validation row keeps its design column even
  # when absent from that fold's fitted observations.
  x$group <- as.character(x$group); x$group[1] <- "rare"
  x$group <- factor(x$group)
  rare <- NNS.stack(x, y, x[1:4, ], method = c(1, 2), folds = 2,
                    ncores = 1, seed = 9, status = FALSE)
  expect_equal(rowSums(rare$probabilities), rep(1, 4))
})

test_that("Stack breaks binary ties in training order for each holdout", {
  y <- rep(c(1, 0), 12)
  x <- data.frame(v = rep(0, length(y)), w = 1)
  fit <- NNS.stack(x, y, x[1, , drop = FALSE], CV.size = .25, folds = 1,
    seed = 16, optimize.threshold = FALSE, ncores = 1, status = FALSE,
    obj.fn = expression(mean(predicted[actual == 1] == 1)), objective = "max")
  # Row 1 is held out. Class 0 appears first in the fold and class 1 in the refit.
  expect_equal(c(fit$OBJfn.reg, fit$OBJfn.dim.red), c(0, 0))
  expect_equal(c(fit$reg, fit$dim.red, fit$stack), c(1, 1, 1))

  balanced <- NNS.stack(x, y, x[1, , drop = FALSE], CV.size = .25, folds = 1,
    balance = TRUE, seed = 4, pred.int = .5, optimize.threshold = FALSE, ncores = 1, status = FALSE,
    obj.fn = expression(mean(predicted[actual == 1] == 1)), objective = "max")
  # The balanced training samples for both the holdout and the refit start with class 0.
  expect_equal(c(balanced$OBJfn.reg, balanced$OBJfn.dim.red), c(0, 0))
  expect_equal(c(balanced$reg, balanced$dim.red, balanced$stack), c(0, 0, 0))
  expect_equal(balanced$Prediction.set[[1]], 0)
})

test_that("pooled class ties follow the earliest training row in either holdout order", {
  y <- rep(c(1, 0), 12)
  x <- data.frame(v = rep(0, length(y)))
  # Expected decisions for validation rows 1, 5, 8, 9, 12, 15, 16, 17, 19, 21, 22, 24.
  fit <- function() NNS.stack(x, y, method = 1, CV.size = .25, folds = 3,
    seed = 16, optimize.threshold = FALSE, ncores = 1, status = FALSE,
    obj.fn = expression(mean(predicted == c(0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1))))
  expect_equal(fit()$OBJfn.reg, 1)

  classify <- NNS:::.nns_stack_classify
  local_mocked_bindings(.nns_stack_classify = function(x, y, z, class.values, splits, ...) {
    classify(x, y, z, class.values, rev(splits), ...)
  }, .package = "NNS")
  expect_equal(fit()$OBJfn.reg, 1)
})

test_that("overlapping holdouts preserve an exact average at the class cutoff", {
  y <- rep(c(1, 0), 25)
  set.seed(3893)
  x <- data.frame(v = sample(0:1, length(y), replace = TRUE))
  classify <- NNS:::.nns_stack_classify
  score <- function(permutation) {
    local_mocked_bindings(.nns_stack_classify = function(x, y, z, class.values, splits, ...) {
      classify(x, y, z, class.values, splits[permutation], ...)
    }, .package = "NNS")
    fit <- function(x, method) NNS.stack(x, y, method = method, seed = 993,
      CV.size = .44, folds = 3, optimize.threshold = FALSE, ncores = 1, status = FALSE,
      obj.fn = expression({
        # Row 20 is scored 15th. Its class-0 probabilities are 7/12, 1/2 and 5/12.
        expect_equal(predicted[15], 1)
        mean(predicted == actual)
      }))
    c(fit(x, 1)$OBJfn.reg, fit(transform(x, w = 1), 2)$OBJfn.dim.red)
  }
  expect_equal(score(c(1, 2, 3)), c(16/41, 16/41))
  expect_equal(score(c(1, 3, 2)), c(16/41, 16/41))
  expect_equal(score(c(2, 1, 3)), c(16/41, 16/41))
  expect_equal(score(c(2, 3, 1)), c(16/41, 16/41))
  expect_equal(score(c(3, 1, 2)), c(16/41, 16/41))
  expect_equal(score(c(3, 2, 1)), c(16/41, 16/41))
})

test_that("Stack breaks multiclass ties in training order", {
  y <- rep(c(2, 1, 3), 12)
  x <- data.frame(v = rep(0, length(y)))
  fit <- NNS.stack(x, y, x[1, , drop = FALSE], type = "CLASS", method = 1,
    CV.size = .25, folds = 1, seed = 33002, ncores = 1, status = FALSE,
    obj.fn = expression(mean(predicted[actual == 2] == 2)), objective = "max")
  expect_equal(fit$OBJfn.reg, 0)
  expect_equal(fit$stack, 2)
})

test_that("Stack reports the multivariate coordinates separately from the class indices", {
  set.seed(13)
  x <- data.frame(signal = rep(c(0, 1), 40), nuisance = rnorm(80))
  fit <- NNS.stack(x, x$signal, data.frame(signal = c(0, 1), nuisance = 0),
    dim.red.method = c(signal = 0, nuisance = 1), ncores = 1, status = FALSE)
  expect_equal(unname(fit$NNS.dim.red.variables), "nuisance")
  expect_equal(fit$NNS.reg.variables, "signal")
  expect_equal(fit$stack, c(0, 1))
})

test_that("class codes give the same accuracy as labels for an unsorted factor", {
  x <- iris[1:100, 1:2]
  y <- factor(iris$Species[1:100], levels = c("versicolor", "setosa"))
  rows <- c(1, 100)
  query <- x[rows, ]
  actual <- y[rows]

  reg <- NNS.reg(x, y, point.est = query, plot = FALSE)
  expect_identical(reg$Point.est == as.numeric(actual),
                   reg$class.levels[reg$Point.est] == as.character(actual))

  stack <- NNS.stack(x, y, query, method = 1, folds = 2, seed = 17,
                     ncores = 1, status = FALSE)
  expect_identical(stack$stack == as.numeric(actual),
                   stack$class.levels[stack$stack] == as.character(actual))

  boost <- NNS.boost(x, y, query, epochs = 0, folds = 2, learner.trials = 3,
                     seed = 9, feature.importance = FALSE, status = FALSE)
  expect_identical(boost$results == as.numeric(actual),
                   boost$class.levels[boost$results] == as.character(actual))
})

test_that("logical class predictions keep their truth values", {
  x <- rep(1:2, each = 8)
  y <- x == 2
  fit <- NNS.reg(x, y, point.est = 1:2, plot = FALSE, confidence.interval = .8)
  expect_identical(fit$Point.est, c(FALSE, TRUE))
  expect_identical(fit$Fitted.xy$y.hat, y)
  expect_identical(unlist(fit$Prediction.set), c(FALSE, TRUE))
})

test_that("a declared NA class keeps its own probability column", {
  x <- rep(1:3, each = 5)
  y <- factor(rep(c("b", NA, "a"), each = 5), exclude = NULL)
  run <- function(y) NNS.reg(x, y, order = "max", point.est = 1:3, plot = FALSE)
  fit <- run(y)
  expect_identical(fit$class.levels[fit$Point.est], as.character(y[c(1, 6, 11)]))
  expect_equal(unname(unclass(fit$Point.probabilities)),
               unname(unclass(run(factor(as.integer(y)))$Point.probabilities)))
  is.na(y)[1] <- TRUE
  expect_error(run(y), "missing values")
})

test_that("per-class thresholds match a declared NA class by name", {
  y <- factor(rep(c("b", NA, "a"), each = 10), levels = c("b", NA, "a"), exclude = NULL)
  set.seed(201)
  x <- data.frame(signal = as.integer(y) + rnorm(30), noise = rnorm(30))
  fit <- function(th) NNS.reg(x, y, point.est = x[c(1, 11, 21), ], threshold = th,
                              dim.red.method = "cor", plot = FALSE)
  th <- setNames(c(0, .2, .7), levels(y))
  expect_equal(fit(rev(th))$Point.probabilities, fit(unname(th))$Point.probabilities,
               tolerance = 1e-12)
})
