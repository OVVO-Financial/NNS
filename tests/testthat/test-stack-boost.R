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

boost_binary <- function(lo, hi) {
  y <- ifelse(iris$Species[1:140] == "setosa", hi, lo)
  NNS.boost(iris[1:140, 1:4], y, IVs.test = iris[c(1:3, 141:143), 1:4],
            epochs = 2, learner.trials = 3, feature.importance = FALSE,
            status = FALSE, seed = 1, pred.int = 0.9)
}

test_that("NNS.boost returns 0/1 response values for a 0/1 numeric response", {
  f <- boost_binary(0, 1)
  expect_equal(f$results, c(1, 1, 1, 0, 0, 0))
  expect_equal(unname(as.matrix(f$pred.int)), cbind(c(1, 1, 1, 0, 0, 0), c(1, 1, 1, 0, 0, 0)))
})

test_that("NNS.boost returns -5/5 response values for a -5/5 numeric response", {
  # Near-zero epoch scores can trigger the best-trial fallback by roundoff.
  # Accept only that documented warning, without requiring it on every platform.
  fallback_warning <- paste0(
    "No feature combination re-passed the objective cutoff during epochs; ",
    "using the best learner-trial combination. ",
    "Consider a lower [threshold] probability."
  )
  f <- withCallingHandlers(
    boost_binary(-5, 5),
    warning = function(w) {
      expect_identical(conditionMessage(w), fallback_warning)
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(f$results, c(5, 5, 5, -5, -5, -5))
  expect_equal(unname(as.matrix(f$pred.int)), cbind(c(5, 5, 5, -5, -5, -5), c(5, 5, 5, -5, -5, -5)))
})

test_that("NNS.stack codes a character response like the equivalent factor", {
  idx <- c(101:140, 1:100)
  x <- iris[idx, 1:4]
  y <- as.character(iris$Species[idx])
  run <- function(dv) NNS.stack(x, dv, IVs.test = iris[141:150, 1:4], method = 1,
                                folds = 2, status = FALSE, ncores = 1, seed = 1)
  chr <- run(y)
  expect_identical(chr$class.levels, c("setosa", "versicolor", "virginica"))
  expect_identical(chr, run(factor(y)))
})

test_that("NNS.stack encodes a character predictor like the equivalent factor", {
  set.seed(4)
  x <- data.frame(x = rnorm(60), g = rep(c("zeta", "alpha", "mid"), 20))
  y <- x$x + rep(c(2, 3, 0), 20) + rnorm(60)
  xf <- x
  xf$g <- factor(xf$g)
  run <- function(iv) NNS.stack(iv, y, IVs.test = iv[1:5, ], method = 2, folds = 2,
                                status = FALSE, ncores = 1, seed = 1)
  chr <- run(x)
  expect_identical(names(chr$NNS.dim.red.coefficients), c("x", "g__alpha", "g__mid", "g__zeta"))
  expect_identical(chr, run(xf))
})

balance_y <- c(1, 1, 1, 2, 2, 2, 3)

test_that("NNS.stack balancing keeps a single-member class", {
  set.seed(9)
  x <- data.frame(a = rnorm(60), b = rnorm(60))
  frame <- capture_frame("NNS.stack", NNS::NNS.stack(
    x, factor(rep(1:3, c(29, 29, 2))), IVs.test = x[1:3, ], method = 1,
    balance = TRUE, folds = 2, seed = 1, ncores = 1, status = FALSE
  ))
  set.seed(1)
  idx <- frame$.balance_indices(balance_y)
  expect_identical(as.integer(table(balance_y[idx])), c(4L, 4L, 4L))
})

test_that("NNS.boost balancing keeps a single-member class", {
  set.seed(9)
  x <- data.frame(a = rnorm(60), b = rnorm(60))
  frame <- capture_frame("NNS.boost", NNS::NNS.boost(
    x, factor(rep(1:3, c(29, 29, 2))), IVs.test = x[1:3, ], type = "CLASS",
    balance = TRUE, epochs = 0, learner.trials = 3, feature.importance = FALSE,
    status = FALSE, seed = 1
  ))
  set.seed(1)
  z <- frame$.balance_training(data.frame(i = seq_along(balance_y)), balance_y)
  expect_identical(as.integer(table(z$y)), c(4L, 4L, 4L))
  expect_identical(balance_y[z$x$i], z$y)
})

test_that("NNS.stack scores the multivariate ALL candidate with NNS.reg(n.best = 'all')", {
  set.seed(5)
  x <- data.frame(a = rnorm(160), b = rnorm(160))
  y <- rnorm(160)
  frame <- capture_frame("NNS.stack", NNS::NNS.stack(
    x, y, IVs.test = x[1:3, ], method = 1, seed = 1, ncores = 1, status = FALSE
  ))
  oof <- rep(NA_real_, length(y))
  for (s in frame$splits) {
    design <- frame$.method1_design_for_split(s$train, s$validation)
    fit <- NNS.reg(design$train, frame$y_fit_all[s$train], point.est = design$test,
                   n.best = "all", order = frame$order, dist = frame$dist, type = "REG",
                   factor.2.dummy = FALSE, plot = FALSE, residual.plot = FALSE,
                   ncores = 1, point.only = TRUE)
    oof[s$validation] <- fit$Point.est
  }
  expect_equal(frame$candidate_raws$all, oof, tolerance = 1e-12)
  expect_equal(unname(frame$candidate_scores["all"]), sum((oof - y)^2), tolerance = 1e-12)
})

test_that("NNS.stack scores the univariate ALL candidate with the deployed NNS.reg fit", {
  set.seed(5)
  x <- data.frame(a = rnorm(160))
  y <- rnorm(160)
  frame <- capture_frame("NNS.stack", NNS::NNS.stack(
    x, y, IVs.test = x[1:3, , drop = FALSE], method = 1, seed = 1, ncores = 1, status = FALSE
  ))
  oof <- rep(NA_real_, length(y))
  for (s in frame$splits) {
    design <- frame$.method1_design_for_split(s$train, s$validation)
    fit <- NNS.reg(design$train[[1L]], frame$y_fit_all[s$train], point.est = design$test[[1L]],
                   n.best = "all", order = frame$order, dist = frame$dist, type = "REG",
                   factor.2.dummy = FALSE, plot = FALSE, residual.plot = FALSE,
                   ncores = 1, point.only = TRUE)
    oof[s$validation] <- fit$Point.est
  }
  expect_equal(frame$candidate_raws$all, oof, tolerance = 1e-12)
  expect_equal(unname(frame$candidate_scores["all"]), sum((oof - y)^2), tolerance = 1e-12)
})

test_that("boost's objective sees response values for a two-valued numeric response", {
  set.seed(1)
  x <- matrix(rnorm(240), 80, 3)
  y <- ifelse(x[, 1] + rnorm(80, sd = .5) > 0, 5, -5)
  seen <- numeric(0)
  record <- function(a) seen <<- c(seen, a)
  loss <- as.expression(bquote({
    .(record)(actual)
    sum((predicted - actual)^2)
  }))
  NNS.boost(x, y, x[1:10, ], obj.fn = loss, learner.trials = 3, epochs = 0,
            threshold = 1, CV.size = .25, seed = 1, status = FALSE,
            feature.importance = FALSE)
  expect_setequal(unique(seen), c(-5, 5))
  seen <- numeric(0)
  NNS.boost(x, ifelse(y > 0, .2, -.1), x[1:10, ], obj.fn = loss, learner.trials = 3,
            epochs = 0, threshold = 1, CV.size = .25, seed = 1, status = FALSE,
            feature.importance = FALSE)
  expect_true(all(seen %in% c(-.1, .2)))
})