# Follow-ups from the review of the completion branch (#150).

test_that("stack keeps the same method 1 columns when a fold has a constant predictor", {
  # b is constant in the first chronological training fold, so its fold
  # coefficient is 0. Dropping it gave that fold a univariate design with no
  # neighbour counts to tune, and the stack stopped.
  x <- data.frame(a = rep(1:10, 6), b = c(rep(0, 40), rep(c(-10, 10), 10)))
  y <- x$a + 2 * x$b
  fit <- NNS.stack(x, y, x[51:60, ], ts.test = 10, folds = 2, order = "max",
                   dim.red.method = "cor", ncores = 1, status = FALSE)
  expect_length(fit$stack, 10)
  expect_true(all(is.finite(fit$stack)))
  expect_lt(max(abs(fit$stack - y[51:60])), 2)
})

test_that("a declared factor level with no training rows is a valid test value", {
  x <- data.frame(g = factor(c(rep(c("a", "b"), 19), "a", "rare")))
  y <- factor(rep(c("L", "R"), 20))
  fit <- NNS.boost(x, y, ts.test = 8, learner.trials = 1, epochs = 0,
                   features.only = TRUE, feature.importance = FALSE,
                   status = FALSE, seed = 9)
  expect_equal(unname(unlist(fit$feature.weights)), 1)

  g <- factor(c("a", "b", "a", "b"), levels = c("a", "b", "c"))
  at <- function(level) NNS.reg(data.frame(g = g), c(1, 2, 1, 2),
    point.est = data.frame(g = factor(level, levels = levels(g))),
    plot = FALSE, residual.plot = FALSE)$Point.est
  expect_equal(at("a"), 1)
  expect_equal(at("b"), 2)
  expect_true(is.finite(at("c")))
  # Level order still does not matter.
  g2 <- factor(g, levels = c("c", "b", "a"))
  expect_equal(NNS.reg(data.frame(g = g2), c(1, 2, 1, 2),
    point.est = data.frame(g = factor("c", levels = levels(g2))),
    plot = FALSE, residual.plot = FALSE)$Point.est, at("c"))

  set.seed(1)
  gg <- factor(sample(c("a", "b"), 60, TRUE), levels = c("a", "b", "c"))
  zz <- rnorm(60)
  s <- NNS.stack(data.frame(g = gg, z = zz), zz + (gg == "b"),
    data.frame(g = factor(c("a", "c"), levels = levels(gg)), z = c(0, 0)),
    folds = 2, ncores = 1, status = FALSE)
  expect_true(all(is.finite(s$stack)))
  # An undeclared character value is still unseen.
  expect_error(NNS.reg(data.frame(g = c("a", "b", "a", "b")), c(1, 2, 1, 2),
    point.est = data.frame(g = "c"), plot = FALSE), "unseen level")
})

test_that("class probabilities over a neighbour-count path equal separate fits", {
  set.seed(11)
  X <- cbind(u = round(runif(120), 1), v = round(runif(120), 1))  # tied distances
  cls <- sample(1:3, 120, TRUE)
  P <- outer(cls, 1:3, "==") * 1
  Q <- cbind(u = round(runif(25), 1), v = round(runif(25), 1))
  ks <- list(1L, 2L, 3L, 7L, "all")
  path <- NNS:::.nns_class_base(X, P, Q, NULL, NULL, "NNS", 1L,
                                need.fitted = FALSE, k.path = ks)$predicted.path
  for (a in seq_along(ks)) {
    one <- NNS:::.nns_class_base(X, P, Q, NULL, ks[[a]], "NNS", 1L)$predicted
    expect_identical(path[[a]], one)
  }
  # Point-only classification skips the fitted classes but not the predictions.
  full <- NNS.reg(X, factor(cls), point.est = Q, plot = FALSE, residual.plot = FALSE)
  po <- NNS.reg(X, factor(cls), point.est = Q, plot = FALSE, residual.plot = FALSE,
                point.only = TRUE)
  expect_identical(po$Point.est, full$Point.est)
  expect_identical(po$Point.probabilities, full$Point.probabilities)
  expect_null(po$Fitted.xy)
})

test_that("stack accepts a single test row with a factor predictor", {
  # vapply() returned the one-row dummy block as a vector, which was then
  # reshaped into one column per row instead of one row.
  g <- factor(rep(c("a", "b"), 20))
  z <- seq(0, 1, length.out = 40)
  fit <- NNS.stack(data.frame(g = g, z = z), z + (g == "b"),
    data.frame(g = factor("b", levels = c("a", "b")), z = .5),
    folds = 2, ncores = 1, status = FALSE)
  expect_length(fit$stack, 1)
  expect_true(abs(fit$stack - 1.5) < .25)
})

test_that("categorical NNS.reg honours an explicit n.best", {
  x <- data.frame(a = factor(c("a", "a", "b")), b = factor(c("u", "v", "u")))
  q <- data.frame(a = "b", b = "v")   # an unobserved joint combination
  est <- function(nb) NNS.reg(x, c(0, 10, 20), point.est = q, n.best = nb,
    noise.reduction = "mean", plot = FALSE, residual.plot = FALSE)$Point.est
  expect_equal(est(NULL), 15)         # the two nearest cells, tied
  expect_equal(est("all"), est(3))
  expect_false(isTRUE(all.equal(est("all"), est(NULL))))
})

test_that("dy.d_ accepts a one-column data frame of evaluation points", {
  set.seed(1)
  x <- cbind(a = runif(200, -1, 1), b = runif(200, -1, 1))
  y <- x[, 1]^2 + x[, 2]
  a <- dy.d_(x, y, 1, eval.points = data.frame(a = c(-.5, .5)), messages = FALSE)
  b <- dy.d_(x, y, 1, eval.points = c(-.5, .5), messages = FALSE)
  expect_equal(a, b)
})

test_that("bulk distance wrappers match an unnamed matrix by position", {
  rpm <- data.frame(X1 = c(0, 1, 0, 1), X2 = c(0, 0, 1, 1), y.hat = c(1, 2, 3, 4))
  q <- matrix(c(.2, .9, .1, .8), 2)
  expect_equal(NNS:::NNS.distance.bulk(rpm, q, 2),
               NNS:::NNS.distance.bulk(rpm, setNames(as.data.frame(q), c("X1", "X2")), 2))
  expect_equal(NNS:::NNS.distance.path.bulk(rpm, q, 3),
               NNS:::NNS.distance.path.bulk(rpm, setNames(as.data.frame(q), c("X1", "X2")), 3))
})

test_that("NNS.ARMA.optim adds no bias shift by default", {
  expect_identical(eval(formals(NNS.ARMA.optim)$bias.correction), c("none", "split"))
  # The same-window rule (estimate and accept on one window) is not offered.
  expect_error(NNS.ARMA.optim(as.numeric(AirPassengers[1:96]), h = 12, seasonal.factor = 12,
    print.trace = FALSE, ncores = 1, pred.int = NULL, bias.correction = "legacy"))
  fit <- NNS.ARMA.optim(as.numeric(AirPassengers[1:96]), h = 12,
    seasonal.factor = c(6, 12), print.trace = FALSE, ncores = 1, pred.int = NULL)
  expect_identical(fit$bias.shift, 0)
  expect_identical(fit$bias.validation$policy, "none")
  expect_equal(fit$results, fit$raw.results)
})
