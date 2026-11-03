test_that("a shifted binary pair keeps perfect dependence", {
  b <- rep(0:1, each = 50)
  expect_equal(unlist(NNS.dep(b + .1, b + .1)), c(Correlation = 1, Dependence = 1))
})

test_that("y determined by a decimal-valued discrete x has dependence 1", {
  x <- rep(c(.1, .2, .3, .7), 25)
  expect_equal(NNS.dep(x, x^2)$Dependence, 1)
})

test_that("a modulo period beyond the n/5 scan keeps its ARMA weight", {
  set.seed(1)
  x <- 10 + sin(2 * pi * (1:96) / 24) + rnorm(96, sd = .1)
  s <- NNS.seas(x, modulo = 24, mod.only = FALSE, plot = FALSE)
  w <- NNS:::ARMA.seas.weighting(FALSE, s$all.periods)
  expect_gt(w$Weights[match(24, w$lag)], 0)
  expect_lt(s$all.periods$Coefficient.of.Variation[s$all.periods$Period == 24],
            s$all.periods$Variable.Coefficient.of.Variation[1])
})

test_that("ARMA.optim out-of-sample bands sit on the clamped forecasts", {
  # A fixed positive bias in the validation prefix persists in its acceptance
  # half. The same shift must clamp the lower out-of-sample forecast to zero.
  local_mocked_bindings(
    NNS.ARMA = function(variable, h, ...) rep(c(1, 9), length.out = h),
    .package = "NNS"
  )
  x <- rep(1, 80)
  o <- NNS.ARMA.optim(x, h = 6, training.set = 64, seasonal.factor = c(2, 4),
                      obj.fn = expression(mean((predicted - actual)^2)),
                      ncores = 1, print.trace = FALSE, plot = FALSE)
  expect_true(o$bias.validation$accepted)
  expect_true(any(o$results == 0))
  expect_true(all(o$upper.pred.int >= o$results))
  expect_true(all(o$upper.pred.int[o$results == 0] > 0))
  expect_equal(o$results, pmax(0, o$raw.results + o$bias.shift))
})

test_that("constant cells keep the direction of a discrete relationship", {
  b <- rep(0:1, each = 50)
  expect_equal(unlist(NNS.dep(b, -b)), c(Correlation = -1, Dependence = 1))
  expect_equal(NNS.dep(b + .1, -(b + .1))$Correlation, -1)
  u <- rep(1:3, 30)
  expect_equal(NNS.dep(u, (u - 2)^2)$Correlation, 0)
  g <- rep(1:3, each = 30)
  expect_equal(NNS.dep(g, (g - 2)^2)$Correlation, 0)
  expect_equal(NNS.dep(rev(g), rev((g - 2)^2))$Correlation, 0)
  x <- 1e6 + rep(c(-.01, 0, .01), each = 1000); y <- 1e6 + rep(c(.01, 0, .01), each = 1000)
  expect_equal(c(NNS.dep(x, y)$Correlation, NNS.dep(rev(x), rev(y))$Correlation), c(0, 0))
  expect_equal(NNS.dep(u, u)$Correlation, 1)
  expect_equal(NNS.dep(rep(c(0, 1, 2, 10), 30), rep(0:3, 30))$Correlation, 1)
})
