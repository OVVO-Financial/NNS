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
  set.seed(17)
  x <- rexp(60)
  o <- NNS.ARMA.optim(x, h = 6, training.set = 56, seasonal.factor = c(2, 4),
                      ncores = 1, print.trace = FALSE, plot = FALSE)
  expect_true(any(o$results == 0))
  expect_true(all(o$upper.pred.int >= o$results))
  expect_true(all(o$upper.pred.int[o$results == 0] > 0))
})
