arma <- function(method, shrink = FALSE, h = 1) {
  NNS.ARMA(AirPassengers, h = h, seasonal.factor = 12, method = method, shrink = shrink, plot = FALSE)
}

test_that("shrink blends the nonlin estimate with the means estimate", {
  expect_equal(arma("nonlin", TRUE), (arma("nonlin") + arma("means")) / 2)
})

test_that("shrink with method = \"both\" shrinks both components", {
  expect_equal(arma("both", TRUE), (arma("lin", TRUE) + arma("nonlin", TRUE)) / 2)
})

test_that("shrink moves every nonlin forecast towards the means forecast", {
  nonlin <- arma("nonlin", h = 12)
  shrunk <- arma("nonlin", TRUE, h = 12)
  means <- arma("means", h = 12)
  expect_true(all(nonlin != means))
  expect_true(all(shrunk != nonlin))
  expect_true(all(shrunk > pmin(nonlin, means) & shrunk < pmax(nonlin, means)))
})

test_that("nonlinear shrink preserves positive affine changes of level", {
  f <- function(x) NNS.ARMA(x, h = 1, seasonal.factor = 12,
                            method = "nonlin", shrink = TRUE, plot = FALSE)
  expect_equal(f(2 * AirPassengers + 7), 2 * f(AirPassengers) + 7, tolerance = 1e-6)
})
