test_that("pred.int works for a one-point forecast path (#98)", {
  set.seed(123)
  r <- NNS.ARMA(AirPassengers, h = 1, seasonal.factor = 12, pred.int = .95, plot = FALSE)
  expect_equal(nrow(r), 1L)
  expect_lte(r[[2]], r[[1]]); expect_gte(r[[3]], r[[1]])
  set.seed(123)
  r <- NNS.ARMA(AirPassengers, h = 1, seasonal.factor = c(12, 24), method = "lin", pred.int = .95, plot = FALSE)
  expect_lt(r[[2]], r[[3]])   # the per-lag spread widens it
})

test_that("pred.int works for a flat forecast path (#98)", {
  set.seed(123)
  r <- NNS.ARMA(rep(2, 30), h = 12, seasonal.factor = 2, method = "means", pred.int = .95, plot = FALSE)
  expect_equal(unname(unlist(r[, 1])), rep(2, 12))
  expect_equal(unname(unlist(r[, 2])), rep(2, 12))
  expect_equal(unname(unlist(r[, 3])), rep(2, 12))
})
