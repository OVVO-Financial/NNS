test_that("local derivative stencils stay in support at boundaries and interior points", {
  x <- seq(-2, 3, length.out = 101)
  for (at in c(-2, -1.99, 0, 2.99, 3)) for (width in c(.03, .1, .25)) {
    s <- NNS:::.nns_dy_stencil(x, at, width)
    expect_true(all(s$nodes >= min(x) & s$nodes <= max(x)))
    # These are derivative identities for a polynomial, including unequal
    # spacing and one-sided endpoint formulas, not comparisons to the code.
    expect_equal(sum(s$first * (s$nodes^2 + 3 * s$nodes + 4)), 2 * at + 3,
                 tolerance = 1e-9)
    expect_equal(sum(s$second * (s$nodes^2 + 3 * s$nodes + 4)), 2,
                 tolerance = 1e-8)
  }
})

test_that("discrete and constant coordinates have explicit derivative conventions", {
  two <- NNS:::.nns_dy_stencil(c(0, 0, 1, 1), .5, .1)
  expect_equal(sum(two$first * c(1, 4)), 3)
  expect_true(all(is.na(two$second)))
  constant <- NNS:::.nns_dy_stencil(rep(2, 10), 2, .1)
  expect_true(is.na(constant$first))
  expect_error(NNS:::.nns_dy_stencil(1:5, 6, .1), "inside the observed")
})

test_that("multivariate derivative boundaries remain finite and reject extrapolation", {
  x <- as.matrix(expand.grid(a = seq(-1, 1, length.out = 15),
                            b = seq(-1, 1, length.out = 15)))
  y <- 2 * x[, 1] + x[, 2]^2
  z <- dy.d_(x, y, wrt = 1:2, eval.points = rbind(c(-1, 0), c(1, 0)), messages = FALSE)
  expect_true(all(is.finite(unlist(z))))
  expect_error(dy.d_(x, y, 1, eval.points = matrix(c(0, 2), 1), messages = FALSE),
               "inside the observed")
})
