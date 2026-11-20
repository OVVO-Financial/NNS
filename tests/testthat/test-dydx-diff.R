test_that("dy.dx stops with a clear message when eval.point is missing", {
  x <- seq(-2, 2, length.out = 101)
  for (ep in list(NULL, numeric(0))) {
    expect_error(dy.dx(x, x^2, eval.point = ep),
                 "eval.point must be provided", fixed = TRUE)
  }
  expect_error(dy.dx(x, x^2), "eval.point must be provided", fixed = TRUE)
})

test_that("dy.dx does not depend on the row order of x", {
  set.seed(17)
  x <- runif(300, -2, 2)
  y <- x^3
  o <- order(x)
  ep <- c(-1.5, 0.5, 1.5)
  expect_equal(dy.dx(x, y, ep), dy.dx(x[o], y[o], ep))
  expect_equal(dy.dx(x, y, ep)$first.derivative, 3 * ep^2, tolerance = 0.1)
})

test_that("dy.d_ gives an irrelevant balanced regressor a slope near zero", {
  # A complete Cartesian design avoids incidental sample association and
  # sparse-neighbour variation at the evaluated tuples.
  x <- as.matrix(expand.grid(seq(-2, 2, length.out = 21),
                            seq(-2, 2, length.out = 21)))
  ep <- rbind(c(-1, 0), c(0, 1), c(1, -1))
  r <- dy.d_(x, 2 * x[, 1], wrt = 1:2, eval.points = ep, messages = FALSE)
  expect_true(all(r["First", 1][[1]] > 1))
  expect_true(all(abs(r["First", 2][[1]]) < 0.1))
})

test_that("NNS.diff converges to smooth derivatives within 1e-6", {
  cases <- list(list(log, 5, 0.2),
                list(tanh, 0.5, 1 - tanh(0.5)^2),
                list(function(x) x^4 - 3 * x, -2, -35))
  for (cs in cases) {
    r <- NNS.diff(cs[[1]], cs[[2]], digits = 16)
    expect_equal(unname(r["converged", 1]), 1)
    # Secant-projection root finding is numerically platform dependent; the
    # convergence flag does not certify 1e-7 absolute derivative accuracy.
    # Retain a strict absolute numerical check aligned with the neighboring
    # constant-offset test, including on CRAN macOS R old-release.
    expect_lt(abs(r["DERIVATIVE", 1] - cs[[3]]), 1e-6)
  }
})

test_that("NNS.diff is unaffected by a constant added to f", {
  d <- function(k, ...) NNS.diff(function(x) sin(x) + k, 1, digits = 16, ...)
  for (k in c(0, 1000)) {
    r <- d(k)
    expect_equal(unname(r["converged", 1]), 1)
    expect_lt(abs(r["DERIVATIVE", 1] - cos(1)), 1e-6)
  }
  expect_equal(d(1000, tol = 1e-4)["DERIVATIVE", 1],
               d(0, tol = 1e-4)["DERIVATIVE", 1], tolerance = 1e-8)
})

test_that("NNS.diff refines a step when projected intercepts coincide", {
  origin <- NNS.diff(exp, 0)
  expect_equal(unname(origin["converged", 1]), 1)
  expect_gt(origin["iterations", 1], 0)
  expect_lt(abs(origin["DERIVATIVE", 1] - 1), 1e-6)

  inflection <- NNS.diff(sin, 2 * pi)
  expect_equal(unname(inflection["converged", 1]), 1)
  expect_lt(abs(inflection["DERIVATIVE", 1] - 1), 1e-6)
})

test_that("NNS.diff respects the iteration limit when the intercepts coincide", {
  limited <- NNS.diff(function(x) (x - 1)^3, 1, max.iter = 1)
  expect_equal(unname(limited["converged", 1]), 0)
  expect_equal(unname(limited["termination.code", 1]), 1)
  expect_equal(unname(limited["Inferred h", 1]), 0.055)
})

test_that("NNS.diff flags unresolved one-sided slopes and roundoff at zero", {
  cusp <- NNS.diff(abs, 0, max.iter = 5)
  expect_equal(unname(cusp["converged", 1]), 0)
  offset <- NNS.diff(function(x) 1e12 + exp(x), 0)
  expect_equal(unname(offset["converged", 1]), 0)
  expect_equal(unname(offset["termination.code", 1]), 2)
})

test_that("NNS.diff refines nearly coincident projections without root finding", {
  x <- -2.10611725906715e-8
  expect_no_warning(result <- NNS.diff(function(t) exp(-t^2), x))
  expect_equal(unname(result["converged", 1]), 1)
  expect_lt(abs(result["DERIVATIVE", 1] + 2 * x * exp(-x^2)), 1e-10)
})

test_that("NNS.diff resolves linear and smooth slopes with an offset", {
  line <- NNS.diff(function(x) 1 + x, 0, h = 2^-18)
  expect_identical(unname(line[c("DERIVATIVE", "termination.code"), 1]), c(1, 0))
  expect_gt(line["iterations", 1], 0)
  shifted <- NNS.diff(function(x) 10 + exp(x), 0)
  expect_equal(unname(shifted["converged", 1]), 1)
  expect_lt(abs(shifted["DERIVATIVE", 1] - 1), 1e-10)
})
