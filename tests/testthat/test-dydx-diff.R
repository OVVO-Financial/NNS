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

test_that("dy.dx accepts only the overall option for text queries", {
  x <- seq(0, 6, length.out = 70)
  expect_error(dy.dx(x, sin(x), "1.75"), "numeric.*overall")
  expect_error(dy.dx(x, sin(x), c("overall", "typo")), "numeric.*overall")
  fit <- NNS.reg(x, sin(x), plot = FALSE, ncores = 1)
  expect_equal(dy.dx(x, sin(x), "overall"), mean(fit$Fitted.xy$gradient))
})

test_that("dy.dx returns both derivatives in query order, including repeats", {
  x <- seq(0, 6, length.out = 70)
  points <- c(3, 1, 3, 2)
  expected <- dy.dx(x, sin(x), 1:3)
  actual <- dy.dx(x, sin(x), points)
  expect_equal(actual$eval.point, points)
  expect_equal(actual$first.derivative, expected$first.derivative[points])
  expect_equal(actual$second.derivative, expected$second.derivative[points])
})

test_that("dy.dx uses the extrapolated fit outside the observed range", {
  x <- seq(0, 2 * pi, pi / 100)
  points <- c(1, 2 * pi + 0.025, 8)
  d <- dy.dx(x, sin(x), points)
  one <- lapply(points, function(p) dy.dx(x, sin(x), p))
  expect_true(all(is.finite(d$first.derivative)))
  expect_equal(d$first.derivative, vapply(one, `[[`, numeric(1), "first.derivative"))
  expect_equal(d$second.derivative, vapply(one, `[[`, numeric(1), "second.derivative"))
})
