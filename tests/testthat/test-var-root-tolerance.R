var_root_data <- function() {
  set.seed(42)
  list(
    symmetric = c(-.003, -.002, -.001, .001, .002, .003),
    skewed = c(rexp(40) - 1, 5),
    returns = rnorm(500, 0, .01)
  )
}

test_that("integer-degree LPM.VaR finds the symmetric median at small scale", {
  x <- c(-.003, -.002, -.001, .001, .002, .003)
  expect_equal(LPM.VaR(.5, 4, x), 0, tolerance = 1e-12)
  expect_equal(UPM.VaR(.5, 4, x), 0, tolerance = 1e-12)
})

test_that("LPM.VaR and UPM.VaR scale with the data", {
  p <- c(.01, .05, .5, .95)
  for (x in var_root_data()) for (d in 1:5) {
    spread <- diff(range(x))
    lpm <- LPM.VaR(p, d, x)
    upm <- UPM.VaR(p, d, x)
    for (k in c(1e-4, 1, 1e4)) {
      expect_lt(max(abs(LPM.VaR(p, d, k * x) - k * lpm)), 1e-10 * k * spread)
      expect_lt(max(abs(UPM.VaR(p, d, k * x) - k * upm)), 1e-10 * k * spread)
    }
  }
})

test_that("LPM.VaR and UPM.VaR invert their ratios at any scale", {
  p <- c(.01, .05, .5, .95)
  for (x in var_root_data()) for (d in 1:5) for (k in c(1e-4, 1, 1e4)) {
    y <- k * x
    lr <- vapply(p, function(q) LPM.ratio(d, LPM.VaR(q, d, y), y), numeric(1))
    ur <- vapply(p, function(q) UPM.ratio(d, UPM.VaR(q, d, y), y), numeric(1))
    # Degree 5 uses the general solver, which the fractional-degree tests hold to 1e-7.
    tol <- if (d <= 4) 1e-8 else 1e-7
    expect_lt(max(abs(lr - p)), tol)
    expect_lt(max(abs(ur - p)), tol)
  }
})

test_that("integer-degree VaR inverts the ratio where it is steep", {
  x <- c(rep(0, 99999), 1)
  for (d in 1:4) {
    expect_lt(abs(LPM.ratio(d, LPM.VaR(.036, d, x), x) - .036), 1e-10)
    expect_lt(abs(UPM.ratio(d, UPM.VaR(.964, d, x), x) - .964), 1e-10)
  }
})
