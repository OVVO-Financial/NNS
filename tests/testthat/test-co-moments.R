lpm_leg <- function(v, t, d) if (d == 0) as.numeric(v <= t) else pmax(t - v, 0)^d
upm_leg <- function(v, t, d) if (d == 0) as.numeric(v > t) else pmax(v - t, 0)^d

test_that("co-moments apply each degree to its own variable when one is fractional", {
  expect_equal(D.UPM(.5, 1, c(-4, -9), c(2, 3), 0, 0), 6.5)
  expect_equal(D.UPM(1, .5, c(-2, -3), c(4, 9), 0, 0), 6.5)

  set.seed(11)
  x <- rnorm(60); y <- rnorm(60)
  degrees <- c(0, .5, 1, 1.5, 2, 3)
  for (tx in c(0.1, -0.3)) for (ty in c(-0.2, 0.4)) {
    for (a in degrees) for (b in degrees) {
      expect_equal(Co.LPM(a, x, y, tx, ty, b), mean(lpm_leg(x, tx, a) * lpm_leg(y, ty, b)))
      expect_equal(Co.UPM(a, x, y, tx, ty, b), mean(upm_leg(x, tx, a) * upm_leg(y, ty, b)))
      expect_equal(D.LPM(a, b, x, y, tx, ty), mean(upm_leg(x, tx, b) * lpm_leg(y, ty, a)))
      expect_equal(D.UPM(a, b, x, y, tx, ty), mean(lpm_leg(x, tx, a) * upm_leg(y, ty, b)))
    }
  }
})

test_that("co-moments reject unequal lengths, empty input and recycled targets", {
  set.seed(1)
  x <- rnorm(100); y <- rnorm(99)
  fns <- list(
    function(x, y, tx, ty) Co.LPM(1, x, y, tx, ty),
    function(x, y, tx, ty) Co.UPM(1, x, y, tx, ty),
    function(x, y, tx, ty) D.LPM(1, 1, x, y, tx, ty),
    function(x, y, tx, ty) D.UPM(1, 1, x, y, tx, ty)
  )
  for (f in fns) {
    expect_error(f(x, y, 0, 0), "same length")
    expect_error(f(x, y, NULL, NULL), "same length")
    expect_error(f(x, y, seq(-1, 1, length.out = 200), 0), "same length")
    expect_error(f(numeric(0), numeric(0), 0, 0), "must not be empty")
    expect_error(f(numeric(0), numeric(0), NULL, NULL), "must not be empty")
    expect_error(f(1:5, 1:5, numeric(0), 0), "must not be empty")
    expect_error(f(1:5, 1:5, 0, numeric(0)), "must not be empty")
    expect_error(f(1:5, 1:5, numeric(0), numeric(0)), "must not be empty")
    expect_error(f(x, x, c(0, 1, 2), c(0, 1)), "length 1")
    expect_equal(f(x, x, c(-1, 0, 1), 0.5),
                 c(f(x, x, -1, 0.5), f(x, x, 0, 0.5), f(x, x, 1, 0.5)))
  }
})
