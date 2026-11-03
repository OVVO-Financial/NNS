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
