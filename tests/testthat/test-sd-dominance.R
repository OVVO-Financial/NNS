test_that(".uni tests walk samples of different lengths", {
  expect_equal(NNS.FSD.uni(4, c(3, 0)), 1)
  expect_equal(NNS.SSD.uni(c(1, 4), c(3, 1, 3, 0, 2)), 1)
  expect_equal(NNS.TSD.uni(c(4, 1), c(0, 3, 3, 4)), 1)
  expect_equal(NNS.FSD(4, c(3, 0), plot = FALSE), "X FSD Y")
  expect_equal(NNS.SSD(c(1, 4), c(3, 1, 3, 0, 2), plot = FALSE), "X SSD Y")
  expect_equal(NNS.TSD(c(4, 1), c(0, 3, 3, 4), plot = FALSE), "X TSD Y")
})

test_that(".uni tests report no dominance when the curves coincide", {
  expect_equal(NNS.FSD.uni(c(1, 2), c(1, 1, 2, 2), "continuous"), 0)
  expect_equal(NNS.SSD.uni(c(1, 2), c(1, 1, 2, 2)), 0)
  expect_equal(NNS.TSD.uni(c(1, 2), c(1, 1, 2, 2)), 0)
})

test_that("continuous FSD ratios that agree at the knots but cross are not dominance", {
  x <- c(0, 0, 3, 3); y <- c(1, 1, 3, 0)
  expect_equal(NNS.FSD.uni(x, y, "continuous"), 0)
  expect_equal(NNS.FSD.uni(y, x, "continuous"), 0)
  expect_setequal(NNS.SD.efficient.set(cbind(x, y), 1, "continuous", status = FALSE), c("x", "y"))
})

test_that("continuous FSD uses CDF 1 when every observation equals the threshold", {
  expect_equal(NNS.FSD.uni(rep(4, 6), c(5, 5, 0, 4, 0, 2), "continuous"), 0)
  expect_equal(NNS.FSD.uni(c(0, 3), c(0, 0), "continuous"), 1)
  expect_equal(NNS.FSD(c(0, 0), c(0, 3), "continuous", plot = FALSE), "Y FSD X")
})

test_that("the missing-value check does not recycle samples of different lengths", {
  expect_silent(NNS.FSD(c(1, 2, 5), c(0, 1, 1, 2, 3), plot = FALSE))
  expect_silent(NNS.SSD(c(1, 2, 5), c(0, 1, 1, 2, 3), plot = FALSE))
  expect_silent(NNS.TSD(c(1, 2, 5), c(0, 1, 1, 2, 3), plot = FALSE))
  expect_silent(NNS.SSD.uni(c(1, 2, 5), c(0, 1, 1, 2, 3)))
  expect_error(NNS.SSD(c(1, NA), 1:3, plot = FALSE), "missing values")
  expect_error(NNS.TSD.uni(c(1, NaN), 1:3), "missing values")
})

test_that("TSD checks LPM2 between the sample points", {
  x <- c(-0.85, -1.52, 0.68); y <- c(-0.21, -0.33, -1.67)
  expect_equal(mean(pmax(0.16 - x, 0)^2) - mean(pmax(0.16 - y, 0)^2), 0.03886667, tolerance = 1e-6)
  expect_equal(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
  expect_equal(NNS.TSD.uni(x, y), 0)
  expect_setequal(NNS.SD.efficient.set(cbind(x, y), 3, status = FALSE), c("x", "y"))
})

test_that("continuous FSD checks the ratios between the sample points", {
  x <- c(0.15, 1.55, -0.86, -0.8, 0.25, 0.33); y <- c(-0.15, 1.17, 0.26, -0.04, -2.44, 1.19)
  ratio <- function(v, t) mean(pmax(t - v, 0)) / mean(abs(t - v))
  expect_equal(ratio(x, 0.5) - ratio(y, 0.5), 0.002972295, tolerance = 1e-6)
  expect_equal(NNS.FSD(x, y, "continuous", plot = FALSE), "NO FSD EXISTS")
  expect_equal(NNS.FSD.uni(x, y, "continuous"), 0)
  expect_setequal(NNS.SD.efficient.set(cbind(x, y), 1, "continuous", status = FALSE), c("x", "y"))
})

test_that("SSD treats rounding-level curve differences as ties", {
  expect_equal(NNS.SSD(c(0.1, 0.7), c(0, 0.8), plot = FALSE), "X SSD Y")
  expect_equal(NNS.SSD(c(1, 7), c(0, 8), plot = FALSE), "X SSD Y")
  expect_equal(NNS.SSD.uni(c(0.1, 0.7), c(0, 0.8)), 1)
  z <- qnorm(ppoints(100))
  expect_equal(NNS.SSD(z, 2 * z, plot = FALSE), "X SSD Y")
  expect_equal(NNS.SSD.uni(z + 0.37, 2 * z + 0.37), 1)
  expect_equal(NNS.SSD(3 * z + 0.37, 6 * z + 0.37, plot = FALSE), "X SSD Y")
})

test_that("a mean-preserving spread of a large sample is still dominated", {
  z <- qnorm(ppoints(1e5))
  expect_equal(NNS.SSD.uni(z + 100, 2 * z + 100), 1)
  expect_equal(NNS.TSD.uni(z + 100, 2 * z + 100), 1)
  expect_equal(NNS.SSD(z + 100, 2 * z + 100, plot = FALSE), "X SSD Y")
  expect_equal(NNS.TSD(z + 100, 2 * z + 100, plot = FALSE), "X TSD Y")
})

test_that("an empty sample is an error", {
  expect_error(NNS.FSD(numeric(0), c(0, 1), plot = FALSE), "at least one observation")
  expect_error(NNS.SSD(c(0, 1), numeric(0), plot = FALSE), "at least one observation")
  expect_error(NNS.TSD(numeric(0), c(0, 1), plot = FALSE), "at least one observation")
})

test_that("continuous FSD sees ratios that differ only between the knots", {
  x <- c(0, 1, 1); y <- c(0, 0, 1)
  expect_equal(NNS.FSD(x, y, "continuous", plot = FALSE), "X FSD Y")
  expect_equal(c(NNS.FSD.uni(x, y, "continuous"), NNS.FSD.uni(y, x, "continuous")), c(1, 0))
  expect_equal(NNS.SD.efficient.set(cbind(x, y), 1, "continuous", status = FALSE), "x")
  expect_equal(NNS.FSD(c(1, 1), c(0, 1), "continuous", plot = FALSE), "X FSD Y")
  expect_equal(c(NNS.FSD.uni(c(1, 1), c(0, 1), "continuous"), NNS.FSD.uni(c(0, 1), c(1, 1), "continuous")), c(1, 0))
})
