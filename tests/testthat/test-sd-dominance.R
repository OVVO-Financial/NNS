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
  expect_equal(NNS.FSD(c(0, 0), c(0, 3), "continuous", plot = FALSE), "Y OMEGA X")
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
  expect_equal(NNS.FSD(x, y, "continuous", plot = FALSE), "NO OMEGA DOMINANCE EXISTS")
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
  expect_equal(NNS.FSD(x, y, "continuous", plot = FALSE), "X OMEGA Y")
  expect_equal(c(NNS.FSD.uni(x, y, "continuous"), NNS.FSD.uni(y, x, "continuous")), c(1, 0))
  expect_equal(NNS.SD.efficient.set(cbind(x, y), 1, "continuous", status = FALSE), "x")
  expect_equal(NNS.FSD(c(1, 1), c(0, 1), "continuous", plot = FALSE), "X OMEGA Y")
  expect_equal(c(NNS.FSD.uni(c(1, 1), c(0, 1), "continuous"), NNS.FSD.uni(c(0, 1), c(1, 1), "continuous")), c(1, 0))
})

test_that("adding a constant to both samples keeps the TSD result", {
  x <- c(-.0085, -.0152, .0068); y <- c(-.0021, -.0033, -.0167)
  expect_equal(NNS.TSD.uni(x + 20000, y + 20000), NNS.TSD.uni(x, y))
  x <- c(.01, .07); y <- c(0, .08)
  expect_equal(NNS.TSD.uni(x + 1e6, y + 1e6), 1)
  expect_equal(NNS.SD.efficient.set(cbind(x, y) + 1e6, degree = 3, status = FALSE), "x")
  x <- c(.01, .03); y <- c(0, .04)
  expect_equal(c(NNS.SSD.uni(x + 20000, y + 20000), NNS.TSD.uni(x + 20000, y + 20000)), c(1, 1))
  x <- c(0, .01, .03, 1:996, 1000); y <- c(0, .02, .02, 1:996, 999)
  for (s in c(20000, 1e6)) {
    expect_equal(NNS.TSD.uni(x + s, y + s), 0)
    expect_equal(NNS.TSD.uni(c(0, .02, .02, 1:996, 1000) + s, x + s), 1)
  }
})

test_that("an infinite value is an error for SSD, TSD and continuous FSD", {
  a <- c(0, 1, Inf); b <- c(0, 1, 2)
  expect_error(NNS.SSD.uni(a, b), "infinite values")
  expect_error(NNS.TSD.uni(b, -a), "infinite values")
  expect_error(NNS.FSD.uni(a, b, "continuous"), "infinite values")
  expect_equal(c(NNS.FSD.uni(a, b), NNS.FSD.uni(b, a)), c(1, 0))
})

test_that("the efficient set does not depend on column order", {
  x <- c(-1, -1, 2); y <- c(-2, 1, 1)
  expect_equal(NNS.FSD.uni(x, y, "continuous"), 1)
  expect_identical(NNS.SD.efficient.set(cbind(Y = y, X = x), 1, "continuous", status = FALSE), "X")
  expect_identical(NNS.SD.efficient.set(cbind(X = x, Y = y), 1, "continuous", status = FALSE), "X")
  expect_identical(NNS.SD.cluster(cbind(Y = y, X = x), 1, "continuous")$Clusters,
                   NNS.SD.cluster(cbind(X = x, Y = y), 1, "continuous")$Clusters)

  a <- 0:5; h <- 3 * 2^-44
  M <- cbind(A = a, B = a + h * c(0, -2, 3, 0, -1, 0), C = a + h * c(0, -1, 0, 3, -2, 0))
  D <- NNS:::sd_dom_matrix_prefix_parallel(M, 2L, "discrete")
  expect_identical(D[cbind(c(1, 2, 3), c(2, 3, 1))], rep(1L, 3))
  for (p in list(1:3, c(2, 3, 1), c(3, 2, 1), c(2, 1, 3)))
    expect_setequal(NNS.SD.efficient.set(M[, p], 2, status = FALSE), c("A", "B", "C"))

  differ <- character(0)
  for (s in 1:200) {
    set.seed(s)
    nc <- sample(4:6, 1)
    M <- matrix(sample(-3:3, 5 * nc, replace = TRUE), 5, nc, dimnames = list(NULL, paste0("C", 1:nc)))
    p <- sample(nc)
    for (d in 1:3) for (type in if (d == 1) c("discrete", "continuous") else "discrete") {
      if (!setequal(NNS.SD.efficient.set(M[, p], d, type, status = FALSE),
                    NNS.SD.efficient.set(M, d, type, status = FALSE)))
        differ <- c(differ, paste(s, d, type))
    }
  }
  expect_identical(differ, character(0))
})

test_that("the efficient set and the clusters reject infinite values, except for discrete FSD", {
  healthy <- diff(log(c(100, 102, 99, 101, 100, 97, 99, 101)))
  bust <- diff(log(c(100, 80, 60, 40, 20, 10, 5, 0)))
  M <- cbind(healthy, bust)
  for (d in 2:3) expect_error(NNS.SD.efficient.set(M, d, status = FALSE), "infinite values")
  expect_error(NNS.SD.efficient.set(M, 1, "continuous", status = FALSE), "infinite values")
  expect_error(NNS.SD.cluster(M, 2), "infinite values")
  expect_equal(NNS.SD.efficient.set(M, 1, status = FALSE), "healthy")
  expect_equal(NNS.SD.cluster(M, 1)$Clusters, list(Cluster_1 = "healthy", Cluster_2 = "bust"))
})

test_that("NNS.TSD plots the curves against the sample values", {
  pdf(NULL); on.exit(dev.off())
  NNS.TSD(c(10, 11, 12), c(10.5, 11.5, 13))
  expect_gt(par("usr")[1], 9)
})
