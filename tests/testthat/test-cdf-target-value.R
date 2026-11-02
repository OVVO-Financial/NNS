cdf_types <- c("CDF", "survival", "hazard", "cumulative hazard")

cdf_target_data <- function() {
  set.seed(11)
  list(grid = 1:20, skewed = sort(round(rlnorm(30), 2)))
}

cdf_at <- function(x, degree, type, target) {
  NNS.CDF(x, degree, target = target, type = type, plot = FALSE)$target.value
}

test_that("univariate target.value equals the curve at each data point", {
  for (x in cdf_target_data()) for (degree in 0:1) for (type in cdf_types) {
    curve <- NNS.CDF(x, degree, type = type, plot = FALSE)$Function[[2]]
    at <- vapply(sort(x), function(t) cdf_at(x, degree, type, t), numeric(1))
    expect_equal(at, curve, tolerance = 1e-12, info = paste(type, degree))
  }
})

test_that("univariate target.value lies between the neighbouring curve values", {
  for (x in cdf_target_data()) for (degree in 0:1) for (type in cdf_types[-3]) {
    xs <- sort(unique(x))
    curve <- NNS.CDF(x, degree, type = type, plot = FALSE)$Function[[2]]
    curve <- curve[!duplicated(sort(x))]
    mid <- (xs[-1] + xs[-length(xs)]) / 2
    at <- vapply(mid, function(t) cdf_at(x, degree, type, t), numeric(1))
    lo <- pmin(curve[-1], curve[-length(curve)])
    hi <- pmax(curve[-1], curve[-length(curve)])
    expect_true(all(at >= lo - 1e-12 & at <= hi + 1e-12), info = paste(type, degree))
  }
})

test_that("hazard target.value between data points is finite and clamped like the curve", {
  for (x in cdf_target_data()) for (degree in 0:1) {
    xs <- sort(unique(x))
    mid <- (xs[-1] + xs[-length(xs)]) / 2
    at <- vapply(mid, function(t) cdf_at(x, degree, "hazard", t), numeric(1))
    expect_true(all(is.finite(at) & at >= 0 & at <= 1e6))
  }
})

test_that("cumulative hazard at a target is -log of the survival there", {
  expect_equal(cdf_at(1:20, 0, "cumulative hazard", 19), -log(.05))
  expect_equal(cdf_at(1:20, 0, "cumulative hazard", 10), -log(.5))
  h <- NNS.CDF(1:20, 0, target = 20, type = "cumulative hazard", plot = FALSE)
  expect_true(is.finite(h$target.value))
  expect_equal(h$target.value, h$Function[[2]][20])
})
