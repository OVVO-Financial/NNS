historical_dep_11_6_5 <- function(x, y, asym = FALSE) {
  l <- length(x)
  obs <- as.integer(max(8, l / 8))

  part_xy <- suppressWarnings(
    NNS.part(x, y, order = NULL, obs.req = obs, min.obs.stop = FALSE,
             type = "XONLY", Voronoi = FALSE)
  )
  part_yx <- suppressWarnings(
    NNS.part(y, x, order = NULL, obs.req = obs, min.obs.stop = FALSE,
             type = "XONLY", Voronoi = FALSE)
  )

  if (nrow(part_xy$regression.points) == 0L) {
    return(list(Correlation = 0, Dependence = 0))
  }

  part_xy <- part_xy$dt
  part_xy <- part_xy[complete.cases(part_xy), , drop = FALSE]
  part_yx <- part_yx$dt
  part_yx <- part_yx[complete.cases(part_yx), , drop = FALSE]

  dep_fn <- function(xx, yy) {
    NNS.copula(cbind(xx, yy)) * sign(NNS:::fast_lm(xx, yy)$coef[2L])
  }

  grouped <- function(part) {
    group_name <- if ("quadrant" %in% names(part)) {
      "quadrant"
    } else {
      "prior.quadrant"
    }
    ids <- part[[group_name]]
    groups <- unique(ids)
    values <- vapply(groups, function(id) {
      idx <- which(ids == id)
      dep_fn(part[idx, 1L], part[idx, 2L])
    }, numeric(1L))
    weights <- vapply(groups, function(id) {
      sum(ids == id) / l
    }, numeric(1L))
    list(values = values, weights = weights)
  }

  grouped_xy <- suppressWarnings(grouped(part_xy))
  grouped_yx <- suppressWarnings(grouped(part_yx))

  global_dep <- dep_fn(x, y)
  grouped_xy$values[is.na(grouped_xy$values)] <- global_dep
  grouped_yx$values[is.na(grouped_yx$values)] <- global_dep

  dep_xy <- sum(abs(grouped_xy$values) * grouped_xy$weights)
  dep_yx <- sum(abs(grouped_yx$values) * grouped_yx$weights)
  dependence <- if (asym) dep_xy else max(c(dep_yx, dep_xy))

  lx <- length(unique(part_xy[[1L]]))
  ly <- length(unique(part_xy[[2L]]))
  degree_x <- min(10, max(1, lx - 1), max(1, ly - 1))

  if ((lx < sqrt(l)) * (ly < sqrt(l)) == 1) {
    poly_base <- suppressWarnings(tryCatch(
      NNS:::fast_lm_mult(poly(x, degree_x), abs(y))$r.squared,
      warning = function(w) dependence,
      error = function(e) dependence
    ))

    dependence <- NNS:::gravity(c(
      dependence,
      NNS.copula(cbind(x, y), plot = FALSE),
      poly_base
    ))
  }

  corr_xy <- sum(grouped_xy$values * grouped_xy$weights)
  corr_yx <- sum(grouped_yx$values * grouped_yx$weights)
  correlation <- if (asym) corr_xy else max(c(corr_yx, corr_xy))

  list(Correlation = correlation, Dependence = dependence)
}

expect_historical_dep <- function(x, y, asym = FALSE, tolerance = 1e-12) {
  expected <- historical_dep_11_6_5(x, y, asym = asym)
  actual <- NNS.dep(x, y, asym = asym, print.map = FALSE)

  expect_equal(actual$Correlation, expected$Correlation, tolerance = tolerance)
  expect_equal(actual$Dependence, expected$Dependence, tolerance = tolerance)
  expect_gte(actual$Dependence, 0)
  expect_lte(actual$Dependence, 1)
}

test_that("NNS.dep reproduces NNS 11.6.5 for independent Gaussian samples", {
  set.seed(123)
  x <- rnorm(100)
  y <- rnorm(100)

  expect_historical_dep(x, y)
  expect_lt(NNS.dep(x, y)$Dependence, 0.4)
})

test_that("NNS.dep reproduces NNS 11.6.5 for nonlinear continuous data", {
  set.seed(42)
  x <- seq(-3, 3, length.out = 160)
  y <- x^2 + 0.15 * rnorm(length(x))

  expect_historical_dep(x, y)
})

test_that("NNS.dep reproduces the NNS 11.6.5 discrete fallback", {
  x <- rep(1:4, each = 30)
  y <- rep(c(1, 1, 2, 2), each = 30)

  expect_historical_dep(x, y)
})

test_that("NNS.dep preserves NNS 11.6.5 directional semantics", {
  set.seed(7)
  x <- runif(140, -2, 2)
  y <- exp(x) + rnorm(140, sd = 0.2)

  expect_historical_dep(x, y, asym = TRUE)
  expect_historical_dep(y, x, asym = TRUE)
})

test_that("NNS.dep matrix entries use the restored pair estimator", {
  set.seed(99)
  x <- rnorm(90)
  y <- x^2 + rnorm(90, sd = 0.3)
  z <- rnorm(90)
  m <- cbind(x = x, y = y, z = z)

  result <- NNS.dep(m)
  xy <- NNS.dep(x, y)
  xz <- NNS.dep(x, z)
  yz <- NNS.dep(y, z)

  expect_equal(result$Correlation["x", "y"], xy$Correlation)
  expect_equal(result$Dependence["x", "y"], xy$Dependence)
  expect_equal(result$Correlation["x", "z"], xz$Correlation)
  expect_equal(result$Dependence["x", "z"], xz$Dependence)
  expect_equal(result$Correlation["y", "z"], yz$Correlation)
  expect_equal(result$Dependence["y", "z"], yz$Dependence)
})
