# The compiled kernels that replaced interpreted loops reproduce the R
# definitions they were written from.

test_that("the compiled LPM.VaR solver matches the R implementation", {
  set.seed(101)
  for (i in 1:40) {
    n <- sample(c(1, 2, 5, 37, 300), 1)
    x <- switch(sample(4, 1), rnorm(n), round(rnorm(n), 1), rnorm(n) + 1e6,
                sample(1:4, n, TRUE) + 0)
    p <- c(runif(6), 0, 1, .5)
    for (degree in 1:4) {
      expect_equal(NNS_LPM_VaR_integer_cpp(p, degree, x),
                   .NNS_LPM_VaR_integer_R(p, degree, x), tolerance = 1e-12)
    }
  }
  # Values the kernel leaves to the R implementation.
  expect_null(NNS_LPM_VaR_integer_cpp(.5, 1L, c(1, Inf, 3)))
  expect_null(NNS_LPM_VaR_integer_cpp(NA_real_, 1L, c(1, 2, 3)))
  x <- rnorm(50)
  expect_equal(UPM.VaR(.2, 2, x), LPM.VaR(.8, 2, x))
  expect_equal(LPM.VaR(c(.1, .6), 1, x), c(LPM.VaR(.1, 1, x), LPM.VaR(.6, 1, x)))
})

test_that("the compiled class-probability loop matches the R loop", {
  reference <- function(centres, shares, points, ks, dist, x.min, x.max) {
    lapply(ks, function(k) {
      t(vapply(seq_len(nrow(points)), function(i) {
        d <- .nns_mreg_distance_vector(centres, points[i, ], dist, x.min, x.max)
        exact <- which(d == 0)
        if (length(exact)) return(colMeans(shares[exact, , drop = FALSE]))
        o <- order(d)
        ds <- d[o]
        m <- sum(ds <= ds[k])
        near <- ds[seq_len(m)]
        weights <- .nns_mreg_ensemble_weights(near)
        for (v in unique(near[duplicated(near)])) {
          tied <- which(near == v)
          weights[tied] <- mean(weights[tied])
        }
        as.numeric(crossprod(weights, shares[o[seq_len(m)], , drop = FALSE]))
      }, numeric(ncol(shares))))
    })
  }
  set.seed(102)
  for (dist in c("NNS", "L2", "L1", "FACTOR")) {
    for (i in 1:6) {
      cells <- sample(c(3, 12, 40), 1); p <- sample(2:4, 1); classes <- sample(2:4, 1)
      centres <- matrix(round(rnorm(cells * p), sample(1:3, 1)), cells, p)
      shares <- matrix(runif(cells * classes), cells, classes)
      shares <- shares / rowSums(shares)
      points <- rbind(centres[sample(cells, 3), , drop = FALSE],
                      matrix(rnorm(8 * p), 8, p))
      x.min <- apply(centres, 2, min); x.max <- apply(centres, 2, max)
      ks <- unique(pmin(cells, c(1L, 2L, 5L, cells)))
      code <- match(dist, c("NNS", "L2", "L1", "FACTOR")) - 1L
      expect_equal(NNS_class_predict_cpp(centres, shares, points, ks, code, x.min, x.max),
                   reference(centres, shares, points, ks, dist, x.min, x.max),
                   tolerance = 1e-12)
    }
  }
})

test_that("cell means and column sorts match colMeans and sort", {
  set.seed(103)
  x <- matrix(rnorm(60), 20, 3)
  cells <- as.factor(sample(c("1.10", "1.2", "2.1"), 20, TRUE))
  expected <- do.call(rbind, lapply(split(seq_len(20), cells),
                                    function(g) colMeans(x[g, , drop = FALSE])))
  expect_equal(NNS_group_means_cpp(x, as.integer(cells), nlevels(cells)),
               unname(expected), tolerance = 1e-12)
  expect_identical(NNS_sort_columns_cpp(x), apply(x, 2, sort))
  x[3, 2] <- NA
  expect_null(NNS_sort_columns_cpp(x))
})

test_that("NNS.meboot's compiled rho targeting matches optimize() on ranks", {
  set.seed(104)
  n <- 40; reps <- 6
  original <- rnorm(n)
  replicates <- matrix(rnorm(n * reps), n, reps)
  r.o <- rank(original, ties.method = "average")
  aligned <- as.integer(r.o); anti <- as.integer(n + 1 - r.o)
  target <- r.o - mean(r.o)
  norm <- sqrt(sum(target^2))
  rho <- c(.9, .5, 0, -.5, -1, 1)
  solved <- NNS_meboot_rho_cpp(replicates, aligned, anti, target, norm, rho, TRUE)
  expect_false(any(solved$redo))
  for (i in seq_len(reps)) {
    sorted <- sort(replicates[, i])
    e <- sorted[aligned]; m <- sorted[anti]
    objective <- function(t) {
      ranks <- rank(t * m + (1 - t) * e, ties.method = "average")
      centred <- ranks - mean(ranks)
      abs(sum(centred * target) / (sqrt(sum(centred * centred)) * norm) - rho[i])
    }
    t.hat <- optimize(objective, interval = c(0, 1), tol = 0.01)$minimum
    expect_equal(solved$replicates[, i], t.hat * m + (1 - t.hat) * e, tolerance = 1e-12)
  }
  # A replicate with no spread is left to the R loop, where optimize() warns.
  flat <- NNS_meboot_rho_cpp(matrix(1, n, 1), aligned, anti, target, norm, .5, TRUE)
  expect_true(flat$redo)
})

test_that("dy.d_ stencils built together match stencils built one by one", {
  set.seed(105)
  for (column in list(rnorm(80), round(rnorm(80), 1), rep(c(0, 1), 40))) {
    at <- c(min(column), max(column), sample(column, 5))
    widths <- c(.05, .2)
    together <- .nns_dy_stencils(column, at, widths)
    for (b in seq_along(widths)) for (i in seq_along(at)) {
      expect_equal(together[[b]][[i]], .nns_dy_stencil(column, at[i], widths[b]))
    }
  }
  expect_error(.nns_dy_stencils(1:5, c(2, 6), .1), "inside the observed")
})

test_that("NNS.dep on a matrix equals NNS.dep on each pair", {
  set.seed(106)
  x <- cbind(a = rnorm(120), b = round(rnorm(120), 1), c = rexp(120), d = rnorm(120))
  x[, "d"] <- x[, "a"]^2 + rnorm(120, sd = .1)
  for (asym in c(FALSE, TRUE)) {
    m <- NNS.dep(x, asym = asym)
    for (i in 1:4) for (j in 1:4) if (i != j) {
      pair <- if (asym) NNS.dep(x[, i], x[, j], asym = TRUE) else
        NNS.dep(x[, min(i, j)], x[, max(i, j)])
      expect_equal(m$Dependence[i, j], pair$Dependence, tolerance = 1e-12)
      expect_equal(m$Correlation[i, j], pair$Correlation, tolerance = 1e-12)
    }
  }
})

test_that("the default thread count does not change multivariate estimates", {
  set.seed(107)
  n <- 3000
  x <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("a", "b", "c")))
  y <- x[, 1] - x[, 2]^2 + rnorm(n, sd = .3)
  one <- NNS.reg(x, y, point.est = x[1:50, ], plot = FALSE, ncores = 1)
  default <- NNS.reg(x, y, point.est = x[1:50, ], plot = FALSE)
  expect_identical(default$Point.est, one$Point.est)
  expect_identical(default$Fitted.xy, one$Fitted.xy)
  withr::with_envvar(c("_R_CHECK_LIMIT_CORES_" = "TRUE"),
                     expect_lte(.nns_default_cores(), 2L))
  expect_gte(.nns_default_cores(), 1L)
})
