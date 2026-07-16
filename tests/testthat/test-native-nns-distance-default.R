test_that("NNS.reg NULL distance is native NNS alias", {
  x <- cbind(a = c(0, 1, 0.2, 0.9, 0.4, 0.8), b = c(0, 0.2, 1, 0.9, 0.7, 0.1))
  y <- c(0, 2, 4, 8, 5, 3)
  pts <- cbind(a = c(0.35, 0.75), b = c(0.85, 0.15))

  fit_null <- NNS.reg(x, y, point.est = pts, n.best = 3, dist = NULL,
                      plot = FALSE, residual.plot = FALSE, ncores = 1)
  fit_nns <- NNS.reg(x, y, point.est = pts, n.best = 3, dist = "nns",
                     plot = FALSE, residual.plot = FALSE, ncores = 1)

  expect_identical(fit_null$dist, "NNS")
  expect_equal(fit_null$Point.est, fit_nns$Point.est, tolerance = 1e-12)
})

test_that("native NNS distance differs from L1 and L2 when rankings differ", {
  x <- cbind(a = c(0, 0.2, 0.9, 1), b = c(0, 1, 0.2, 1))
  y <- c(0, 10, 20, 30)
  pts <- cbind(a = 0.55, b = 0.55)

  nns <- NNS.reg(x, y, point.est = pts, n.best = 2, dist = NULL,
                 plot = FALSE, residual.plot = FALSE, ncores = 1)$Point.est
  l1 <- NNS.reg(x, y, point.est = pts, n.best = 2, dist = "L1",
                plot = FALSE, residual.plot = FALSE, ncores = 1)$Point.est
  l2 <- NNS.reg(x, y, point.est = pts, n.best = 2, dist = "L2",
                plot = FALSE, residual.plot = FALSE, ncores = 1)$Point.est

  expect_false(isTRUE(all.equal(nns, l1)))
  expect_false(isTRUE(all.equal(nns, l2)))
})

test_that("native R and C++ multivariate prediction paths agree", {
  x <- cbind(a = c(0, 1, 0.2, 0.9, 0.4, 0.8), b = c(0, 0.2, 1, 0.9, 0.7, 0.1))
  y <- c(0, 2, 4, 8, 5, 3)
  pts <- cbind(a = c(0.35, 0.75), b = c(0.85, 0.15))

  fit <- NNS.reg(x, y, point.est = pts, n.best = 3, dist = NULL,
                 plot = FALSE, residual.plot = FALSE, ncores = 1)
  ref <- .nns_mreg_predict_reference(pts, fit$RPM, 3, "NNS",
                                     apply(x, 2, min), apply(x, 2, max), FALSE)
  expect_equal(fit$Point.est, ref, tolerance = 1e-12)
})

test_that("explicit distance modes and invalid values are handled", {
  x <- cbind(a = rnorm(12), b = rep(1:3, 4))
  y <- seq_len(12)
  for (d in c("NNS", "L1", "L2", "FACTOR")) {
    expect_no_error(NNS.reg(x, y, point.est = x[1:2, ], n.best = 2, dist = d,
                            plot = FALSE, residual.plot = FALSE, ncores = 1))
  }
  expect_error(NNS.reg(x, y, dist = "DTW", plot = FALSE), "dist")
})

test_that("NNS.stack NULL distance matches explicit NNS", {
  x <- data.frame(a = c(0, 1, 0.2, 0.9, 0.4, 0.8, 0.3, 0.7),
                  b = c(0, 0.2, 1, 0.9, 0.7, 0.1, 0.6, 0.4))
  y <- c(0, 2, 4, 8, 5, 3, 6, 7)
  args <- list(IVs.train = x, DV.train = y, IVs.test = x[1:2, ], method = 1,
               folds = 2, ncores = 1, status = FALSE, stack = FALSE, seed = 42)
  a <- do.call(NNS.stack, c(args, list(dist = NULL)))
  b <- do.call(NNS.stack, c(args, list(dist = "NNS")))
  expect_equal(a$reg, b$reg, tolerance = 1e-12)
})

test_that("NNS.boost NULL distance matches explicit NNS and propagates L2", {
  x <- data.frame(a = c(0, 1, 0.2, 0.9, 0.4, 0.8, 0.3, 0.7),
                  b = c(0, 0.2, 1, 0.9, 0.7, 0.1, 0.6, 0.4))
  y <- c(0, 2, 4, 8, 5, 3, 6, 7)
  args <- list(IVs.train = x, DV.train = y, IVs.test = x[1:2, ], learner.trials = 2,
               epochs = 1, status = FALSE, seed = 42)
  a <- do.call(NNS.boost, c(args, list(dist = NULL)))
  b <- do.call(NNS.boost, c(args, list(dist = "NNS")))
  expect_equal(a$results, b$results, tolerance = 1e-12)
  expect_no_error(do.call(NNS.boost, c(args, list(dist = "L2"))))
})
