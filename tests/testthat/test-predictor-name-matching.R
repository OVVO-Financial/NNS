test_that("test columns sharing no names with training are matched positionally", {
  # cbind() names its columns after the supplied expressions, so a test set
  # built as cbind(test.x_1, test.x_2) carries names that describe nothing
  # about the training predictors x_1 / x_2. Those names cannot express an
  # alignment, so the columns are taken in the order supplied.
  set.seed(123)
  x_1 <- rnorm(120); x_2 <- rnorm(120)
  y <- 10 * x_1 + 10 * x_2 + rnorm(120)

  set.seed(321)
  test.x_1 <- rnorm(20); test.x_2 <- rnorm(20)

  named <- NNS.stack(cbind(x_1, x_2), y, IVs.test = cbind(test.x_1, test.x_2),
                     method = 1, folds = 2, ncores = 1, status = FALSE)$stack
  bare <- NNS.stack(cbind(x_1, x_2), y,
                    IVs.test = unname(cbind(test.x_1, test.x_2)),
                    method = 1, folds = 2, ncores = 1, status = FALSE)$stack

  expect_length(named, 20L)
  expect_equal(named, bare, tolerance = 1e-12)
})

test_that("foreign test names do not silently reorder columns", {
  set.seed(4)
  x <- data.frame(a = rnorm(80), b = rnorm(80))
  y <- x$a - 2 * x$b

  ordered <- x[1:5, c("a", "b")]
  reversed <- x[1:5, c("b", "a")]

  foreign <- setNames(ordered, c("zz", "yy"))
  foreign.reversed <- setNames(reversed, c("zz", "yy"))

  base <- NNS.reg(x, y, point.est = ordered, plot = FALSE,
                  residual.plot = FALSE, ncores = 1)$Point.est

  # Positional, so the same values in the same order agree ...
  expect_equal(
    NNS.reg(x, y, point.est = foreign, plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    base, tolerance = 1e-12
  )
  # ... and swapping the columns changes the answer, since the names are not
  # consulted at all.
  expect_false(isTRUE(all.equal(
    NNS.reg(x, y, point.est = foreign.reversed, plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    base
  )))
})

test_that("matching names are still aligned by name, not by position", {
  set.seed(5)
  x <- data.frame(a = rnorm(80), b = rnorm(80))
  y <- x$a - 2 * x$b

  expect_equal(
    NNS.reg(x, y, point.est = x[1:5, c("a", "b")], plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    NNS.reg(x, y, point.est = x[1:5, c("b", "a")], plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    tolerance = 1e-12
  )
})

test_that("partially overlapping names remain an error", {
  # A test set that names some training predictors and not others is
  # ambiguous, and in practice a wrong or misspelled column.
  set.seed(6)
  x <- data.frame(a = rnorm(60), b = rnorm(60))
  y <- x$a + x$b

  bad <- setNames(x[1:5, ], c("a", "c"))

  expect_error(
    NNS.reg(x, y, point.est = bad, plot = FALSE, residual.plot = FALSE),
    "exactly match"
  )
  expect_error(
    NNS.stack(x, y, IVs.test = bad, method = 1, folds = 2, ncores = 1,
              status = FALSE),
    "exactly match"
  )
})

test_that("named test vectors follow the same rule", {
  set.seed(7)
  x <- data.frame(a = rnorm(60), b = rnorm(60))
  y <- x$a - 2 * x$b

  base <- NNS.reg(x, y, point.est = c(a = 0.5, b = -0.5), plot = FALSE,
                  residual.plot = FALSE, ncores = 1)$Point.est

  # Permutation of the training names is aligned by name.
  expect_equal(
    NNS.reg(x, y, point.est = c(b = -0.5, a = 0.5), plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    base, tolerance = 1e-12
  )
  # Names shared with nothing in training fall back to position.
  expect_equal(
    NNS.reg(x, y, point.est = c(q = 0.5, z = -0.5), plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    base, tolerance = 1e-12
  )
  expect_error(
    NNS.reg(x, y, point.est = c(a = 0.5, z = -0.5), plot = FALSE,
            residual.plot = FALSE),
    "exactly match"
  )
})

test_that("a single predictor accepts any test column name", {
  set.seed(8)
  x <- rnorm(60); y <- x^2 + rnorm(60)

  expect_equal(
    NNS.reg(x, y, point.est = cbind(test.x = x[1:4]), plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    NNS.reg(x, y, point.est = x[1:4], plot = FALSE,
            residual.plot = FALSE, ncores = 1)$Point.est,
    tolerance = 1e-12
  )
})

test_that("NNS.boost applies the same test-column rule", {
  set.seed(9)
  x <- data.frame(a = rnorm(80), b = rnorm(80), c = rnorm(80))
  y <- as.numeric(x$a + x$b > 0)

  foreign <- setNames(x[1:6, ], c("t.a", "t.b", "t.c"))
  expect_error(
    NNS.boost(IVs.train = x, DV.train = y, IVs.test = foreign,
              epochs = 8, learner.trials = 5, folds = 1, status = FALSE),
    NA
  )

  partial <- setNames(x[1:6, ], c("a", "t.b", "t.c"))
  expect_error(
    NNS.boost(IVs.train = x, DV.train = y, IVs.test = partial,
              epochs = 8, learner.trials = 5, folds = 1, status = FALSE),
    "exactly match"
  )

  by.name <- NNS.boost(IVs.train = x, DV.train = y,
                       IVs.test = x[1:6, c("c", "a", "b")],
                       epochs = 8, learner.trials = 5, folds = 1,
                       status = FALSE, seed = 7)$results
  in.order <- NNS.boost(IVs.train = x, DV.train = y,
                        IVs.test = x[1:6, c("a", "b", "c")],
                        epochs = 8, learner.trials = 5, folds = 1,
                        status = FALSE, seed = 7)$results
  expect_equal(by.name, in.order, tolerance = 1e-12)
})
