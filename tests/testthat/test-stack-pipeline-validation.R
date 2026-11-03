.capture_pipeline_stack <- function(...) {
  captured <- new.env()
  suppressMessages(trace("NNS.stack", where = asNamespace("NNS"), print = FALSE,
    exit = bquote(assign("frame", environment(), envir = .(captured)))))
  on.exit(suppressMessages(untrace("NNS.stack", where = asNamespace("NNS"))))
  result <- NNS::NNS.stack(...)
  list(result = result, frame = captured$frame)
}

test_that("forecast validation uses only the supplied prefix and generated predictors", {
  set.seed(131)
  x <- data.frame(a = rnorm(45), b = rnorm(45), c = rnorm(45))
  y <- x$a - 2 * x$b + .1 * rnorm(45)
  validation <- list(IVs.train = x[1:40, ], DV.train = y[1:40],
                     IVs.test = x[41:43, ] + 1, DV.test = y[41:43])
  run <- function(xx, yy) .capture_pipeline_stack(xx, yy, IVs.test = x[41:45, ],
    .validation = validation, ts.test = 5, folds = 1, ncores = 1, status = FALSE)
  a <- run(x, y)
  changed <- x; changed[41:45, ] <- changed[41:45, ] + 100
  yy <- y; yy[41:45] <- yy[41:45] - 200
  b <- run(changed, yy)
  expect_equal(a$result$OBJfn.reg, b$result$OBJfn.reg)
  expect_equal(a$result$OBJfn.dim.red, b$result$OBJfn.dim.red)
  expect_equal(a$frame$dim_best_count, b$frame$dim_best_count)
  expect_equal(a$frame$candidate_raws, b$frame$candidate_raws)
  expect_equal(a$frame$splits[[1]]$train, 1:40)
  expect_equal(a$frame$splits[[1]]$validation, 41:43)
  expect_equal(a$frame$x.cv[41:43, ], validation$IVs.test)
  # The final fit still receives all original training rows.
  expect_equal(b$frame$x, changed)
  expect_equal(b$frame$y_fit_all, yy)
  expect_false(isTRUE(all.equal(a$result$stack, b$result$stack)))
})

test_that("stack method one retains selected coordinates instead of duplicate indices", {
  set.seed(128)
  x <- data.frame(a = rnorm(80), b = rnorm(80), c = rnorm(80))
  y <- x$a * x$b + .5 * x$c
  fit <- .capture_pipeline_stack(x, y, IVs.test = x[1:3, ], folds = 2,
    ncores = 1, status = FALSE, dim.red.method = c(a = 1, b = -1, c = .1))
  fr <- fit$frame
  selected <- which(abs(fr$dim_full_coef) > 0)
  expect_equal(fr$.method1_full_design()$train,
               fr$full_design$train[, selected, drop = FALSE])
  if (length(selected) > 1) expect_false(isTRUE(all.equal(
    fr$.method1_full_design()$train[[1]], fr$.method1_full_design()$train[[2]])))
})

test_that("pipeline validation rejects incompatible split modes", {
  x <- data.frame(a = 1:10, b = (1:10)^2)
  y <- seq(.1, 1, .1)
  v <- list(IVs.train = x[1:8, ], DV.train = y[1:8],
            IVs.test = x[9:10, ], DV.test = y[9:10])
  expect_error(NNS.stack(x, y, .validation = v, ts.test = 2, folds = 2,
    ncores = 1, status = FALSE), "folds = 1")
  v$DV.test <- NA_real_
  expect_error(NNS.stack(x, y, .validation = v, ts.test = 2, folds = 1,
    ncores = 1, status = FALSE), "finite held-out")
})

test_that("numeric two-level forecast responses remain regression with a validation fold", {
  x <- data.frame(a = seq_len(30), b = sin(seq_len(30)))
  y <- rep(c(2, 5), 15)
  v <- list(IVs.train = x[1:25, ], DV.train = y[1:25],
            IVs.test = x[26:30, ], DV.test = y[26:30])
  run <- function(response) NNS.stack(x, response, IVs.test = x[26:30, ],
    .validation = v, ts.test = 5, folds = 1, method = 1,
    ncores = 1, status = FALSE)
  a <- run(y)
  changed <- y; changed[26:30] <- 7
  b <- run(changed)
  expect_null(a$class.levels)
  expect_equal(a$OBJfn.reg, b$OBJfn.reg)
  expect_error(NNS.stack(x, y, .validation = v, ts.test = 5, folds = 1,
    type = "CLASS", ncores = 1, status = FALSE), "requires regression")
})
