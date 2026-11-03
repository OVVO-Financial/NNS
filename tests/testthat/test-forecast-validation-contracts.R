test_that("bias is estimated before the rows used to accept it", {
  score <- expression(mean((predicted - actual)^2))
  p <- rep(10, 8)
  choose <- function(y, policy = "split", objective = "min", obj = score)
    NNS:::.nns_arma_bias_choice(p, y, obj, objective, TRUE, policy)
  persistent <- choose(rep(7, 8))
  transient <- choose(c(rep(7, 4), rep(13, 4)))
  expect_equal(persistent$estimate.rows, 1:4)
  expect_equal(persistent$acceptance.rows, 5:8)
  expect_equal(persistent$estimated.shift, -3)
  expect_equal(transient$estimated.shift, -3)
  expect_true(persistent$accepted)
  expect_equal(persistent$applied.shift, -3)
  expect_false(transient$accepted)
  expect_equal(transient$applied.shift, 0)
  expect_equal(choose(rep(7, 8), objective = "max",
                      obj = expression(-mean((predicted - actual)^2)))$applied.shift, -3)
  expect_equal(choose(rep(7, 8), "none")$status, "disabled")
})

test_that("short validation windows do not claim separate bias validation", {
  x <- NNS:::.nns_arma_bias_choice(4:6, 1:3, expression(mean((predicted - actual)^2)),
                                  "min", TRUE, "split")
  expect_equal(x$applied.shift, 0)
  expect_equal(x$status, "insufficient_validation")
  legacy <- NNS:::.nns_arma_bias_choice(4:6, 1:3, expression(mean((predicted - actual)^2)),
                                       "min", TRUE, "legacy")
  expect_equal(legacy$applied.shift, -3)
  expect_equal(legacy$estimate.rows, legacy$acceptance.rows)
})

test_that("absolute seasonality scoring is invariant to interval-scale offsets", {
  t <- seq_len(120)
  x <- 15 + sin(2 * pi * t / 12) + .2 * cos(t)
  a <- NNS.seas(x, plot = FALSE, dispersion = "absolute")
  b <- NNS.seas(x + 273.15, plot = FALSE, dispersion = "absolute")
  expect_equal(a$periods, b$periods)
  expect_equal(a$all.periods, b$all.periods, tolerance = 1e-10)
})

test_that("Monte Carlo exponent and interval validation rejects invalid grids", {
  for (e in list(0, -1, Inf, NA_real_, c(1, 2)))
    expect_error(NNS.MC(1:10, exp = e), "exp must")
  expect_error(NNS.MC(1:10, lower_rho = .8, upper_rho = .2), "lower_rho")
  expect_error(NNS.MC(1:10, by = 0), "by must")
})

test_that("residual-scale expansion can be disabled independently", {
  set.seed(5)
  x <- rnorm(60)
  run <- function(expand) {
    set.seed(8)
    NNS.meboot(x, reps = 8, rho = 0, drift = FALSE,
               expand.sd = FALSE, force.clt = FALSE,
               expand.residual.sd = expand)["replicates", ]$replicates
  }
  expect_gt(mean(apply(run(TRUE), 2, sd)), mean(apply(run(FALSE), 2, sd)))
})

test_that("VAR validation covariates cannot see observations after their origin", {
  captured <- list()
  testthat::local_mocked_bindings(
    NNS.ARMA.optim = function(variable, h, ...) {
      list(results = tail(variable, 1) + seq_len(h) / 10)
    },
    .nns_var_relevant = function(IVs.train, DV.train, ..., .validation = NULL) {
      colnames(IVs.train)
    },
    NNS.stack = function(IVs.train, DV.train, IVs.test, ..., .validation = NULL) {
      captured[[length(captured) + 1L]] <<- .validation
      list(stack = rep(mean(DV.train), nrow(IVs.test)))
    },
    .package = "NNS"
  )
  t <- seq_len(40)
  x <- cbind(a = sin(t / 5), b = cos(t / 7))
  NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
  first <- captured
  x[33:40, ] <- x[33:40, ] + 100
  captured <- list()
  NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
  expect_length(first, 2)
  for (i in 1:2) {
    expect_equal(captured[[i]]$IVs.train, first[[i]]$IVs.train)
    expect_equal(captured[[i]]$IVs.test, first[[i]]$IVs.test)
    expect_equal(captured[[i]]$DV.train, first[[i]]$DV.train)
    expect_equal(captured[[i]]$DV.test - first[[i]]$DV.test, rep(100, 8))
  }
  # Later validation lags must also use forecasts, not realized heldout values.
  expect_equal(first[[1]]$IVs.test$b_tau_0, x[32, "b"] + (1:8) / 10)
  expect_equal(tail(first[[1]]$IVs.test$b_tau_2, 6), x[32, "b"] + (1:6) / 10)
})


test_that("VAR completes a mixed-frequency prefix before seeing later releases", {
  captured <- list()
  testthat::local_mocked_bindings(
    NNS.ARMA.optim = function(variable, h, ...) {
      list(results = tail(variable, 1) + seq_len(h) / 10)
    },
    NNS.reg = function(x, y, point.est, ...) {
      list(Point.est = rep(mean(y), length(point.est)))
    },
    .nns_var_relevant = function(IVs.train, DV.train, ..., .validation = NULL) {
      colnames(IVs.train)
    },
    NNS.stack = function(IVs.train, DV.train, IVs.test, ..., .validation = NULL) {
      if (!is.null(.validation)) captured[[length(captured) + 1L]] <<- .validation
      list(stack = rep(mean(DV.train), nrow(IVs.test)))
    },
    .package = "NNS"
  )
  t <- seq_len(40)
  x <- cbind(a = 20 + sin(t / 5), b = 20 + cos(t / 7))
  x[32, "b"] <- NA_real_  # trailing missing release at the validation origin
  x[34, "a"] <- NA_real_
  x[c(35, 38), "b"] <- NA_real_
  out <- NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
  first <- captured
  expect_equal(unname(out$validation$observed.targets), c(7L, 6L))
  expect_length(first[[1]]$DV.test, 7L)
  expect_length(first[[2]]$DV.test, 6L)
  expect_equal(tail(first[[2]]$DV.train, 1), mean(x[1:31, "b"]))
  x[33:40, ] <- x[33:40, ] + 100
  captured <- list()
  NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
  for (i in 1:2) {
    expect_equal(captured[[i]]$IVs.train, first[[i]]$IVs.train)
    expect_equal(captured[[i]]$DV.train, first[[i]]$DV.train)
    expect_equal(captured[[i]]$IVs.test, first[[i]]$IVs.test)
  }
})

test_that("VAR relevance coefficients never fit heldout rows in pipeline validation", {
  trained_rows <- integer()
  testthat::local_mocked_bindings(
    NNS.reg = function(x, y, point.est = NULL, ...) {
      trained_rows <<- c(trained_rows, length(y))
      list(equation = data.frame(Coefficient = rep(1, ncol(x))),
           Point.est = rep(mean(y), NROW(point.est)))
    },
    .package = "NNS"
  )
  t <- seq_len(40)
  x <- data.frame(a = sin(t), b = cos(t), c = t)
  y <- sin(t / 5)
  fold <- list(IVs.train = x[1:32, ], DV.train = y[1:32],
               IVs.test = x[33:40, ] * 0, DV.test = y[33:40])
  NNS:::.nns_var_relevant(x, y, 8,
    expression(mean((predicted - actual)^2)), "min", "cor", .validation = fold)
  expect_true(length(trained_rows) > 0)
  expect_true(all(trained_rows == 32L))
})


test_that("numeric two-level VAR targets remain forecast regressions", {
  local_mocked_bindings(
    NNS.ARMA.optim = function(variable, h, ...) {
      list(results = rep(mean(variable), h))
    },
    .package = "NNS"
  )
  for (levels in list(c(0, 1), c(-5, 5))) {
    t <- seq_len(40)
    x <- cbind(a = rep(levels, 20), b = sin(t / 5))
    out <- NNS.VAR(x, h = 2, tau = 1, ncores = 1, status = FALSE,
                   obj.fn = expression(mean((predicted - actual)^2)))
    expect_identical(dim(out$ensemble), c(2L, 2L))
    expect_true(all(is.finite(as.matrix(out$ensemble))))
    expect_equal(unname(out$validation$observed.targets), c(8L, 8L))
  }
})
