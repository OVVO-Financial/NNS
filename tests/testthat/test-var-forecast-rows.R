var_fixture <- function(n = 40) {
  set.seed(251)
  t <- seq_len(n)
  cbind(x = sin(t / 5) + rnorm(n, sd = .2),
        y = cos(t / 7) + rnorm(n, sd = .2))
}

# Runs NNS.VAR and records the IVs.test passed to each multivariate NNS.stack.
var_with_test_rows <- function(variables, h, tau) {
  tests <- list()
  stack <- NNS.stack
  local_mocked_bindings(
    NNS.stack = function(IVs.train, DV.train, IVs.test = NULL, ...) {
      tests[[length(tests) + 1L]] <<- IVs.test
      stack(IVs.train, DV.train, IVs.test, ...)
    },
    .package = "NNS"
  )
  res <- NNS.VAR(variables, h = h, tau = tau, ncores = 1, status = FALSE)
  list(result = res, tests = tests)
}

# Lag-matrix rows n+1..n+h: column "<var>_tau_<k>" at period t is the value at t - k.
forecast_rows <- function(res, cols, h) {
  extended <- rbind(as.data.frame(res$interpolated_and_extrapolated),
                    as.data.frame(res$univariate))
  n <- nrow(extended) - h
  vapply(cols, function(col) {
    k <- as.integer(sub(".*_tau_", "", col))
    extended[n + seq_len(h) - k, sub("_tau_[0-9]+$", "", col)]
  }, numeric(h))
}

expect_forecast_rows <- function(run, h, n_stacks) {
  expect_length(run$tests, n_stacks)
  for (iv in run$tests) {
    expected <- forecast_rows(run$result, colnames(iv), h)
    expect_equal(unname(as.matrix(iv)), unname(matrix(expected, nrow = h)))
  }
}

test_that("multi-step NNS.VAR evaluates the multivariate stage on the forecast rows", {
  run <- var_with_test_rows(var_fixture(), h = 4, tau = 1)
  expect_forecast_rows(run, h = 4, n_stacks = 2)
})

test_that("one-step NNS.VAR with per-variable lags uses the forecast row", {
  run <- var_with_test_rows(var_fixture(), h = 1, tau = c(1, 3))
  expect_forecast_rows(run, h = 1, n_stacks = 2)
})

test_that("list lags longer than h mix observed and forecast values in the forecast rows", {
  run <- var_with_test_rows(var_fixture(), h = 2, tau = list(c(1, 3), 2))
  expect_forecast_rows(run, h = 2, n_stacks = 2)
})

test_that("a single series uses its own extended lags as forecast rows", {
  run <- var_with_test_rows(var_fixture()[, "x", drop = FALSE], h = 3, tau = 2)
  expect_forecast_rows(run, h = 3, n_stacks = 1)
})

test_that("the multivariate forecast responds to the univariate forecasts", {
  shift <- 0
  local_mocked_bindings(
    NNS.ARMA.optim = function(variable, h, ...) {
      list(results = tail(variable, 1) + shift + seq_len(h) / 10)
    },
    .package = "NNS"
  )
  x <- var_fixture()
  low <- NNS.VAR(x, h = 3, tau = 1, ncores = 1, status = FALSE)
  shift <- 5
  high <- NNS.VAR(x, h = 3, tau = 1, ncores = 1, status = FALSE)
  expect_equal(as.matrix(high$univariate) - as.matrix(low$univariate),
               matrix(5, 3, 2), ignore_attr = TRUE)
  expect_gt(max(abs(as.matrix(high$multivariate) - as.matrix(low$multivariate))), 0.01)
})
