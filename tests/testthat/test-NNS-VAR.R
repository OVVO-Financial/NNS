.var_driven_data <- function(n = 80, seed = 7) {
  set.seed(seed)
  x <- cumsum(rnorm(n)); z <- rnorm(n)
  y <- c(0, 0.8 * head(x, -1)) + rnorm(n, sd = 0.3)
  cbind(x = x, y = y, z = z)
}

test_that("NNS.stack reports exactly the regressors its estimate uses", {
  A <- .var_driven_data()
  IV <- data.frame(x_tau_1 = head(A[, "x"], -1), z_tau_1 = head(A[, "z"], -1),
                   y_tau_1 = head(A[, "y"], -1))
  DV <- tail(A[, "y"], -1)
  train <- 1:70; test <- 71:74
  fit <- function(test_rows) {
    NNS.stack(IV[train, ], DV[train], IVs.test = test_rows, ts.test = 14, folds = 1,
              method = c(1, 2), order = NULL, ncores = 1, status = FALSE)
  }

  base <- fit(IV[test, ])
  selected <- base$NNS.dim.red.variables
  expect_true(length(selected) >= 1L)
  expect_true(all(selected %in% names(IV)))
  expect_named(base$NNS.dim.red.coefficients, names(IV))
  expect_setequal(names(which(base$NNS.dim.red.coefficients != 0)), selected)

  # Regressors outside the selection cannot move the estimate...
  unused <- setdiff(names(IV), selected)
  if (length(unused)) {
    shifted <- IV[test, ]
    shifted[, unused] <- shifted[, unused] + 100
    expect_equal(fit(shifted)$stack, base$stack)
  }
  # ...while a selected one does.
  shifted <- IV[test, ]
  shifted[, selected[1L]] <- shifted[, selected[1L]] + 5
  expect_false(isTRUE(all.equal(fit(shifted)$stack, base$stack)))
})


.var_two_driver_data <- function(n = 80, seed = 2) {
  set.seed(seed)
  x <- as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive"))
  w <- as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive"))
  y <- c(0, head(x, -1) + head(w, -1)) + rnorm(n, sd = 0.2)
  cbind(x = x, y = y, w = w)
}

# The lagged training matrix NNS.VAR builds for tau = 1 on complete data
.var_lag1_train <- function(A) {
  n <- nrow(A)
  lagged <- cbind(A[-1, , drop = FALSE], A[-n, , drop = FALSE])
  colnames(lagged) <- c(paste0(colnames(A), "_tau_0"), paste0(colnames(A), "_tau_1"))
  data.frame(lagged, check.names = FALSE)
}

.var_default_obj <- eval(formals(NNS.VAR)$obj.fn)

test_that("NNS.VAR relevant variables are the regressors the cross-validated threshold keeps", {
  train <- .var_lag1_train(.var_two_driver_data())
  IV <- train[, -2]; DV <- train[, 2]
  rel <- NNS:::.nns_var_relevant(IV, DV, ts.test = 16, obj.fn = .var_default_obj,
                                 objective = "min", dim.red.method = "cor")
  expect_setequal(rel, c("x_tau_1", "w_tau_1"))
  expect_true(all(rel %in% names(IV)))
})

test_that("NNS.VAR multivariate is a full NNS.stack (stack = FALSE) on the relevant variables", {
  A <- .var_two_driver_data()
  out <- NNS.VAR(A, h = 3, tau = 1, status = FALSE, ncores = 1)
  train <- .var_lag1_train(A)
  ts <- ceiling(0.2 * nrow(train))
  for (i in seq_len(ncol(A))) {
    IV <- train[, -i]
    rel <- NNS:::.nns_var_relevant(IV, train[, i], ts.test = ts, obj.fn = .var_default_obj,
                                   objective = "min", dim.red.method = "cor")
    expect_setequal(na.omit(out$relevant_variables[[i]]), rel)
    expect_false(paste0(colnames(A)[i], "_tau_0") %in% out$relevant_variables[[i]])
    fit <- suppressWarnings(NNS.stack(IV[, rel, drop = FALSE], train[, i],
                                      IVs.test = tail(IV, 3)[, rel, drop = FALSE], ts.test = ts,
                                      folds = 1, obj.fn = .var_default_obj, objective = "min",
                                      method = c(1, 2), dim.red.method = "cor", order = NULL,
                                      ncores = 1, stack = FALSE, status = FALSE))
    expect_equal(out$multivariate[[i]], as.numeric(fit$stack))
  }
})

test_that("NNS.VAR relevant variables drive the ensemble weights when naive.weights = FALSE", {
  A <- .var_two_driver_data()
  out <- NNS.VAR(A, h = 2, tau = 2, naive.weights = FALSE, status = FALSE, ncores = 1)
  for (v in colnames(A)) {
    rel <- na.omit(out$relevant_variables[[v]])
    own <- mean(sub("_tau_.*$", "", rel) == v)
    expect_equal(out$ensemble[[v]], own * out$univariate[[v]] + (1 - own) * out$multivariate[[v]])
  }
})

test_that("NNS.VAR selects relevant variables with every dimension reduction method", {
  A <- .var_two_driver_data()
  for (method in c("nns.dep", "nns.caus", "all")) {
    out <- NNS.VAR(A, h = 2, tau = 1, dim.red.method = method, status = FALSE, ncores = 1)
    rel <- na.omit(out$relevant_variables$y)
    expect_true(length(rel) >= 1L)
    expect_true(all(rel %in% c("x_tau_0", "w_tau_0", "x_tau_1", "y_tau_1", "w_tau_1")))
    expect_true(all(is.finite(out$multivariate$y)))
  }
})
