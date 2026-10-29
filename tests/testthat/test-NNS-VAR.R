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

test_that("NNS.VAR relevant_variables are the selected regressors, never the target itself", {
  A <- .var_driven_data()
  for (method in c("cor", "nns.dep")) {
    out <- NNS.VAR(A, h = 3, tau = 2, dim.red.method = method, status = FALSE, ncores = 1)
    for (v in colnames(A)) {
      rel <- na.omit(out$relevant_variables[[v]])
      expect_true(length(rel) >= 1L)
      expect_false(paste0(v, "_tau_0") %in% rel)
    }
    # y is driven by lagged x
    expect_true("x_tau_1" %in% na.omit(out$relevant_variables$y))
  }
})

test_that("NNS.VAR multivariate estimates forecast the horizon, not the last observed rows", {
  L <- cbind(x1 = as.numeric(-2:17), x2 = seq(1, 39, 2))
  out <- NNS.VAR(L, h = 3, tau = 2, status = FALSE, ncores = 1)
  expect_equal(unname(as.matrix(out$univariate)), cbind(18:20, c(41, 43, 45)), tolerance = 1e-6)
  expect_equal(unname(as.matrix(out$multivariate)), cbind(18:20, c(41, 43, 45)), tolerance = 1e-6)
  expect_equal(unname(as.matrix(out$ensemble)), cbind(18:20, c(41, 43, 45)), tolerance = 1e-6)
})

test_that("NNS.VAR validates h and tau before fitting", {
  A <- .var_driven_data(40)
  expect_error(NNS.VAR(A, h = -1, status = FALSE), "\\[h\\]")
  expect_error(NNS.VAR(A, h = 1.5, status = FALSE), "\\[h\\]")
  expect_error(NNS.VAR(A, h = 2, tau = c(1, 2), status = FALSE), "\\[tau\\]")
  expect_error(NNS.VAR(A, h = 2, tau = c(1, 1, 1, 1), status = FALSE), "\\[tau\\]")
  expect_error(NNS.VAR(A, h = 2, tau = 0, status = FALSE), "\\[tau\\]")
  expect_error(NNS.VAR(A, h = 2, tau = list(1:2, 1), status = FALSE), "\\[tau\\]")
  expect_error(NNS.VAR(A[1:10, ], h = 4, tau = 2, status = FALSE), "Not enough observations")
  expect_equal(nrow(NNS.VAR(A, tau = 1, status = FALSE, ncores = 1)$ensemble), 1L)
})

test_that("NNS.VAR fills trailing NAs by extrapolation and keeps observed values", {
  L <- cbind(x1 = as.numeric(-2:17), x2 = seq(1, 39, 2))
  L[5, 1] <- NA
  L[20, 1] <- NA
  filled <- NNS.VAR(L, h = 0, tau = 2, status = FALSE, ncores = 1)
  expect_equal(filled$x1, as.numeric(-2:17), tolerance = 1e-6)
  expect_equal(filled$x2, seq(1, 39, 2))
})

test_that("NNS.VAR labels forecast rows with the following periods", {
  A <- .var_driven_data(60)
  monthly <- NNS.VAR(ts(A, start = c(2015, 1), frequency = 12), h = 2, tau = 1,
                     status = FALSE, ncores = 1)
  expect_equal(tail(rownames(monthly$interpolated_and_extrapolated), 1L), "Dec 2019")
  expect_equal(rownames(monthly$ensemble), c("Jan 2020", "Feb 2020"))

  # Weekly data used to fail with duplicate "Mon YYYY" row names.
  weekly <- NNS.VAR(ts(A, start = c(2020, 1), frequency = 52), h = 2, tau = 1,
                    status = FALSE, ncores = 1)
  expect_equal(rownames(weekly$ensemble), c("2021(9)", "2021(10)"))
})

test_that("NNS.VAR daily xts input is labelled by date", {
  skip_if_not_installed("xts")
  A <- .var_driven_data(60)
  daily <- xts::xts(A, order.by = seq(as.Date("2024-01-01"), by = "day", length.out = 60))
  out <- NNS.VAR(daily, h = 2, tau = 1, status = FALSE, ncores = 1)
  expect_equal(rownames(out$ensemble), c("2024-03-01", "2024-03-02"))
})

test_that("NNS.VAR handles one variable, awkward names and own-lag weights", {
  A <- .var_driven_data(60)
  single <- NNS.VAR(A[, 1, drop = FALSE], h = 2, tau = 2, status = FALSE, ncores = 1)
  single_rel <- na.omit(single$relevant_variables$x)
  expect_true(length(single_rel) >= 1L)
  expect_true(all(single_rel %in% c("x_tau_1", "x_tau_2")))

  B <- A; colnames(B) <- c("a_tau", "a", "a")
  out <- NNS.VAR(B, h = 2, tau = 1, naive.weights = FALSE, status = FALSE, ncores = 1)
  expect_equal(colnames(out$ensemble), c("a_tau", "a", "a.1"))
  expect_true(all(is.finite(as.matrix(out$ensemble))))
})

test_that("NNS.VAR keeps never-negative variables non-negative", {
  set.seed(3)
  P <- cbind(p1 = abs(rnorm(60)) + 0.01, p2 = seq_len(60) / 10)
  out <- NNS.VAR(P, h = 6, tau = 2, status = FALSE, ncores = 1)
  expect_true(all(as.matrix(out$multivariate) >= 0))
  expect_true(all(as.matrix(out$ensemble) >= 0))
})

test_that("NNS.VAR recovers from seasonal periods too long for a short series", {
  A <- .var_driven_data(45)
  out <- NNS.VAR(A, h = 3, tau = 12, status = FALSE, ncores = 1)
  expect_equal(dim(out$ensemble), c(3L, 3L))
  expect_true(all(is.finite(as.matrix(out$ensemble))))
})
