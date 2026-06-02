# NNS Time-Series Prediction Interval Benchmark

This experiment was ported from Python to R in order to facilitate the `NNS` timeseries comparison.  The original Python version is [here](https://github.com/microprediction/conformalprediction/blob/main/benchmark/run_timeseries.py).

## Outline

This benchmark compares prediction intervals for a simulated nonlinear, heteroskedastic time series. The goal is not only to check whether a method reaches the target marginal coverage of 90%, but also whether that coverage remains stable across volatility regimes and rolling windows.

The R version performs the following steps:

1. Simulates five time-series paths with a nonlinear deterministic level, autoregressive persistence, and changing volatility regimes.
2. Fits a lagged baseline model used by the conformal and probabilistic comparison methods.
3. Evaluates conformal methods: fixed split CP, ACI, AgACI, conformal PID, and NexCP.
4. Evaluates probabilistic baselines: EWMA-vol Gaussian, static Gaussian recalibration, true sigma on estimated mean, and the true conditional oracle.
5. Runs `NNS.ARMA.optim` in walk-forward chunks, estimating seasonal periods in each chunk from the available training data using `NNS.seas()`.
6. Scores all methods using marginal coverage, worst rolling-window coverage, volatility-stratified coverage, interval width, and Winkler interval score.

## Fidelity to the Python benchmark

The R script recreates the structure of the Python benchmark as closely as practical, but exact numerical equality should not be expected.

The data-generating process is matched at the equation level: deterministic level, autoregressive persistence, piecewise volatility regimes, and Gaussian innovations by default. However, Python and R use different pseudo-random number generators. Therefore, the same seed labels do not produce identical sample paths. The realized data-generating paths are similar in structure, but not numerically identical.

The baseline ridge model is also not bit-for-bit identical. The Python version uses `sklearn`'s `Ridge(alpha = 1.0)` inside a `StandardScaler` pipeline. The R version uses `glmnet` when available and falls back to `lm` otherwise. These differences may slightly alter the mean forecasts used by the conformal and probabilistic baseline methods.

The original Python benchmark also contains optional Python-specific methods such as MAPIE and `timemachines` skaters. Those are not reproduced directly in this R port. The R benchmark instead focuses on the shared conformal and probabilistic baselines plus the native `NNS` time-series comparison.

The final R version also corrects the oracle. The earlier oracle used the deterministic level as the mean. Because the DGP contains an autoregressive component, the true conditional mean is:

```text
level_t + 0.55 * (y_{t-1} - level_{t-1})
```

The results below use that corrected true conditional oracle.

## Metrics

| Metric | Meaning |
| --- | --- |
| `marg_cov` | Overall empirical coverage. Target is 0.90. |
| `worst_win_cov` | Worst rolling-window coverage using a 100-step window. Higher is better. |
| `cov_lowvol` | Coverage in the lowest-volatility stratum. |
| `cov_hivol` | Coverage in the highest-volatility stratum. |
| `cond_cov_gap` | Largest absolute deviation from 0.90 across volatility strata. Lower is better. |
| `width` | Mean interval width. Lower is sharper, conditional on adequate coverage. |
| `interval_score` | Winkler interval score. Lower is better. It penalizes both excessive width and missed observations. |
| `CRPS`, `logscore` | Proper scoring rules for methods with Gaussian predictive distributions. Not applicable to pure interval methods. |

`NA` values in `CRPS` and `logscore` are expected for interval-only methods. Conformal methods and `NNS.ARMA.optim` built-in prediction intervals produce interval bounds, not full Gaussian predictive densities.

## Results

Mean over 5 seeds, with `alpha = 0.10` and target coverage equal to `0.90`.

| Rank | Method | Family | Marginal coverage | Worst rolling coverage | Low-vol coverage | High-vol coverage | Conditional coverage gap | Width | Interval score | CRPS | Log score |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | oracle (true conditional mu,sigma) | oracle | 0.897 | 0.820 | 0.895 | 0.900 | 0.021 | 4.472 | 5.622 | 0.774 | 1.605 |
| 2 | NNS.ARMA.optim (built-in PI) | nns | 0.914 | 0.784 | 0.926 | 0.897 | 0.041 | 5.582 | 6.756 | NA | NA |
| 3 | EWMA-vol Gaussian | prob | 0.891 | 0.838 | 0.904 | 0.887 | 0.021 | 5.409 | 6.892 | 0.939 | 1.864 |
| 4 | NexCP (weighted) | cp | 0.895 | 0.782 | 0.922 | 0.891 | 0.030 | 5.467 | 6.937 | NA | NA |
| 5 | AgACI | cp | 0.905 | 0.802 | 0.948 | 0.877 | 0.048 | 5.597 | 7.035 | NA | NA |
| 6 | ACI | cp | 0.897 | 0.844 | 0.910 | 0.890 | 0.012 | 5.682 | 7.130 | NA | NA |
| 7 | true sigma on est. mu | oracle | 0.788 | 0.546 | 0.661 | 0.862 | 0.239 | 4.472 | 7.277 | 0.949 | 1.966 |
| 8 | static Gaussian (recal) | prob | 0.908 | 0.690 | 0.999 | 0.780 | 0.120 | 5.957 | 7.942 | 0.969 | 1.988 |
| 9 | fixed split (CP) | cp | 0.912 | 0.696 | 0.999 | 0.787 | 0.113 | 6.066 | 7.949 | NA | NA |
| 10 | conformal PID | cp | 0.892 | 0.572 | 1.000 | 0.746 | 0.154 | 6.091 | 8.508 | NA | NA |

## Interpretation

The corrected oracle is the expected best method. It knows both the true conditional mean and the true conditional volatility. Its interval score of 5.622 provides a useful lower bound for the other methods. Its coverage is also close to the 0.90 target, with high-volatility coverage exactly at 0.900.

Among non-oracle methods, `NNS.ARMA.optim` has the best Winkler interval score at 6.756. It also maintains strong marginal coverage at 0.914 and high-volatility coverage at 0.897. This is the central NNS result: it produces sharp intervals without sacrificing the high-volatility regime.

ACI has the strongest calibration profile among the conformal methods. Its marginal coverage is 0.897, its worst rolling-window coverage is 0.844, and its conditional coverage gap is only 0.012. This means ACI is very close to the target coverage across volatility strata. Its tradeoff is a wider average interval and a higher interval score than NNS.

NexCP is the most efficient conformal competitor in this run. It has an interval score of 6.937, marginal coverage of 0.895, and high-volatility coverage of 0.891. It is less tightly calibrated than ACI, but more efficient by the Winkler score.

The fixed split conformal baseline illustrates the global-pooling problem. Its marginal coverage is 0.912, which looks acceptable, but this average hides severe misallocation: low-volatility coverage is 0.999 while high-volatility coverage is only 0.787. In other words, the global residual pool overcovers calm regimes and undercovers volatile regimes.

The `true sigma on est. mu` row is also instructive. This method knows the true volatility but uses the estimated ridge mean. Its marginal coverage falls to 0.788. This shows that even perfect volatility information cannot rescue a structurally biased mean forecast. Prediction intervals need both a good center and a good scale.

The conformal PID method performs poorly in this run. Its worst rolling-window coverage is 0.572 and its high-volatility coverage is 0.746. The default PID controls appear poorly matched to the abrupt step changes in the simulated volatility regimes.

## Main conclusion

The benchmark supports a balanced conclusion:

> Adaptive conformal methods, especially ACI and NexCP, substantially improve on fixed split conformal in nonstationary time-series settings. However, `NNS.ARMA.optim` achieves the best non-oracle interval score while maintaining near-target marginal and high-volatility coverage.

The result is not that conformal prediction is useless. Rather, the result is that global conformal calibration can be badly misallocated in time-series data, and adaptive conformal methods must work hard to repair that. `NNS.ARMA.optim` performs strongly because the time-series model itself is adaptive before any conformal wrapper is applied.

## Appendix: Full R script

```r
# run_nns_arma_timeseries_benchmark.R
#
# Time-series benchmark in R.
#
# Methods:
#   Conformal:
#     fixed split CP, ACI, AgACI, conformal PID, NexCP weighted
#
#   Probabilistic:
#     true conditional oracle, true sigma on estimated mu,
#     EWMA-vol Gaussian, GARCH(1,1) Gaussian, static Gaussian
#
#   NNS:
#     NNS.ARMA.optim walk-forward with built-in prediction intervals
#
# Important:
#   seasonal.factor is NOT hard-coded.
#   Each NNS.ARMA.optim walk-forward chunk estimates seasonal factors from
#   the current training data using NNS.seas(training_series, plot = FALSE).
#
# Required packages:
#   NNS, data.table
#
# Optional packages:
#   glmnet   for ridge baseline
#   rugarch  for GARCH baseline
#
# Run:
#   Rscript run_nns_arma_timeseries_benchmark.R
#
# Writes:
#   results/ts_results.csv
#   results/ts_results_all.csv
#   figures/ts_coverage.png
#   figures/ts_plane.png
#   figures/ts_width.png

library(NNS)
library(data.table)

HAS_RUGARCH <- requireNamespace("rugarch", quietly = TRUE)
if (!HAS_RUGARCH) {
  message("[INFO] rugarch not available - GARCH method will be skipped.")
}

HAS_GLMNET <- requireNamespace("glmnet", quietly = TRUE)
if (!HAS_GLMNET) {
  message("[INFO] glmnet not available - ridge baseline will fall back to lm.")
}

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Constants

ALPHA <- 0.10
TARGET_COV <- 1 - ALPHA

N_LAGS <- 12L
FIT_END <- 700L
CAL_END <- 1000L
WINDOW <- 100L
N_SEEDS <- 5L

TRAINING_FRAC <- 0.90
MAX_H <- 250L

NNS_NCORES <- 1L

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Data generating process

make_timeseries <- function(T = 3500L, seed = 0L, heavy_tail = FALSE) {
  set.seed(seed + 1L)

  tt <- seq_len(T)

  level <- 0.002 * tt +
    1.50 * sin(2 * pi * tt / 50) +
    0.75 * sin(2 * pi * tt / 200)

  sigma <- rep(1.0, T)
  sigma[tt > 900 & tt <= 1400] <- 2.5
  sigma[tt > 1900 & tt <= 2450] <- 0.55
  sigma[tt > 2800] <- 1.8

  eps <- if (heavy_tail) {
    rt(T, df = 5) / sqrt(5 / 3)
  } else {
    rnorm(T)
  }

  y <- numeric(T)
  y[1] <- level[1] + sigma[1] * eps[1]

  for (i in 2:T) {
    y[i] <- level[i] + 0.55 * (y[i - 1] - level[i - 1]) + sigma[i] * eps[i]
  }

  data.table(
    t = tt,
    y = as.numeric(y),
    level = as.numeric(level),
    sigma = as.numeric(sigma)
  )
}

true_conditional_mean <- function(d, raw_idx) {
  d$level[raw_idx] + 0.55 * (d$y[raw_idx - 1L] - d$level[raw_idx - 1L])
}

# Lag-feature matrix

lag_features <- function(y, n_lags = N_LAGS) {
  n <- length(y)
  yy <- y[(n_lags + 1):n]

  X <- matrix(NA_real_, nrow = length(yy), ncol = n_lags)

  for (k in seq_len(n_lags)) {
    X[, k] <- y[(n_lags + 1 - k):(n - k)]
  }

  colnames(X) <- paste0("lag", seq_len(n_lags))

  list(
    X = X,
    yy = yy
  )
}

# Ridge baseline

ridge_forecast <- function(X, yy, fit_end) {
  if (HAS_GLMNET) {
    n_tr <- fit_end
    lambda <- 1.0 / n_tr

    fit <- glmnet::glmnet(
      X[1:n_tr, , drop = FALSE],
      yy[1:n_tr],
      alpha = 0,
      lambda = lambda,
      standardize = TRUE
    )

    mu <- as.numeric(
      glmnet::predict.glmnet(
        fit,
        newx = X,
        s = lambda
      )
    )

  } else {
    df_tr <- as.data.frame(X[1:fit_end, , drop = FALSE])
    df_tr$y <- yy[1:fit_end]

    fit <- lm(y ~ ., data = df_tr)

    mu <- as.numeric(
      predict(
        fit,
        newdata = as.data.frame(X)
      )
    )
  }

  mu
}

# Scoring helpers

coverage <- function(lo, hi, y) {
  mean(y >= lo & y <= hi, na.rm = TRUE)
}

mean_width <- function(lo, hi) {
  mean(hi - lo, na.rm = TRUE)
}

frac_infinite <- function(lo, hi) {
  mean(!is.finite(lo) | !is.finite(hi))
}

rolling_coverage <- function(lo, hi, y, window = WINDOW) {
  n <- length(y)

  if (n < window) {
    return(numeric(0))
  }

  vapply(seq_len(n - window + 1L), function(i) {
    idx <- i:(i + window - 1L)
    coverage(lo[idx], hi[idx], y[idx])
  }, numeric(1))
}

worst_window_coverage <- function(lo, hi, y, window = WINDOW) {
  rc <- rolling_coverage(lo, hi, y, window)

  if (length(rc) == 0L) {
    return(NA_real_)
  }

  min(rc, na.rm = TRUE)
}

interval_score <- function(lo, hi, y, alpha = ALPHA) {
  mean(
    (hi - lo) +
      (2 / alpha) * pmax(lo - y, 0) +
      (2 / alpha) * pmax(y - hi, 0),
    na.rm = TRUE
  )
}

coverage_by_stratum <- function(lo, hi, y, sigma, k = 4L) {
  r <- rank(sigma, ties.method = "first")

  grp <- cut(
    r,
    breaks = k,
    labels = FALSE,
    include.lowest = TRUE
  )

  vapply(seq_len(k), function(j) {
    idx <- which(grp == j)

    if (length(idx) == 0L) {
      NA_real_
    } else {
      coverage(lo[idx], hi[idx], y[idx])
    }
  }, numeric(1))
}

z_alpha <- function(alpha = ALPHA) {
  qnorm(1 - alpha / 2)
}

gaussian_interval <- function(mu, sigma, alpha = ALPHA) {
  sigma <- pmax(as.numeric(sigma), 1e-8)
  z <- z_alpha(alpha)

  list(
    lo = mu - z * sigma,
    hi = mu + z * sigma
  )
}

crps_gaussian <- function(mu, sigma, y) {
  sigma <- pmax(as.numeric(sigma), 1e-8)
  z <- (y - mu) / sigma

  mean(
    sigma * (
      z * (2 * pnorm(z) - 1) +
        2 * dnorm(z) -
        1 / sqrt(pi)
    ),
    na.rm = TRUE
  )
}

log_score_gaussian <- function(mu, sigma, y) {
  sigma <- pmax(as.numeric(sigma), 1e-8)

  mean(
    -dnorm(y, mean = mu, sd = sigma, log = TRUE),
    na.rm = TRUE
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

score_method <- function(method, family, lo, hi, y_te, sig_te,
                         mu_ = NULL, s_ = NULL) {
  lo_raw <- as.numeric(lo)
  hi_raw <- as.numeric(hi)
  y_te <- as.numeric(y_te)
  sig_te <- as.numeric(sig_te)

  if (length(lo_raw) != length(hi_raw) ||
      length(lo_raw) != length(y_te) ||
      length(lo_raw) != length(sig_te)) {
    stop(
      method,
      ": length mismatch. lo=", length(lo_raw),
      ", hi=", length(hi_raw),
      ", y=", length(y_te),
      ", sigma=", length(sig_te)
    )
  }

  lo2 <- pmin(lo_raw, hi_raw)
  hi2 <- pmax(lo_raw, hi_raw)

  cbs <- coverage_by_stratum(lo2, hi2, y_te, sig_te, k = 4L)

  row <- data.table(
    method = method,
    family = family,
    marg_cov = coverage(lo2, hi2, y_te),
    worst_win_cov = worst_window_coverage(lo2, hi2, y_te, WINDOW),
    cov_lowvol = cbs[1],
    cov_hivol = cbs[length(cbs)],
    cond_cov_gap = max(abs(cbs - TARGET_COV), na.rm = TRUE),
    width = mean_width(lo2, hi2),
    frac_inf = frac_infinite(lo2, hi2),
    interval_score = interval_score(lo2, hi2, y_te, ALPHA),
    CRPS = NA_real_,
    logscore = NA_real_
  )

  if (!is.null(mu_) && !is.null(s_)) {
    row$CRPS <- crps_gaussian(mu_, s_, y_te)
    row$logscore <- log_score_gaussian(mu_, s_, y_te)
  }

  row
}

# Conformal methods

fixed_split_cp <- function(mu_te, resid_cal, alpha = ALPHA) {
  scores <- sort(abs(resid_cal))
  k <- ceiling((length(scores) + 1L) * (1 - alpha))
  q <- if (k > length(scores)) Inf else scores[k]

  list(
    lo = mu_te - q,
    hi = mu_te + q
  )
}

aci <- function(mu_te, y_te, alpha = ALPHA, gamma = 0.03, warm = NULL) {
  n <- length(y_te)

  lo <- numeric(n)
  hi <- numeric(n)

  alpha_t <- alpha
  hist_scores <- if (!is.null(warm)) abs(warm) else numeric(0)

  for (t in seq_len(n)) {
    k <- ceiling((length(hist_scores) + 1L) * (1 - alpha_t))

    q_t <- if (length(hist_scores) == 0L || k > length(hist_scores)) {
      Inf
    } else {
      sort(hist_scores)[k]
    }

    lo[t] <- mu_te[t] - q_t
    hi[t] <- mu_te[t] + q_t

    err_t <- as.integer(y_te[t] < lo[t] || y_te[t] > hi[t])

    alpha_t <- alpha_t + gamma * (alpha - err_t)
    alpha_t <- pmax(0.001, pmin(0.999, alpha_t))

    hist_scores <- c(hist_scores, abs(y_te[t] - mu_te[t]))
  }

  list(
    lo = lo,
    hi = hi
  )
}

agaci <- function(mu_te, y_te, alpha = ALPHA, warm = NULL,
                  gammas = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1)) {
  experts <- lapply(gammas, function(g) {
    aci(mu_te, y_te, alpha = alpha, gamma = g, warm = warm)
  })

  lo <- Reduce("+", lapply(experts, `[[`, "lo")) / length(experts)
  hi <- Reduce("+", lapply(experts, `[[`, "hi")) / length(experts)

  list(
    lo = lo,
    hi = hi
  )
}

conformal_pid <- function(mu_te, y_te, alpha = ALPHA, warm = NULL,
                          Kp = 0.1, Ki = 0.01, Kd = 0.001) {
  n <- length(y_te)

  lo <- numeric(n)
  hi <- numeric(n)

  hist_scores <- if (!is.null(warm)) abs(warm) else numeric(0)

  err_prev <- 0
  integral <- 0

  for (t in seq_len(n)) {
    k <- ceiling((length(hist_scores) + 1L) * (1 - alpha))

    q_t <- if (length(hist_scores) == 0L || k > length(hist_scores)) {
      Inf
    } else {
      sort(hist_scores)[k]
    }

    lo[t] <- mu_te[t] - q_t
    hi[t] <- mu_te[t] + q_t

    err_t <- as.integer(y_te[t] < lo[t] || y_te[t] > hi[t]) - alpha

    integral <- integral + err_t
    deriv <- err_t - err_prev
    delta <- Kp * err_t + Ki * integral + Kd * deriv

    err_prev <- err_t

    hist_scores <- c(
      hist_scores,
      abs(y_te[t] - mu_te[t]) * max(1e-6, 1 + delta)
    )
  }

  list(
    lo = lo,
    hi = hi
  )
}

nexcp <- function(mu_te, y_te, alpha = ALPHA, warm = NULL, decay = 0.99) {
  n <- length(y_te)

  lo <- numeric(n)
  hi <- numeric(n)

  hist_scores <- if (!is.null(warm)) abs(warm) else numeric(0)
  hist_weights <- if (!is.null(warm)) {
    decay ^ (rev(seq_along(warm)) - 1)
  } else {
    numeric(0)
  }

  for (t in seq_len(n)) {
    if (length(hist_scores) == 0L) {
      q_t <- Inf
    } else {
      w_norm <- hist_weights / sum(hist_weights)

      ord <- order(hist_scores)
      cum_w <- cumsum(w_norm[ord])

      idx <- which(cum_w >= (1 - alpha))[1]
      q_t <- if (is.na(idx)) Inf else hist_scores[ord[idx]]
    }

    lo[t] <- mu_te[t] - q_t
    hi[t] <- mu_te[t] + q_t

    hist_scores <- c(hist_scores, abs(y_te[t] - mu_te[t]))
    hist_weights <- c(hist_weights * decay, 1)
  }

  list(
    lo = lo,
    hi = hi
  )
}

# Probabilistic baselines

recal_const <- function(mu_te, resid_cal, alpha = ALPHA) {
  s <- sd(resid_cal, na.rm = TRUE)
  z <- z_alpha(alpha)

  list(
    lo = mu_te - z * s,
    hi = mu_te + z * s,
    mu = mu_te,
    sigma = rep(s, length(mu_te))
  )
}

ewma_vol <- function(mu_te, y_te, alpha = ALPHA, warm = NULL, lam = 0.94) {
  all_resid <- c(
    if (!is.null(warm)) warm else numeric(0),
    y_te - mu_te
  )

  n_warm <- if (!is.null(warm)) length(warm) else 0L

  var_vec <- numeric(length(all_resid))
  var_vec[1] <- all_resid[1]^2

  for (i in 2:length(all_resid)) {
    var_vec[i] <- lam * var_vec[i - 1] + (1 - lam) * all_resid[i - 1]^2
  }

  sig_te <- sqrt(var_vec[(n_warm + 1):length(all_resid)])
  sig_te <- pmax(sig_te, 1e-6)

  z <- z_alpha(alpha)

  list(
    lo = mu_te - z * sig_te,
    hi = mu_te + z * sig_te,
    mu = mu_te,
    sigma = sig_te
  )
}

oracle_sigma_method <- function(mu_te, sig_te, alpha = ALPHA) {
  z <- z_alpha(alpha)

  list(
    lo = mu_te - z * sig_te,
    hi = mu_te + z * sig_te,
    mu = mu_te,
    sigma = sig_te
  )
}

garch_vol <- function(resid_tr, resid_te, mu_te, alpha = ALPHA) {
  if (!HAS_RUGARCH) {
    return(NULL)
  }

  tryCatch({
    spec <- rugarch::ugarchspec(
      variance.model = list(
        model = "sGARCH",
        garchOrder = c(1, 1)
      ),
      mean.model = list(
        armaOrder = c(0, 0),
        include.mean = FALSE
      ),
      distribution.model = "norm"
    )

    fit <- rugarch::ugarchfit(
      spec,
      data = resid_tr,
      solver = "hybrid"
    )

    fc <- rugarch::ugarchforecast(
      fit,
      n.ahead = length(resid_te)
    )

    sig_te <- as.numeric(rugarch::sigma(fc))
    sig_te <- pmax(sig_te, 1e-6)

    z <- z_alpha(alpha)

    list(
      lo = mu_te - z * sig_te,
      hi = mu_te + z * sig_te,
      mu = mu_te,
      sigma = sig_te
    )

  }, error = function(e) {
    message("  [GARCH failed: ", conditionMessage(e), "]")
    NULL
  })
}

# NNS seasonal period detection

get_nns_seas_periods <- function(training_series) {
  seas <- NNS::NNS.seas(
    variable = training_series,
    plot = FALSE
  )

  if (is.character(seas)) {
    stop(
      "NNS.seas returned character: ",
      paste(seas, collapse = " ")
    )
  }

  periods <- NULL

  if (is.list(seas)) {
    periods <- seas$Periods %||%
      seas$periods %||%
      seas$all.periods %||%
      seas$best.period
  }

  if (is.null(periods)) {
    stop(
      "NNS.seas: cannot find Periods or periods. Names: ",
      paste(names(seas), collapse = ", ")
    )
  }

  if (is.matrix(periods) ||
      is.data.frame(periods) ||
      data.table::is.data.table(periods)) {
    periods <- as.numeric(periods[, 1])
  }

  periods <- sort(unique(as.integer(na.omit(as.numeric(periods)))))
  periods <- periods[is.finite(periods)]
  periods <- periods[periods > 1]
  periods <- periods[periods < length(training_series)]

  if (length(periods) == 0L) {
    stop("NNS.seas returned no usable periods.")
  }

  periods
}

# NNS.ARMA.optim walk-forward

run_nns_arma_walkforward <- function(d,
                                     training_frac = TRAINING_FRAC,
                                     max_h = MAX_H) {
  T_raw <- nrow(d)

  current_train <- N_LAGS + CAL_END

  all_pred <- numeric(0)
  all_lo <- numeric(0)
  all_hi <- numeric(0)
  all_y <- numeric(0)
  all_sig <- numeric(0)

  chunks <- list()
  chunk_id <- 0L

  while (current_train < T_raw) {
    chunk_id <- chunk_id + 1L

    remaining <- T_raw - current_train

    implied_h <- floor(
      current_train * (1 - training_frac) / training_frac
    )

    h_i <- min(
      max_h,
      remaining,
      max(1L, implied_h)
    )

    end_i <- current_train + h_i

    training_series_i <- d$y[1:current_train]

    seas_i <- get_nns_seas_periods(training_series_i)

    message(
      "  NNS chunk ", chunk_id,
      ": train=", current_train,
      " h=", h_i,
      " seas=", paste(seas_i, collapse = ",")
    )

    fit <- NNS::NNS.ARMA.optim(
      variable = d$y[1:end_i],
      h = NULL,
      training.set = current_train,
      seasonal.factor = seas_i,
      lin.only = FALSE,
      negative.values = TRUE,
      obj.fn = expression(mean((predicted - actual)^2)),
      objective = "min",
      linear.approximation = TRUE,
      ncores = NNS_NCORES,
      pred.int = TARGET_COV,
      print.trace = FALSE,
      plot = FALSE
    )

    pred_i <- as.numeric(fit$results)
    lo_i <- as.numeric(fit$lower.pred.int)
    hi_i <- as.numeric(fit$upper.pred.int)

    if (length(pred_i) != h_i ||
        length(lo_i) != h_i ||
        length(hi_i) != h_i) {
      stop(
        "NNS.ARMA.optim length mismatch in chunk ", chunk_id,
        ". Expected h=", h_i,
        ", got results=", length(pred_i),
        ", lower=", length(lo_i),
        ", upper=", length(hi_i)
      )
    }

    pred_idx <- (current_train + 1L):end_i

    all_pred <- c(all_pred, pred_i)
    all_lo <- c(all_lo, pmin(lo_i, hi_i))
    all_hi <- c(all_hi, pmax(lo_i, hi_i))
    all_y <- c(all_y, d$y[pred_idx])
    all_sig <- c(all_sig, d$sigma[pred_idx])

    chunks[[chunk_id]] <- data.table(
      chunk = chunk_id,
      train_end = current_train,
      h = h_i,
      end = end_i,
      n_seas_periods_input = length(seas_i),
      seas_periods_input = paste(seas_i, collapse = ","),
      period = paste(fit$period, collapse = ","),
      weights = paste(fit$weights, collapse = ","),
      method = as.character(fit$method),
      shrink = as.character(fit$shrink),
      nns_regress = as.character(fit$nns.regress),
      obj_fn = as.numeric(fit$obj.fn),
      bias_shift = as.numeric(fit$bias.shift)
    )

    current_train <- end_i
  }

  list(
    pred = all_pred,
    lo = all_lo,
    hi = all_hi,
    y = all_y,
    sigma = all_sig,
    chunks = rbindlist(chunks, fill = TRUE)
  )
}

# Main per-seed run

run_once <- function(seed = 0L, heavy_tail = FALSE) {
  d <- make_timeseries(
    T = 3500L,
    seed = seed,
    heavy_tail = heavy_tail
  )

  lf <- lag_features(d$y, N_LAGS)

  X <- lf$X
  yy <- lf$yy

  mu <- ridge_forecast(X, yy, FIT_END)

  sig_all <- d$sigma[(N_LAGS + 1):nrow(d)]

  resid <- yy - mu

  te_idx <- (CAL_END + 1L):length(yy)

  mu_te <- mu[te_idx]
  y_te <- yy[te_idx]
  sig_te <- sig_all[te_idx]

  raw_te_idx <- te_idx + N_LAGS
  true_mu_te <- true_conditional_mean(d, raw_te_idx)

  resid_cal <- resid[(FIT_END + 1L):CAL_END]
  resid_tr <- resid[1:FIT_END]
  warm <- resid[1:CAL_END]

  methods <- list()

  # Conformal baselines

  fs <- fixed_split_cp(mu_te, resid_cal)
  methods[["fixed split (CP)"]] <- c(
    fs,
    list(mu_ = NULL, s_ = NULL, family = "cp")
  )

  ac <- aci(mu_te, y_te, ALPHA, gamma = 0.03, warm = warm)
  methods[["ACI"]] <- c(
    ac,
    list(mu_ = NULL, s_ = NULL, family = "cp")
  )

  ag <- agaci(mu_te, y_te, ALPHA, warm = warm)
  methods[["AgACI"]] <- c(
    ag,
    list(mu_ = NULL, s_ = NULL, family = "cp")
  )

  pid <- conformal_pid(mu_te, y_te, ALPHA, warm = warm)
  methods[["conformal PID"]] <- c(
    pid,
    list(mu_ = NULL, s_ = NULL, family = "cp")
  )

  nx <- nexcp(mu_te, y_te, ALPHA, warm = warm)
  methods[["NexCP (weighted)"]] <- c(
    nx,
    list(mu_ = NULL, s_ = NULL, family = "cp")
  )

  # Oracle baselines

  oi <- gaussian_interval(true_mu_te, sig_te, ALPHA)
  methods[["oracle (true conditional mu,sigma)"]] <- list(
    lo = oi$lo,
    hi = oi$hi,
    mu_ = true_mu_te,
    s_ = sig_te,
    family = "oracle"
  )

  os <- oracle_sigma_method(mu_te, sig_te)
  methods[["true sigma on est. mu"]] <- list(
    lo = os$lo,
    hi = os$hi,
    mu_ = mu_te,
    s_ = sig_te,
    family = "oracle"
  )

  # Probabilistic baselines

  ew <- ewma_vol(mu_te, y_te, ALPHA, warm = warm)
  methods[["EWMA-vol Gaussian"]] <- list(
    lo = ew$lo,
    hi = ew$hi,
    mu_ = ew$mu,
    s_ = ew$sigma,
    family = "prob"
  )

  rc <- recal_const(mu_te, resid_cal)
  methods[["static Gaussian (recal)"]] <- list(
    lo = rc$lo,
    hi = rc$hi,
    mu_ = rc$mu,
    s_ = rc$sigma,
    family = "prob"
  )

  gv <- garch_vol(resid_tr, resid[te_idx], mu_te)

  if (!is.null(gv)) {
    methods[["GARCH(1,1) Gaussian"]] <- list(
      lo = gv$lo,
      hi = gv$hi,
      mu_ = gv$mu,
      s_ = gv$sigma,
      family = "prob"
    )
  }

  # NNS walk-forward

  nns_wf <- run_nns_arma_walkforward(d)

  methods[["NNS.ARMA.optim (built-in PI)"]] <- list(
    lo = nns_wf$lo,
    hi = nns_wf$hi,
    mu_ = NULL,
    s_ = NULL,
    family = "nns"
  )

  rows <- lapply(names(methods), function(nm) {
    m <- methods[[nm]]
    fam <- m$family

    if (fam == "nns") {
      lo_v <- m$lo
      hi_v <- m$hi
      y_v <- nns_wf$y
      sig_v <- nns_wf$sigma
      mu_v <- m$mu_
      s_v <- m$s_
    } else {
      lo_v <- m$lo
      hi_v <- m$hi
      y_v <- y_te
      sig_v <- sig_te
      mu_v <- m$mu_
      s_v <- m$s_
    }

    score_method(
      nm,
      fam,
      lo_v,
      hi_v,
      y_v,
      sig_v,
      mu_v,
      s_v
    )
  })

  list(
    scores = rbindlist(rows),
    methods = methods,
    y_te = y_te,
    sig_te = sig_te,
    nns_wf = nns_wf,
    mu_te = mu_te,
    true_mu_te = true_mu_te
  )
}

# Figures

make_figures <- function(keep, agg) {
  methods <- keep$methods
  y_te <- keep$y_te
  sig_te <- keep$sig_te
  nns_wf <- keep$nns_wf

  t_vec <- seq_along(y_te)

  # Rolling coverage over time.
  sel <- c(
    "fixed split (CP)",
    "ACI",
    "conformal PID",
    "NNS.ARMA.optim (built-in PI)",
    "EWMA-vol Gaussian"
  )

  png("figures/ts_coverage.png", width = 1100, height = 550)

  plot(
    NULL,
    xlim = c(1, length(y_te) - WINDOW),
    ylim = c(0.4, 1.02),
    xlab = paste0("test step, rolling coverage window = ", WINDOW),
    ylab = "coverage",
    main = "Rolling coverage under drift"
  )

  cols <- c(
    "#1f4ed8",
    "#dc2626",
    "#16a34a",
    "#7e22ce",
    "#15803d"
  )

  for (i in seq_along(sel)) {
    nm <- sel[i]

    if (!nm %in% names(methods)) {
      next
    }

    m <- methods[[nm]]

    if (m$family == "nns") {
      lo_v <- m$lo
      hi_v <- m$hi
      y_v <- nns_wf$y
    } else {
      lo_v <- m$lo
      hi_v <- m$hi
      y_v <- y_te
    }

    rc <- rolling_coverage(lo_v, hi_v, y_v, WINDOW)

    lines(
      seq_along(rc),
      rc,
      col = cols[i],
      lwd = 1.4
    )
  }

  abline(h = TARGET_COV, lty = 2, lwd = 1)

  legend(
    "bottomleft",
    legend = sel,
    col = cols[seq_along(sel)],
    lwd = 1.4,
    cex = 0.75,
    ncol = 2
  )

  dev.off()

  # Coverage-sharpness plane.
  png("figures/ts_plane.png", width = 950, height = 700)

  fam_cols <- c(
    cp = "#1f4ed8",
    prob = "#15803d",
    oracle = "#c2410c",
    nns = "#7e22ce"
  )

  with(agg, {
    plot(
      worst_win_cov,
      interval_score,
      col = fam_cols[family],
      pch = 19,
      cex = 0.9,
      xlab = "worst rolling-window coverage",
      ylab = "interval score, lower is better",
      main = "Time-series efficiency versus worst-case coverage"
    )

    abline(v = TARGET_COV, lty = 2, lwd = 1)

    text(
      worst_win_cov,
      interval_score,
      labels = method,
      cex = 0.55,
      pos = 3,
      col = fam_cols[family]
    )

    legend(
      "bottomleft",
      legend = c("conformal", "probabilistic", "oracle", "NNS"),
      col = unname(fam_cols),
      pch = 19,
      cex = 0.8
    )
  })

  dev.off()

  # Interval width over time.
  png("figures/ts_width.png", width = 1100, height = 500)

  z <- z_alpha(ALPHA)

  plot(
    t_vec,
    pmin(2 * z * sig_te, 30),
    type = "l",
    lwd = 1.3,
    xlab = "test step",
    ylab = "interval width",
    main = "Does interval width track volatility?",
    ylim = c(0, 30)
  )

  nns_key <- "NNS.ARMA.optim (built-in PI)"

  if (nns_key %in% names(methods)) {
    m <- methods[[nns_key]]
    w <- pmin(m$hi - m$lo, 30)

    lines(
      seq_along(w),
      w,
      col = "#7e22ce",
      lwd = 1.1
    )
  }

  ew_key <- "EWMA-vol Gaussian"

  if (ew_key %in% names(methods)) {
    m <- methods[[ew_key]]

    lines(
      t_vec,
      pmin(m$hi - m$lo, 30),
      col = "#15803d",
      lwd = 1.1
    )
  }

  fs_key <- "fixed split (CP)"

  if (fs_key %in% names(methods)) {
    m <- methods[[fs_key]]

    lines(
      t_vec,
      pmin(m$hi - m$lo, 30),
      col = "#1f4ed8",
      lwd = 1.1
    )
  }

  legend(
    "topright",
    legend = c(
      "oracle 2*z*sigma_t",
      nns_key,
      ew_key,
      fs_key
    ),
    col = c(
      "black",
      "#7e22ce",
      "#15803d",
      "#1f4ed8"
    ),
    lwd = c(1.3, 1.1, 1.1, 1.1),
    cex = 0.75
  )

  dev.off()
}

# Aggregate over seeds

run_all <- function() {
  all_scores <- list()
  keep <- NULL

  for (seed in 0:(N_SEEDS - 1L)) {
    message("\n=== seed ", seed, " ===")

    res <- run_once(
      seed = seed,
      heavy_tail = FALSE
    )

    all_scores[[length(all_scores) + 1L]] <- res$scores

    if (seed == 0L) {
      keep <- res
    }
  }

  scores_dt <- rbindlist(all_scores, fill = TRUE)

  metric_cols <- c(
    "marg_cov",
    "worst_win_cov",
    "cov_lowvol",
    "cov_hivol",
    "cond_cov_gap",
    "width",
    "frac_inf",
    "interval_score",
    "CRPS",
    "logscore"
  )

  agg <- scores_dt[
    ,
    lapply(.SD, safe_mean),
    by = .(method, family),
    .SDcols = metric_cols
  ][order(interval_score)]

  col_order <- c(
    "method",
    "family",
    "marg_cov",
    "worst_win_cov",
    "cov_lowvol",
    "cov_hivol",
    "cond_cov_gap",
    "width",
    "frac_inf",
    "interval_score",
    "CRPS",
    "logscore"
  )

  agg <- agg[, .SD, .SDcols = intersect(col_order, names(agg))]

  fwrite(scores_dt, "results/ts_results_all.csv")
  fwrite(agg, "results/ts_results.csv")

  agg_p <- copy(agg)
  num_cols <- names(agg_p)[sapply(agg_p, is.numeric)]
  agg_p[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]

  cat(
    "\n=== TIME-SERIES BENCHMARK",
    " mean over ", N_SEEDS,
    " seeds, alpha = ", ALPHA,
    ", target coverage = ", TARGET_COV,
    " ===\n\n",
    sep = ""
  )

  print(agg_p)

  cat("\nWrote:\n")
  cat("  results/ts_results.csv\n")
  cat("  results/ts_results_all.csv\n")

  make_figures(keep, agg)

  cat("  figures/ts_coverage.png\n")
  cat("  figures/ts_plane.png\n")
  cat("  figures/ts_width.png\n")

  invisible(
    list(
      scores = scores_dt,
      summary = agg
    )
  )
}

# Entry point

results <- run_all()
```
