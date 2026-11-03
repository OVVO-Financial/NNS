# Matched one-step evaluation; see nns_arma_conformal_benchmark_report.md.
# All outcomes used for tuning/calibration precede the scored forecast origin.
# Native NNS intervals do not define an entire predictive density.

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

ALPHA <- 0.10
TARGET_COV <- 1 - ALPHA
N_LAGS <- 12L
FIT_END <- 700L
CAL_END <- 1000L
WINDOW <- 100L
N_SEEDS <- 5L
TRAINING_FRAC <- 0.90
MAX_H <- 1L  # match the one-step information set of every comparator
NNS_NCORES <- 1L

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

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
  eps <- if (heavy_tail) rt(T, df = 5) / sqrt(5 / 3) else rnorm(T)
  y <- numeric(T)
  y[1] <- level[1] + sigma[1] * eps[1]
  for (i in 2:T) {
    y[i] <- level[i] + 0.55 * (y[i - 1] - level[i - 1]) + sigma[i] * eps[i]
  }
  data.table(t = tt, y = as.numeric(y), level = as.numeric(level), sigma = as.numeric(sigma))
}

true_conditional_mean <- function(d, raw_idx) {
  d$level[raw_idx] + 0.55 * (d$y[raw_idx - 1L] - d$level[raw_idx - 1L])
}

lag_features <- function(y, n_lags = N_LAGS) {
  n <- length(y)
  yy <- y[(n_lags + 1):n]
  X <- matrix(NA_real_, nrow = length(yy), ncol = n_lags)
  for (k in seq_len(n_lags)) {
    X[, k] <- y[(n_lags + 1 - k):(n - k)]
  }
  colnames(X) <- paste0("lag", seq_len(n_lags))
  list(X = X, yy = yy)
}

ridge_forecast <- function(X, yy, fit_end) {
  if (HAS_GLMNET) {
    n_tr <- fit_end
    lambda <- 1.0 / n_tr
    fit <- glmnet::glmnet(X[1:n_tr, , drop = FALSE], yy[1:n_tr], alpha = 0, lambda = lambda, standardize = TRUE)
    mu <- as.numeric(glmnet::predict.glmnet(fit, newx = X, s = lambda))
  } else {
    df_tr <- as.data.frame(X[1:fit_end, , drop = FALSE])
    df_tr$y <- yy[1:fit_end]
    fit <- lm(y ~ ., data = df_tr)
    mu <- as.numeric(predict(fit, newdata = as.data.frame(X)))
  }
  mu
}

coverage <- function(lo, hi, y) mean(y >= lo & y <= hi, na.rm = TRUE)
mean_width <- function(lo, hi) mean(hi - lo, na.rm = TRUE)
frac_infinite <- function(lo, hi) mean(!is.finite(lo) | !is.finite(hi))

rolling_coverage <- function(lo, hi, y, window = WINDOW) {
  n <- length(y)
  if (n < window) return(numeric(0))
  vapply(seq_len(n - window + 1L), function(i) {
    idx <- i:(i + window - 1L)
    coverage(lo[idx], hi[idx], y[idx])
  }, numeric(1))
}

worst_window_coverage <- function(lo, hi, y, window = WINDOW) {
  rc <- rolling_coverage(lo, hi, y, window)
  if (length(rc) == 0L) NA_real_ else min(rc, na.rm = TRUE)
}

interval_score <- function(lo, hi, y, alpha = ALPHA) {
  mean((hi - lo) + (2 / alpha) * pmax(lo - y, 0) + (2 / alpha) * pmax(y - hi, 0), na.rm = TRUE)
}

coverage_by_stratum <- function(lo, hi, y, sigma, k = 4L) {
  r <- rank(sigma, ties.method = "first")
  grp <- cut(r, breaks = k, labels = FALSE, include.lowest = TRUE)
  vapply(seq_len(k), function(j) {
    idx <- which(grp == j)
    if (length(idx) == 0L) NA_real_ else coverage(lo[idx], hi[idx], y[idx])
  }, numeric(1))
}

z_alpha <- function(alpha = ALPHA) qnorm(1 - alpha / 2)

gaussian_interval <- function(mu, sigma, alpha = ALPHA) {
  sigma <- pmax(as.numeric(sigma), 1e-8)
  z <- z_alpha(alpha)
  list(lo = mu - z * sigma, hi = mu + z * sigma)
}

crps_gaussian <- function(mu, sigma, y) {
  sigma <- pmax(as.numeric(sigma), 1e-8)
  z <- (y - mu) / sigma
  mean(sigma * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)), na.rm = TRUE)
}

log_score_gaussian <- function(mu, sigma, y) {
  sigma <- pmax(as.numeric(sigma), 1e-8)
  mean(-dnorm(y, mean = mu, sd = sigma, log = TRUE), na.rm = TRUE)
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

score_method <- function(method, family, lo, hi, y_te, sig_te, mu_ = NULL, s_ = NULL) {
  lo_raw <- as.numeric(lo); hi_raw <- as.numeric(hi)
  y_te <- as.numeric(y_te); sig_te <- as.numeric(sig_te)
  if (length(lo_raw) != length(hi_raw) || length(lo_raw) != length(y_te) || length(lo_raw) != length(sig_te)) {
    stop(method, ": length mismatch. lo=", length(lo_raw), ", hi=", length(hi_raw), ", y=", length(y_te), ", sigma=", length(sig_te))
  }
  lo2 <- pmin(lo_raw, hi_raw); hi2 <- pmax(lo_raw, hi_raw)
  cbs <- coverage_by_stratum(lo2, hi2, y_te, sig_te, k = 4L)
  row <- data.table(
    method = method, family = family,
    marg_cov = coverage(lo2, hi2, y_te),
    worst_win_cov = worst_window_coverage(lo2, hi2, y_te, WINDOW),
    cov_lowvol = cbs[1], cov_hivol = cbs[length(cbs)],
    cond_cov_gap = max(abs(cbs - TARGET_COV), na.rm = TRUE),
    width = mean_width(lo2, hi2), frac_inf = frac_infinite(lo2, hi2),
    interval_score = interval_score(lo2, hi2, y_te, ALPHA),
    CRPS = NA_real_, logscore = NA_real_
  )
  if (!is.null(mu_) && !is.null(s_)) {
    row$CRPS <- crps_gaussian(mu_, s_, y_te)
    row$logscore <- log_score_gaussian(mu_, s_, y_te)
  }
  row
}

fixed_split_cp <- function(mu_te, resid_cal, alpha = ALPHA) {
  scores <- sort(abs(resid_cal))
  k <- ceiling((length(scores) + 1L) * (1 - alpha))
  q <- if (k > length(scores)) Inf else scores[k]
  list(lo = mu_te - q, hi = mu_te + q)
}

aci <- function(mu_te, y_te, alpha = ALPHA, gamma = 0.03, warm = NULL) {
  n <- length(y_te); lo <- numeric(n); hi <- numeric(n)
  alpha_t <- alpha
  hist_scores <- if (!is.null(warm)) abs(warm) else numeric(0)
  for (t in seq_len(n)) {
    k <- ceiling((length(hist_scores) + 1L) * (1 - alpha_t))
    q_t <- if (length(hist_scores) == 0L || k > length(hist_scores)) Inf else sort(hist_scores)[k]
    lo[t] <- mu_te[t] - q_t; hi[t] <- mu_te[t] + q_t
    err_t <- as.integer(y_te[t] < lo[t] || y_te[t] > hi[t])
    alpha_t <- alpha_t + gamma * (alpha - err_t)
    alpha_t <- pmax(0.001, pmin(0.999, alpha_t))
    hist_scores <- c(hist_scores, abs(y_te[t] - mu_te[t]))
  }
  list(lo = lo, hi = hi)
}

agaci <- function(mu_te, y_te, alpha = ALPHA, warm = NULL, gammas = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1)) {
  experts <- lapply(gammas, function(g) aci(mu_te, y_te, alpha = alpha, gamma = g, warm = warm))
  lo <- Reduce("+", lapply(experts, `[[`, "lo")) / length(experts)
  hi <- Reduce("+", lapply(experts, `[[`, "hi")) / length(experts)
  list(lo = lo, hi = hi)
}

conformal_pid <- function(mu_te, y_te, alpha = ALPHA, warm = NULL, Kp = 0.1, Ki = 0.01, Kd = 0.001) {
  n <- length(y_te); lo <- numeric(n); hi <- numeric(n)
  hist_scores <- if (!is.null(warm)) abs(warm) else numeric(0)
  err_prev <- 0; integral <- 0
  for (t in seq_len(n)) {
    k <- ceiling((length(hist_scores) + 1L) * (1 - alpha))
    q_t <- if (length(hist_scores) == 0L || k > length(hist_scores)) Inf else sort(hist_scores)[k]
    lo[t] <- mu_te[t] - q_t; hi[t] <- mu_te[t] + q_t
    err_t <- as.integer(y_te[t] < lo[t] || y_te[t] > hi[t]) - alpha
    integral <- integral + err_t
    deriv <- err_t - err_prev
    delta <- Kp * err_t + Ki * integral + Kd * deriv
    err_prev <- err_t
    hist_scores <- c(hist_scores, abs(y_te[t] - mu_te[t]) * max(1e-6, 1 + delta))
  }
  list(lo = lo, hi = hi)
}

nexcp <- function(mu_te, y_te, alpha = ALPHA, warm = NULL, decay = 0.99) {
  n <- length(y_te); lo <- numeric(n); hi <- numeric(n)
  hist_scores <- if (!is.null(warm)) abs(warm) else numeric(0)
  hist_weights <- if (!is.null(warm)) decay ^ (rev(seq_along(warm)) - 1) else numeric(0)
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
    lo[t] <- mu_te[t] - q_t; hi[t] <- mu_te[t] + q_t
    hist_scores <- c(hist_scores, abs(y_te[t] - mu_te[t]))
    hist_weights <- c(hist_weights * decay, 1)
  }
  list(lo = lo, hi = hi)
}

recal_const <- function(mu_te, resid_cal, alpha = ALPHA) {
  s <- sd(resid_cal, na.rm = TRUE)
  z <- z_alpha(alpha)
  list(lo = mu_te - z * s, hi = mu_te + z * s, mu = mu_te, sigma = rep(s, length(mu_te)))
}

ewma_vol <- function(mu_te, y_te, alpha = ALPHA, warm = NULL, lam = 0.94) {
  all_resid <- c(if (!is.null(warm)) warm else numeric(0), y_te - mu_te)
  n_warm <- if (!is.null(warm)) length(warm) else 0L
  var_vec <- numeric(length(all_resid))
  var_vec[1] <- all_resid[1]^2
  for (i in 2:length(all_resid)) {
    var_vec[i] <- lam * var_vec[i - 1] + (1 - lam) * all_resid[i - 1]^2
  }
  sig_te <- sqrt(var_vec[(n_warm + 1):length(all_resid)])
  sig_te <- pmax(sig_te, 1e-6)
  z <- z_alpha(alpha)
  list(lo = mu_te - z * sig_te, hi = mu_te + z * sig_te, mu = mu_te, sigma = sig_te)
}

oracle_sigma_method <- function(mu_te, sig_te, alpha = ALPHA) {
  z <- z_alpha(alpha)
  list(lo = mu_te - z * sig_te, hi = mu_te + z * sig_te, mu = mu_te, sigma = sig_te)
}

garch_vol <- function(resid_tr, resid_te, mu_te, alpha = ALPHA) {
  if (!HAS_RUGARCH) return(NULL)
  tryCatch({
    spec <- rugarch::ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      distribution.model = "norm"
    )
    fit <- rugarch::ugarchfit(spec, data = resid_tr, solver = "hybrid")
    fc <- rugarch::ugarchforecast(fit, n.ahead = length(resid_te))
    sig_te <- pmax(as.numeric(rugarch::sigma(fc)), 1e-6)
    z <- z_alpha(alpha)
    list(lo = mu_te - z * sig_te, hi = mu_te + z * sig_te, mu = mu_te, sigma = sig_te)
  }, error = function(e) {
    message("  [GARCH failed: ", conditionMessage(e), "]")
    NULL
  })
}

get_nns_seas_periods <- function(training_series) {
  seas <- NNS::NNS.seas(variable = training_series, plot = FALSE)
  if (is.character(seas)) stop("NNS.seas returned character: ", paste(seas, collapse = " "))
  periods <- NULL
  if (is.list(seas)) {
    periods <- seas$Periods %||% seas$periods %||% seas$all.periods %||% seas$best.period
  }
  if (is.null(periods)) stop("NNS.seas: cannot find Periods or periods. Names: ", paste(names(seas), collapse = ", "))
  if (is.matrix(periods) || is.data.frame(periods) || data.table::is.data.table(periods)) {
    periods <- as.numeric(periods[, 1])
  }
  periods <- sort(unique(as.integer(na.omit(as.numeric(periods)))))
  periods <- periods[is.finite(periods)]
  periods <- periods[periods > 1]
  periods <- periods[periods < length(training_series)]
  if (length(periods) == 0L) stop("NNS.seas returned no usable periods.")
  periods
}

run_nns_arma_walkforward <- function(d, training_frac = TRAINING_FRAC, max_h = MAX_H) {
  if (max_h != 1L) stop("This benchmark uses matched one-step forecast origins.")
  T_raw <- nrow(d)
  current_train <- N_LAGS + CAL_END
  all_pred <- numeric(0); all_lo <- numeric(0); all_hi <- numeric(0); all_y <- numeric(0); all_sig <- numeric(0)
  chunks <- list(); chunk_id <- 0L
  while (current_train < T_raw) {
    chunk_id <- chunk_id + 1L
    remaining <- T_raw - current_train
    implied_h <- floor(current_train * (1 - training_frac) / training_frac)
    h_i <- min(max_h, remaining, max(1L, implied_h))
    end_i <- current_train + h_i
    inner_train <- floor(current_train * training_frac)
    seas_i <- get_nns_seas_periods(d$y[seq_len(inner_train)])
    message("  NNS chunk ", chunk_id, ": train=", current_train, " h=", h_i, " seas=", paste(seas_i, collapse = ","))
    fit <- NNS::NNS.ARMA.optim(
      variable = d$y[seq_len(current_train)], h = h_i, training.set = inner_train,
      seasonal.factor = seas_i, lin.only = FALSE, negative.values = TRUE,
      obj.fn = expression(mean((predicted - actual)^2)), objective = "min",
      linear.approximation = TRUE, ncores = NNS_NCORES, pred.int = TARGET_COV,
      print.trace = FALSE, plot = FALSE
    )
    pred_i <- as.numeric(fit$results)
    lo_i <- as.numeric(fit$lower.pred.int)
    hi_i <- as.numeric(fit$upper.pred.int)
    if (length(pred_i) != h_i || length(lo_i) != h_i || length(hi_i) != h_i) {
      stop("NNS.ARMA.optim length mismatch in chunk ", chunk_id)
    }
    pred_idx <- (current_train + 1L):end_i
    all_pred <- c(all_pred, pred_i)
    all_lo <- c(all_lo, pmin(lo_i, hi_i))
    all_hi <- c(all_hi, pmax(lo_i, hi_i))
    all_y <- c(all_y, d$y[pred_idx])
    all_sig <- c(all_sig, d$sigma[pred_idx])
    chunks[[chunk_id]] <- data.table(
      chunk = chunk_id, train_end = current_train, h = h_i, end = end_i,
      n_seas_periods_input = length(seas_i), seas_periods_input = paste(seas_i, collapse = ","),
      period = paste(fit$periods, collapse = ","), weights = paste(fit$weights, collapse = ","),
      method = as.character(fit$method), shrink = as.character(fit$shrink),
      nns_regress = as.character(fit$nns.regress), obj_fn = as.numeric(fit$obj.fn),
      bias_shift = as.numeric(fit$bias.shift)
    )
    current_train <- end_i
  }
  list(pred = all_pred, lo = all_lo, hi = all_hi, y = all_y, sigma = all_sig, chunks = rbindlist(chunks, fill = TRUE))
}

run_once <- function(seed = 0L, heavy_tail = FALSE) {
  d <- make_timeseries(T = 3500L, seed = seed, heavy_tail = heavy_tail)
  lf <- lag_features(d$y, N_LAGS)
  X <- lf$X; yy <- lf$yy
  mu <- ridge_forecast(X, yy, FIT_END)
  sig_all <- d$sigma[(N_LAGS + 1):nrow(d)]
  resid <- yy - mu
  te_idx <- (CAL_END + 1L):length(yy)
  mu_te <- mu[te_idx]; y_te <- yy[te_idx]; sig_te <- sig_all[te_idx]
  raw_te_idx <- te_idx + N_LAGS
  true_mu_te <- true_conditional_mean(d, raw_te_idx)
  resid_cal <- resid[(FIT_END + 1L):CAL_END]
  resid_tr <- resid[1:FIT_END]
  warm <- resid[1:CAL_END]
  methods <- list()

  methods[["fixed split (CP)"]] <- c(fixed_split_cp(mu_te, resid_cal), list(mu_ = NULL, s_ = NULL, family = "cp"))
  methods[["ACI"]] <- c(aci(mu_te, y_te, ALPHA, gamma = 0.03, warm = warm), list(mu_ = NULL, s_ = NULL, family = "cp"))
  methods[["AgACI"]] <- c(agaci(mu_te, y_te, ALPHA, warm = warm), list(mu_ = NULL, s_ = NULL, family = "cp"))
  methods[["conformal PID"]] <- c(conformal_pid(mu_te, y_te, ALPHA, warm = warm), list(mu_ = NULL, s_ = NULL, family = "cp"))
  methods[["NexCP (weighted)"]] <- c(nexcp(mu_te, y_te, ALPHA, warm = warm), list(mu_ = NULL, s_ = NULL, family = "cp"))

  oi <- gaussian_interval(true_mu_te, sig_te, ALPHA)
  methods[["oracle (true conditional mu,sigma)"]] <- list(lo = oi$lo, hi = oi$hi, mu_ = true_mu_te, s_ = sig_te, family = "oracle")
  os <- oracle_sigma_method(mu_te, sig_te)
  methods[["true sigma on est. mu"]] <- list(lo = os$lo, hi = os$hi, mu_ = mu_te, s_ = sig_te, family = "oracle")

  ew <- ewma_vol(mu_te, y_te, ALPHA, warm = warm)
  methods[["EWMA-vol Gaussian"]] <- list(lo = ew$lo, hi = ew$hi, mu_ = ew$mu, s_ = ew$sigma, family = "prob")
  rc <- recal_const(mu_te, resid_cal)
  methods[["static Gaussian (recal)"]] <- list(lo = rc$lo, hi = rc$hi, mu_ = rc$mu, s_ = rc$sigma, family = "prob")
  gv <- garch_vol(resid_tr, resid[te_idx], mu_te)
  if (!is.null(gv)) methods[["GARCH(1,1) Gaussian"]] <- list(lo = gv$lo, hi = gv$hi, mu_ = gv$mu, s_ = gv$sigma, family = "prob")

  nns_wf <- run_nns_arma_walkforward(d)
  methods[["NNS.ARMA.optim (native intervals)"]] <- list(
    lo = nns_wf$lo, hi = nns_wf$hi, mu_ = NULL, s_ = NULL,
    family = "nns"
  )

  rows <- lapply(names(methods), function(nm) {
    m <- methods[[nm]]; fam <- m$family
    if (fam == "nns") {
      lo_v <- m$lo; hi_v <- m$hi; y_v <- nns_wf$y; sig_v <- nns_wf$sigma
    } else {
      lo_v <- m$lo; hi_v <- m$hi; y_v <- y_te; sig_v <- sig_te
    }
    score_method(nm, fam, lo_v, hi_v, y_v, sig_v, m$mu_ %||% NULL, m$s_ %||% NULL)
  })

  list(scores = rbindlist(rows), methods = methods, y_te = y_te, sig_te = sig_te, nns_wf = nns_wf, mu_te = mu_te, true_mu_te = true_mu_te)
}

make_figures <- function(keep, agg) {
  methods <- keep$methods; y_te <- keep$y_te; sig_te <- keep$sig_te; nns_wf <- keep$nns_wf
  t_vec <- seq_along(y_te)
  nns_key <- "NNS.ARMA.optim (native intervals)"
  sel <- c("fixed split (CP)", "ACI", "conformal PID", nns_key, "EWMA-vol Gaussian")
  png("figures/ts_coverage.png", width = 1100, height = 550)
  plot(NULL, xlim = c(1, length(y_te) - WINDOW), ylim = c(0.4, 1.02), xlab = paste0("test step, rolling coverage window = ", WINDOW), ylab = "coverage", main = "Rolling coverage under drift")
  cols <- c("#1f4ed8", "#dc2626", "#16a34a", "#7e22ce", "#15803d")
  for (i in seq_along(sel)) {
    nm <- sel[i]
    if (!nm %in% names(methods)) next
    m <- methods[[nm]]
    if (m$family == "nns") { lo_v <- m$lo; hi_v <- m$hi; y_v <- nns_wf$y } else { lo_v <- m$lo; hi_v <- m$hi; y_v <- y_te }
    rc <- rolling_coverage(lo_v, hi_v, y_v, WINDOW)
    lines(seq_along(rc), rc, col = cols[i], lwd = 1.4)
  }
  abline(h = TARGET_COV, lty = 2, lwd = 1)
  legend("bottomleft", legend = sel, col = cols[seq_along(sel)], lwd = 1.4, cex = 0.75, ncol = 2)
  dev.off()

  png("figures/ts_plane.png", width = 950, height = 700)
  fam_cols <- c(cp = "#1f4ed8", prob = "#15803d", oracle = "#c2410c", nns = "#7e22ce")
  with(agg, {
    plot(worst_win_cov, interval_score, col = fam_cols[family], pch = 19, cex = 0.9, xlab = "worst rolling-window coverage", ylab = "interval score, lower is better", main = "Time-series efficiency versus worst-case coverage")
    abline(v = TARGET_COV, lty = 2, lwd = 1)
    text(worst_win_cov, interval_score, labels = method, cex = 0.55, pos = 3, col = fam_cols[family])
    legend("bottomleft", legend = c("conformal", "probabilistic", "oracle", "NNS"), col = unname(fam_cols), pch = 19, cex = 0.8)
  })
  dev.off()

  png("figures/ts_width.png", width = 1100, height = 500)
  z <- z_alpha(ALPHA)
  plot(t_vec, pmin(2 * z * sig_te, 30), type = "l", lwd = 1.3, xlab = "test step", ylab = "interval width", main = "Does interval width track volatility?", ylim = c(0, 30))
  if (nns_key %in% names(methods)) {
    m <- methods[[nns_key]]; w <- pmin(m$hi - m$lo, 30)
    lines(seq_along(w), w, col = "#7e22ce", lwd = 1.1)
  }
  if ("EWMA-vol Gaussian" %in% names(methods)) {
    m <- methods[["EWMA-vol Gaussian"]]
    lines(t_vec, pmin(m$hi - m$lo, 30), col = "#15803d", lwd = 1.1)
  }
  if ("fixed split (CP)" %in% names(methods)) {
    m <- methods[["fixed split (CP)"]]
    lines(t_vec, pmin(m$hi - m$lo, 30), col = "#1f4ed8", lwd = 1.1)
  }
  legend("topright", legend = c("oracle 2*z*sigma_t", nns_key, "EWMA-vol Gaussian", "fixed split (CP)"), col = c("black", "#7e22ce", "#15803d", "#1f4ed8"), lwd = c(1.3, 1.1, 1.1, 1.1), cex = 0.75)
  dev.off()
}

run_all <- function() {
  all_scores <- list(); keep <- NULL
  for (seed in 0:(N_SEEDS - 1L)) {
    message("\n=== seed ", seed, " ===")
    res <- run_once(seed = seed, heavy_tail = FALSE)
    all_scores[[length(all_scores) + 1L]] <- res$scores
    if (seed == 0L) keep <- res
  }
  scores_dt <- rbindlist(all_scores, fill = TRUE)
  metric_cols <- c("marg_cov", "worst_win_cov", "cov_lowvol", "cov_hivol", "cond_cov_gap", "width", "frac_inf", "interval_score", "CRPS", "logscore")
  agg <- scores_dt[, lapply(.SD, safe_mean), by = .(method, family), .SDcols = metric_cols][order(interval_score)]
  col_order <- c("method", "family", "marg_cov", "worst_win_cov", "cov_lowvol", "cov_hivol", "cond_cov_gap", "width", "frac_inf", "interval_score", "CRPS", "logscore")
  agg <- agg[, .SD, .SDcols = intersect(col_order, names(agg))]
  fwrite(scores_dt, "results/ts_results_all.csv")
  fwrite(agg, "results/ts_results.csv")
  agg_p <- copy(agg)
  num_cols <- names(agg_p)[sapply(agg_p, is.numeric)]
  agg_p[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]
  cat("\n=== TIME-SERIES BENCHMARK mean over ", N_SEEDS, " seeds, alpha = ", ALPHA, ", target coverage = ", TARGET_COV, " ===\n\n", sep = "")
  print(agg_p)
  cat("\nWrote:\n")
  cat("  results/ts_results.csv\n")
  cat("  results/ts_results_all.csv\n")
  make_figures(keep, agg)
  cat("  figures/ts_coverage.png\n")
  cat("  figures/ts_plane.png\n")
  cat("  figures/ts_width.png\n")
  invisible(list(scores = scores_dt, summary = agg))
}

if (sys.nframe() == 0L) results <- run_all()
