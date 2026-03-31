# =============================================================================
# NNS.diff BLACK-BOX / NON-ANALYTIC BENCHMARK
# Focus: projected derivative (DERIVATIVE) vs real-only baselines
# =============================================================================

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

safe_rel_err <- function(est, truth) {
  if (length(est) == 0 || is.na(est) || is.nan(est) || is.infinite(est)) {
    return(NA_real_)
  }
  if (is.na(truth) || is.nan(truth) || is.infinite(truth)) {
    return(NA_real_)
  }
  if (abs(truth) < 1e-12) {
    return(abs(est - truth))
  }
  abs((est - truth) / truth)
}

safe_abs_err <- function(est, truth) {
  if (length(est) == 0 || is.na(est) || is.nan(est) || is.infinite(est)) {
    return(NA_real_)
  }
  if (is.na(truth) || is.nan(truth) || is.infinite(truth)) {
    return(NA_real_)
  }
  abs(est - truth)
}

safe_sign <- function(x, tol = 1e-12) {
  if (is.na(x) || is.nan(x) || is.infinite(x)) return(NA_integer_)
  if (abs(x) <= tol) return(0L)
  if (x > 0) return(1L)
  -1L
}

tangent_pred_error <- function(f, x0, d_hat, delta) {
  fx0 <- tryCatch(f(x0), error = function(e) NA_real_)
  fx1 <- tryCatch(f(x0 + delta), error = function(e) NA_real_)
  if (any(!is.finite(c(fx0, fx1, d_hat)))) return(NA_real_)
  abs(fx1 - (fx0 + d_hat * delta))
}

centered_fd <- function(f, x0, h) {
  tryCatch(
    (f(x0 + h) - f(x0 - h)) / (2 * h),
    error = function(e) NA_real_
  )
}

richardson_diff <- function(f, point, h = 0.1, order = 6) {
  Rmat <- matrix(NA_real_, nrow = order, ncol = order)
  
  for (i in seq_len(order)) {
    hi <- h / (2^(i - 1))
    Rmat[i, 1] <- centered_fd(f, point, hi)
  }
  
  for (j in 2:order) {
    for (i in j:order) {
      Rmat[i, j] <- Rmat[i, j - 1] +
        (Rmat[i, j - 1] - Rmat[i - 1, j - 1]) / (4^(j - 1) - 1)
    }
  }
  
  Rmat[order, order]
}

make_frozen_noisy_function <- function(fun, sigma, rep_seed, N_evals = 1) {
  force(fun)
  force(sigma)
  force(rep_seed)
  force(N_evals)
  
  function(x) {
    # Real-only benchmark by design
    if (is.complex(x)) {
      stop("Complex input intentionally unsupported in this benchmark.")
    }
    
    x_vec <- as.numeric(x)
    
    vals <- vapply(x_vec, function(xx) {
      local_seed <- rep_seed + as.integer(abs(round(xx * 1e8))) %% 100000000L
      
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        had_seed <- TRUE
      } else {
        had_seed <- FALSE
      }
      
      on.exit({
        if (had_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
      
      set.seed(local_seed)
      mean(replicate(N_evals, fun(xx) + rnorm(1, mean = 0, sd = sigma)))
    }, numeric(1))
    
    if (length(vals) == 1) vals[[1]] else vals
  }
}

# -----------------------------------------------------------------------------
# NNS.diff wrappers
# -----------------------------------------------------------------------------

extract_nns_rows <- function(out) {
  blank <- list(
    derivative = NA_real_,
    inferred_h = NA_real_,
    iterations = NA_real_,
    converged = NA_real_,
    termination.code = NA_real_,
    fin_init_avg = NA_real_,
    fin_inf_avg = NA_real_
  )
  
  if (is.null(out) || !(is.matrix(out) || is.data.frame(out))) return(blank)
  
  rn <- rownames(out)
  if (is.null(rn)) return(blank)
  
  vals <- suppressWarnings(as.numeric(out[, 1]))
  names(vals) <- rn
  
  get1 <- function(name) {
    if (name %in% names(vals)) as.numeric(vals[name][1]) else NA_real_
  }
  
  # supports both clean names and the duplicated-name version
  fin_init_name <- if ("Initial h averaged finite step" %in% names(vals)) {
    "Initial h averaged finite step"
  } else {
    "Initial h averaged finite step.Averaged Finite Step"
  }
  
  fin_inf_name <- if ("Inferred h averaged finite step" %in% names(vals)) {
    "Inferred h averaged finite step"
  } else {
    "Inferred h averaged finite step.Averaged Finite Step"
  }
  
  list(
    derivative = get1("DERIVATIVE"),
    inferred_h = get1("Inferred h"),
    iterations = get1("iterations"),
    converged = get1("converged"),
    termination.code = get1("termination.code"),
    fin_init_avg = get1(fin_init_name),
    fin_inf_avg = get1(fin_inf_name)
  )
}

safe_nns_diff <- function(f, point, h = 0.1, tol = 1e-10, max.iter = 1000L,
                          digits = 12, print.trace = FALSE) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    if (file.exists(tf)) unlink(tf)
  }, add = TRUE)
  
  tryCatch(
    NNS.diff(
      f = f,
      point = point,
      h = h,
      tol = tol,
      max.iter = max.iter,
      digits = digits,
      print.trace = print.trace
    ),
    error = function(e) structure(
      list(error_message = conditionMessage(e)),
      class = "nns_diff_error"
    )
  )
}

# -----------------------------------------------------------------------------
# Benchmark core
# -----------------------------------------------------------------------------

run_black_box_nns_benchmark <- function(
    test_cases,
    noise_levels = c(0, 1e-4, 1e-3, 1e-2),
    R = 100,
    N_evals = 1,
    fd_grid = 10^seq(-0.5, -4.5, length.out = 30),
    rich_h = 0.1,
    rich_order = 6,
    nns_h = 0.1,
    nns_tol = 1e-10,
    nns_max_iter = 1000L,
    delta_pred = 1e-2,
    sign_delta = 1e-3,
    seed = 123) {
  
  set.seed(seed)
  
  summary_rows <- list()
  detail_rows <- list()
  
  for (tc in test_cases) {
    for (sigma in noise_levels) {
      
      cat(sprintf("Testing %-28s | point = %-8.4g | sigma = %-8g\n",
                  tc$label, tc$pt, sigma))
      
      x0 <- tc$pt
      true_deriv <- tc$df(x0)  # may be NA at kinks
      true_sign <- tc$dir_sign(x0, sign_delta)
      
      # storage
      nns_proj <- rep(NA_real_, R)
      nns_init <- rep(NA_real_, R)
      nns_inf  <- rep(NA_real_, R)
      rich_est <- rep(NA_real_, R)
      
      nns_hs <- rep(NA_real_, R)
      nns_iter <- rep(NA_real_, R)
      nns_term <- rep(NA_real_, R)
      nns_fail <- rep(FALSE, R)
      
      fd_estimates <- matrix(NA_real_, nrow = R, ncol = length(fd_grid))
      
      for (r in seq_len(R)) {
        rep_seed <- seed + 100000L * r + 1000L * round(1000 * sigma)
        
        f_rep <- make_frozen_noisy_function(
          fun = tc$f,
          sigma = sigma,
          rep_seed = rep_seed,
          N_evals = N_evals
        )
        
        out <- safe_nns_diff(
          f = f_rep,
          point = x0,
          h = nns_h,
          tol = nns_tol,
          max.iter = nns_max_iter,
          digits = 12,
          print.trace = FALSE
        )
        
        if (inherits(out, "nns_diff_error")) {
          nns_fail[r] <- TRUE
        } else {
          ext <- extract_nns_rows(out)
          nns_proj[r] <- ext$derivative
          nns_init[r] <- ext$fin_init_avg
          nns_inf[r]  <- ext$fin_inf_avg
          nns_hs[r]   <- ext$inferred_h
          nns_iter[r] <- ext$iterations
          nns_term[r] <- ext$termination.code
        }
        
        rich_est[r] <- tryCatch(
          richardson_diff(f_rep, x0, h = rich_h, order = rich_order),
          error = function(e) NA_real_
        )
        
        for (j in seq_along(fd_grid)) {
          fd_estimates[r, j] <- centered_fd(f_rep, x0, fd_grid[j])
        }
        
        detail_rows[[length(detail_rows) + 1L]] <- data.frame(
          Function = tc$label,
          Point = x0,
          Sigma = sigma,
          N_evals = N_evals,
          Replication = r,
          
          True_Derivative = true_deriv,
          True_Directional_Sign = true_sign,
          
          NNS_Proj = nns_proj[r],
          NNS_FinInit = nns_init[r],
          NNS_FinInf = nns_inf[r],
          Richardson = rich_est[r],
          
          NNS_Inferred_h = nns_hs[r],
          NNS_Iterations = nns_iter[r],
          NNS_Termination = nns_term[r],
          NNS_Failed = nns_fail[r],
          
          stringsAsFactors = FALSE
        )
      }
      
      # Oracle FD by derivative RMSE when derivative exists,
      # otherwise by tangent prediction error
      if (!is.na(true_deriv) && is.finite(true_deriv)) {
        fd_mse <- colMeans((fd_estimates - true_deriv)^2, na.rm = TRUE)
      } else {
        fd_mse <- sapply(seq_along(fd_grid), function(j) {
          mean(vapply(seq_len(R), function(r) {
            rep_seed <- seed + 100000L * r + 1000L * round(1000 * sigma)
            f_rep <- make_frozen_noisy_function(tc$f, sigma, rep_seed, N_evals)
            pe <- tangent_pred_error(f_rep, x0, fd_estimates[r, j], delta_pred)
            ifelse(is.na(pe), NaN, pe^2)
          }, numeric(1)), na.rm = TRUE)
        })
      }
      
      oracle_idx <- which.min(fd_mse)
      oracle_h <- fd_grid[oracle_idx]
      oracle_fd <- fd_estimates[, oracle_idx]
      
      # Metrics
      metric_block <- function(est_vec, method_name) {
        deriv_rrmse <- if (!is.na(true_deriv) && is.finite(true_deriv)) {
          sqrt(mean(vapply(est_vec, safe_rel_err, numeric(1), truth = true_deriv)^2, na.rm = TRUE))
        } else {
          NA_real_
        }
        
        deriv_rmse <- if (!is.na(true_deriv) && is.finite(true_deriv)) {
          sqrt(mean((est_vec - true_deriv)^2, na.rm = TRUE))
        } else {
          NA_real_
        }
        
        sign_acc <- mean(vapply(est_vec, function(z) {
          sz <- safe_sign(z)
          if (is.na(sz) || is.na(true_sign)) return(NA_real_)
          as.numeric(sz == true_sign)
        }, numeric(1)), na.rm = TRUE)
        
        pred_mae <- mean(vapply(seq_along(est_vec), function(i) {
          rep_seed <- seed + 100000L * i + 1000L * round(1000 * sigma)
          f_rep <- make_frozen_noisy_function(tc$f, sigma, rep_seed, N_evals)
          tangent_pred_error(f_rep, x0, est_vec[i], delta_pred)
        }, numeric(1)), na.rm = TRUE)
        
        est_sd <- sd(est_vec, na.rm = TRUE)
        
        c(
          Deriv_RRMSE = deriv_rrmse,
          Deriv_RMSE = deriv_rmse,
          Sign_Accuracy = sign_acc,
          Tangent_MAE = pred_mae,
          Estimate_SD = est_sd
        )
      }
      
      nns_proj_metrics <- metric_block(nns_proj, "NNS_Proj")
      nns_init_metrics <- metric_block(nns_init, "NNS_FinInit")
      nns_inf_metrics  <- metric_block(nns_inf,  "NNS_FinInf")
      rich_metrics     <- metric_block(rich_est, "Richardson")
      oracle_metrics   <- metric_block(oracle_fd, "OracleFD")
      
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        Function = tc$label,
        Category = tc$category,
        Point = x0,
        Sigma = sigma,
        N_evals = N_evals,
        R = R,
        
        True_Derivative = true_deriv,
        True_Directional_Sign = true_sign,
        
        NNS_Proj_RRMSE = nns_proj_metrics["Deriv_RRMSE"],
        NNS_FinInit_RRMSE = nns_init_metrics["Deriv_RRMSE"],
        NNS_FinInf_RRMSE = nns_inf_metrics["Deriv_RRMSE"],
        Richardson_RRMSE = rich_metrics["Deriv_RRMSE"],
        OracleFD_RRMSE = oracle_metrics["Deriv_RRMSE"],
        
        NNS_Proj_SignAcc = nns_proj_metrics["Sign_Accuracy"],
        NNS_FinInit_SignAcc = nns_init_metrics["Sign_Accuracy"],
        NNS_FinInf_SignAcc = nns_inf_metrics["Sign_Accuracy"],
        Richardson_SignAcc = rich_metrics["Sign_Accuracy"],
        OracleFD_SignAcc = oracle_metrics["Sign_Accuracy"],
        
        NNS_Proj_TangentMAE = nns_proj_metrics["Tangent_MAE"],
        NNS_FinInit_TangentMAE = nns_init_metrics["Tangent_MAE"],
        NNS_FinInf_TangentMAE = nns_inf_metrics["Tangent_MAE"],
        Richardson_TangentMAE = rich_metrics["Tangent_MAE"],
        OracleFD_TangentMAE = oracle_metrics["Tangent_MAE"],
        
        NNS_Proj_SD = nns_proj_metrics["Estimate_SD"],
        NNS_FinInit_SD = nns_init_metrics["Estimate_SD"],
        NNS_FinInf_SD = nns_inf_metrics["Estimate_SD"],
        Richardson_SD = rich_metrics["Estimate_SD"],
        OracleFD_SD = oracle_metrics["Estimate_SD"],
        
        OracleFD_h = oracle_h,
        Median_NNS_h = median(abs(nns_hs), na.rm = TRUE),
        Median_NNS_Iterations = median(nns_iter, na.rm = TRUE),
        Median_Termination_Code = median(nns_term, na.rm = TRUE),
        NNS_Failure_Rate = mean(nns_fail, na.rm = TRUE),
        
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    summary = do.call(rbind, summary_rows),
    details = do.call(rbind, detail_rows),
    fd_grid = fd_grid
  )
}

# -----------------------------------------------------------------------------
# Test cases: non-analytic, piecewise, threshold, quantized, black-box style
# -----------------------------------------------------------------------------

test_cases <- list(
  list(
    label = "abs(x) smooth side",
    category = "piecewise_continuous",
    f = function(x) abs(x),
    df = function(x) if (x == 0) NA_real_ else sign(x),
    dir_sign = function(x, delta) safe_sign(abs(x + delta) - abs(x - delta)),
    pt = 0.3
  ),
  list(
    label = "abs(x) at kink",
    category = "piecewise_continuous",
    f = function(x) abs(x),
    df = function(x) NA_real_,
    dir_sign = function(x, delta) safe_sign(abs(x + delta) - abs(x - delta)),
    pt = 0
  ),
  list(
    label = "ReLU smooth side",
    category = "piecewise_continuous",
    f = function(x) max(0, x),
    df = function(x) if (x < 0) 0 else if (x > 0) 1 else NA_real_,
    dir_sign = function(x, delta) safe_sign(max(0, x + delta) - max(0, x - delta)),
    pt = 0.3
  ),
  list(
    label = "ReLU at kink",
    category = "piecewise_continuous",
    f = function(x) max(0, x),
    df = function(x) NA_real_,
    dir_sign = function(x, delta) safe_sign(max(0, x + delta) - max(0, x - delta)),
    pt = 0
  ),
  list(
    label = "clipped linear interior",
    category = "saturation",
    f = function(x) min(max(x, -1), 1),
    df = function(x) if (x > -1 && x < 1) 1 else if (x < -1 || x > 1) 0 else NA_real_,
    dir_sign = function(x, delta) safe_sign(min(max(x + delta, -1), 1) - min(max(x - delta, -1), 1)),
    pt = 0.2
  ),
  list(
    label = "clipped linear threshold",
    category = "saturation",
    f = function(x) min(max(x, -1), 1),
    df = function(x) NA_real_,
    dir_sign = function(x, delta) safe_sign(min(max(x + delta, -1), 1) - min(max(x - delta, -1), 1)),
    pt = 1
  ),
  list(
    label = "indicator threshold",
    category = "discontinuous",
    f = function(x) as.numeric(x > 0),
    df = function(x) NA_real_,
    dir_sign = function(x, delta) safe_sign(as.numeric(x + delta > 0) - as.numeric(x - delta > 0)),
    pt = 0
  ),
  list(
    label = "rounded quadratic",
    category = "quantized",
    f = function(x) round(x^2, 2),
    df = function(x) NA_real_,   # quantization makes classical derivative unreliable
    dir_sign = function(x, delta) safe_sign(round((x + delta)^2, 2) - round((x - delta)^2, 2)),
    pt = 1.2
  ),
  list(
    label = "MC option payoff style",
    category = "black_box_sim",
    f = function(x) {
      # branching, real-only, discontinuous payoff structure
      mean(pmax(x + rnorm(100), 0))
    },
    df = function(x) NA_real_,
    dir_sign = function(x, delta) NA_real_,
    pt = 0
  )
)

# -----------------------------------------------------------------------------
# Run benchmark
# -----------------------------------------------------------------------------

res <- run_black_box_nns_benchmark(
  test_cases = test_cases,
  noise_levels = c(0, 1e-4, 1e-3, 1e-2),
  R = 100,
  N_evals = 1,
  fd_grid = 10^seq(-0.5, -4.5, length.out = 30),
  rich_h = 0.1,
  rich_order = 6,
  nns_h = 0.1,
  nns_tol = 1e-10,
  nns_max_iter = 1000L,
  delta_pred = 1e-2,
  sign_delta = 1e-3,
  seed = 123
)

# -----------------------------------------------------------------------------
# Display results
# -----------------------------------------------------------------------------

print(res$summary)

cat("\n====================================================\n")
cat("Median summary by category and sigma\n")
cat("====================================================\n")

summary_by_group <- aggregate(
  cbind(
    NNS_Proj_RRMSE,
    NNS_FinInit_RRMSE,
    NNS_FinInf_RRMSE,
    Richardson_RRMSE,
    OracleFD_RRMSE,
    NNS_Proj_SignAcc,
    NNS_FinInit_SignAcc,
    NNS_FinInf_SignAcc,
    Richardson_SignAcc,
    OracleFD_SignAcc,
    NNS_Proj_TangentMAE,
    NNS_FinInit_TangentMAE,
    NNS_FinInf_TangentMAE,
    Richardson_TangentMAE,
    OracleFD_TangentMAE,
    Median_NNS_h,
    Median_NNS_Iterations,
    Median_Termination_Code,
    NNS_Failure_Rate
  ) ~ Category + Sigma,
  data = res$summary,
  FUN = median,
  na.rm = TRUE
)

print(summary_by_group)

cat("\n====================================================\n")
cat("Top rows of replication-level details\n")
cat("====================================================\n")

print(head(res$details, 20))