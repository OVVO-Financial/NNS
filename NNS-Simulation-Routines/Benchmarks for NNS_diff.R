# =============================================================================
# BENCHMARK FOR UPDATED NNS.diff
# =============================================================================

suppressPackageStartupMessages({
  library(NNS)
})

rel_err <- function(est, truth) {
  if (length(est) == 0 || is.na(est) || is.nan(est) || is.infinite(est)) return(NA_real_)
  if (abs(truth) < 1e-12) return(abs(est - truth))
  abs((est - truth) / truth)
}

make_frozen_noisy_function <- function(fun, sigma, rep_seed, N_evals = 1) {
  force(fun); force(sigma); force(rep_seed); force(N_evals)
  
  function(x) {
    if (is.complex(x)) {
      return(fun(x))
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
      mean(replicate(N_evals, fun(xx) + rnorm(1, 0, sigma)))
    }, numeric(1))
    
    if (length(vals) == 1) vals[[1]] else vals
  }
}

richardson_diff <- function(f, point, h = 0.1, order = 6) {
  Rmat <- matrix(NA_real_, nrow = order, ncol = order)
  
  for (i in seq_len(order)) {
    hi <- h / (2^(i - 1))
    Rmat[i, 1] <- (f(point + hi) - f(point - hi)) / (2 * hi)
  }
  
  for (j in 2:order) {
    for (i in j:order) {
      Rmat[i, j] <- Rmat[i, j - 1] +
        (Rmat[i, j - 1] - Rmat[i - 1, j - 1]) / (4^(j - 1) - 1)
    }
  }
  
  Rmat[order, order]
}

extract_nns_rows <- function(out) {
  blank <- list(
    derivative = NA_real_,
    inferred_h = NA_real_,
    iterations = NA_real_,
    converged = NA_real_,
    termination.code = NA_real_,
    fin_init_avg = NA_real_,
    fin_inf_avg = NA_real_,
    cplx_inf = NA_real_
  )
  
  if (is.null(out) || !(is.matrix(out) || is.data.frame(out))) return(blank)
  
  rn <- rownames(out)
  if (is.null(rn)) return(blank)
  
  vals <- suppressWarnings(as.numeric(out[, 1]))
  names(vals) <- rn
  
  get1 <- function(name) {
    if (name %in% names(vals)) as.numeric(vals[name][1]) else NA_real_
  }
  
  list(
    derivative = get1("DERIVATIVE"),
    inferred_h = get1("Inferred h"),
    iterations = get1("iterations"),
    converged = get1("converged"),
    termination.code = get1("termination.code"),
    fin_init_avg = get1("Initial h averaged finite step.Averaged Finite Step"),
    fin_inf_avg = get1("Inferred h averaged finite step.Averaged Finite Step"),
    cplx_inf = get1("Complex Step Derivative (Inferred h)")
  )
}

safe_nns_diff <- function(f, point, h = 0.1, tol = 1e-10, max.iter = 1000L,
                          digits = 12, print.trace = FALSE, plot = FALSE) {
  tryCatch(
    NNS.diff(
      f = f,
      point = point,
      h = h,
      tol = tol,
      max.iter = max.iter,
      digits = digits,
      print.trace = print.trace,
      plot = plot
    ),
    error = function(e) structure(
      list(error_message = conditionMessage(e)),
      class = "nns_diff_error"
    )
  )
}

run_nns_test_suite <- function(functions,
                               noise_levels,
                               R = 100,
                               N_evals = 1,
                               h_grid = 10^seq(-0.5, -4.5, length.out = 30),
                               rich_h = 0.1,
                               rich_order = 6,
                               nns_h = 0.1,
                               nns_tol = 1e-10,
                               nns_max_iter = 1000L,
                               nns_digits = 12,
                               seed = 123) {
  
  set.seed(seed)
  
  summary_rows <- list()
  detail_rows <- list()
  
  for (fn_info in functions) {
    for (sigma in noise_levels) {
      
      cat(sprintf("Testing %s at sigma = %g with N_evals = %d ...\n",
                  fn_info$label, sigma, N_evals))
      
      point <- fn_info$pt
      true_deriv <- fn_info$df(point)
      
      err_nns_proj <- rep(NA_real_, R)
      err_nns_fin_init <- rep(NA_real_, R)
      err_nns_fin_inf <- rep(NA_real_, R)
      err_nns_cplx_inf <- rep(NA_real_, R)
      err_rich <- rep(NA_real_, R)
      
      nns_h_root <- rep(NA_real_, R)
      nns_iter <- rep(NA_real_, R)
      nns_conv <- rep(NA_real_, R)
      nns_term <- rep(NA_real_, R)
      nns_fail <- rep(FALSE, R)
      nns_error_message <- rep(NA_character_, R)
      
      fd_rel_errs <- matrix(NA_real_, nrow = R, ncol = length(h_grid))
      
      for (r in seq_len(R)) {
        rep_seed <- seed + 100000L * r + 1000L * round(1000 * sigma)
        
        f_rep <- make_frozen_noisy_function(
          fun = fn_info$f,
          sigma = sigma,
          rep_seed = rep_seed,
          N_evals = N_evals
        )
        
        out <- safe_nns_diff(
          f = f_rep,
          point = point,
          h = nns_h,
          tol = nns_tol,
          max.iter = nns_max_iter,
          digits = nns_digits,
          print.trace = FALSE,
          plot = FALSE
        )
        
        if (inherits(out, "nns_diff_error")) {
          nns_fail[r] <- TRUE
          nns_error_message[r] <- out$error_message
        } else {
          ext <- extract_nns_rows(out)
          
          err_nns_proj[r]     <- rel_err(ext$derivative, true_deriv)
          err_nns_fin_init[r] <- rel_err(ext$fin_init_avg, true_deriv)
          err_nns_fin_inf[r]  <- rel_err(ext$fin_inf_avg, true_deriv)
          err_nns_cplx_inf[r] <- rel_err(ext$cplx_inf, true_deriv)
          
          nns_h_root[r] <- ext$inferred_h
          nns_iter[r]   <- ext$iterations
          nns_conv[r]   <- ext$converged
          nns_term[r]   <- ext$termination.code
        }
        
        rich_hat <- tryCatch(
          richardson_diff(f_rep, point, h = rich_h, order = rich_order),
          error = function(e) NA_real_
        )
        err_rich[r] <- rel_err(rich_hat, true_deriv)
        
        for (j in seq_along(h_grid)) {
          h <- h_grid[j]
          fd_hat <- tryCatch(
            (f_rep(point + h) - f_rep(point - h)) / (2 * h),
            error = function(e) NA_real_
          )
          fd_rel_errs[r, j] <- rel_err(fd_hat, true_deriv)
        }
        
        detail_rows[[length(detail_rows) + 1L]] <- data.frame(
          Function = fn_info$label,
          Sigma = sigma,
          N_evals = N_evals,
          Replication = r,
          True_Derivative = true_deriv,
          NNS_Proj_RelErr = err_nns_proj[r],
          NNS_FinInit_RelErr = err_nns_fin_init[r],
          NNS_FinInf_RelErr = err_nns_fin_inf[r],
          NNS_CplxInf_RelErr = err_nns_cplx_inf[r],
          Richardson_RelErr = err_rich[r],
          NNS_Inferred_h = nns_h_root[r],
          NNS_Iterations = nns_iter[r],
          NNS_Converged = nns_conv[r],
          NNS_Termination = nns_term[r],
          NNS_Failed = nns_fail[r],
          NNS_Error = nns_error_message[r],
          stringsAsFactors = FALSE
        )
      }
      
      fd_mse <- colMeans(fd_rel_errs^2, na.rm = TRUE)
      oracle_idx <- which.min(fd_mse)
      oracle_h <- h_grid[oracle_idx]
      oracle_rrmse <- sqrt(fd_mse[oracle_idx])
      
      trial_best_h <- apply(fd_rel_errs, 1, function(x) {
        if (all(is.na(x))) return(NA_real_)
        h_grid[which.min(x)]
      })
      
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        Function = fn_info$label,
        Point = point,
        Sigma = sigma,
        N_evals = N_evals,
        R = R,
        NNS_Proj_RRMSE = sqrt(mean(err_nns_proj^2, na.rm = TRUE)),
        NNS_FinInit_RRMSE = sqrt(mean(err_nns_fin_init^2, na.rm = TRUE)),
        NNS_FinInf_RRMSE = sqrt(mean(err_nns_fin_inf^2, na.rm = TRUE)),
        NNS_CplxInf_RRMSE = sqrt(mean(err_nns_cplx_inf^2, na.rm = TRUE)),
        Richardson_RRMSE = sqrt(mean(err_rich^2, na.rm = TRUE)),
        OracleFD_RRMSE = oracle_rrmse,
        Median_NNS_h = median(abs(nns_h_root), na.rm = TRUE),
        Median_NNS_Iterations = median(nns_iter, na.rm = TRUE),
        Median_Trial_Best_h = median(trial_best_h, na.rm = TRUE),
        OracleFD_h = oracle_h,
        NNS_Failure_Rate = mean(nns_fail, na.rm = TRUE),
        Median_Termination_Code = median(nns_term, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    summary = do.call(rbind, summary_rows),
    details = do.call(rbind, detail_rows),
    h_grid = h_grid
  )
}

test_funs <- list(
  list(f = sin, df = cos, pt = 1.0, label = "sin(x) at x = 1"),
  list(f = exp, df = exp, pt = 1.0, label = "exp(x) at x = 1"),
  list(f = function(x) x^3, df = function(x) 3 * x^2, pt = 2.0, label = "x^3 at x = 2")
)

sigmas <- c(0, 1e-4, 1e-3, 1e-2)

res <- run_nns_test_suite(
  functions = test_funs,
  noise_levels = sigmas,
  R = 100,
  N_evals = 1,
  h_grid = 10^seq(-0.5, -4.5, length.out = 30),
  rich_h = 0.1,
  rich_order = 6,
  nns_h = 0.1,
  nns_tol = 1e-10,
  nns_max_iter = 1000L,
  nns_digits = 12,
  seed = 123
)

print(res$summary)
print(head(res$details, 12))