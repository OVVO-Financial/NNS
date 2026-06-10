# ==============================================================================
# NNS PARTIAL MOMENTS SIMPLE BEFORE / AFTER TEST
# ==============================================================================
# Save this OUTSIDE the R/ folder, for example:
#   NNS_before_after_test.R
#
# Run before fix:
#   phase <- "before"
#   source("NNS_before_after_test.R")
#
# Install fixed NNS, restart R, then run:
#   phase <- "after"
#   source("NNS_before_after_test.R")
# ==============================================================================

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  install.packages("microbenchmark")
}
library(microbenchmark)

if (!exists("phase")) {
  phase <- "before"
}

pm_fun <- function(name) {
  getFromNamespace(name, "NNS")
}

LPM       <- pm_fun("LPM")
UPM       <- pm_fun("UPM")
LPM_ratio <- pm_fun("LPM.ratio")
UPM_ratio <- pm_fun("UPM.ratio")
NNS_CDF   <- pm_fun("NNS.CDF")

set.seed(123)

n_accuracy <- 1000
n_speed    <- 5000

x_acc <- rnorm(n_accuracy)
targets_acc <- sort(x_acc)

x_speed <- rnorm(n_speed)
targets_speed <- sort(x_speed)

degrees <- c(0, 1, 2, 3)

verify_and_report <- function(test_name, before, after, tolerance = 1e-7) {
  cat("\n", paste0(rep("-", 60), collapse = ""), "\n")
  cat(" Accuracy:", test_name, "\n")
  cat(paste0(rep("-", 60), collapse = ""), "\n")
  
  check <- all.equal(before, after, tolerance = tolerance, check.attributes = FALSE)
  
  if (isTRUE(check)) {
    cat("  [PASS] Before and after outputs match.\n")
  } else {
    cat("  [FAIL/WARNING] Difference found:\n")
    print(check)
  }
}

# ------------------------------------------------------------------------------
# BEFORE
# ------------------------------------------------------------------------------

if (phase == "before") {
  
  cat("\n============================================================\n")
  cat(" RUNNING BEFORE BASELINE\n")
  cat("============================================================\n")
  
  before_results <- list()
  
  for (d in degrees) {
    before_results[[paste0("LPM_", d)]]       <- LPM(d, targets_acc, x_acc)
    before_results[[paste0("UPM_", d)]]       <- UPM(d, targets_acc, x_acc)
    before_results[[paste0("LPM_ratio_", d)]] <- LPM_ratio(d, targets_acc, x_acc)
    before_results[[paste0("UPM_ratio_", d)]] <- UPM_ratio(d, targets_acc, x_acc)
  }
  
  before_results[["NNS_CDF_0"]] <- NNS_CDF(x_acc, 0, plot = FALSE)
  before_results[["NNS_CDF_1"]] <- NNS_CDF(x_acc, 1, plot = FALSE)
  before_results[["NNS_CDF_2"]] <- NNS_CDF(x_acc, 2, plot = FALSE)
  
  saveRDS(before_results, "NNS_before_accuracy_results.rds")
  
  cat("\nSaved accuracy baseline: NNS_before_accuracy_results.rds\n")
  
  cat("\n============================================================\n")
  cat(" BEFORE TIMINGS\n")
  cat("============================================================\n")
  
  bench_before <- microbenchmark(
    LPM_0       = LPM(0, targets_speed, x_speed),
    LPM_1       = LPM(1, targets_speed, x_speed),
    UPM_1       = UPM(1, targets_speed, x_speed),
    LPM_ratio_1 = LPM_ratio(1, targets_speed, x_speed),
    UPM_ratio_1 = UPM_ratio(1, targets_speed, x_speed),
    NNS_CDF_1   = NNS_CDF(x_speed, 1, plot = FALSE),
    times = 10
  )
  
  print(bench_before)
  saveRDS(bench_before, "NNS_before_timings.rds")
  
  cat("\nSaved timing baseline: NNS_before_timings.rds\n")
  cat("\nNow install the fixed NNS, restart R, set phase <- 'after', and rerun.\n")
}

# ------------------------------------------------------------------------------
# AFTER
# ------------------------------------------------------------------------------

if (phase == "after") {
  
  cat("\n============================================================\n")
  cat(" RUNNING AFTER TEST\n")
  cat("============================================================\n")
  
  before_results <- readRDS("NNS_before_accuracy_results.rds")
  before_timings <- readRDS("NNS_before_timings.rds")
  
  after_results <- list()
  
  for (d in degrees) {
    after_results[[paste0("LPM_", d)]]       <- LPM(d, targets_acc, x_acc)
    after_results[[paste0("UPM_", d)]]       <- UPM(d, targets_acc, x_acc)
    after_results[[paste0("LPM_ratio_", d)]] <- LPM_ratio(d, targets_acc, x_acc)
    after_results[[paste0("UPM_ratio_", d)]] <- UPM_ratio(d, targets_acc, x_acc)
  }
  
  after_results[["NNS_CDF_0"]] <- NNS_CDF(x_acc, 0, plot = FALSE)
  after_results[["NNS_CDF_1"]] <- NNS_CDF(x_acc, 1, plot = FALSE)
  after_results[["NNS_CDF_2"]] <- NNS_CDF(x_acc, 2, plot = FALSE)
  
  cat("\n============================================================\n")
  cat(" ACCURACY CHECKS\n")
  cat("============================================================\n")
  
  for (nm in names(before_results)) {
    verify_and_report(nm, before_results[[nm]], after_results[[nm]])
  }
  
  cat("\n============================================================\n")
  cat(" AFTER TIMINGS\n")
  cat("============================================================\n")
  
  bench_after <- microbenchmark(
    LPM_0       = LPM(0, targets_speed, x_speed),
    LPM_1       = LPM(1, targets_speed, x_speed),
    UPM_1       = UPM(1, targets_speed, x_speed),
    LPM_ratio_1 = LPM_ratio(1, targets_speed, x_speed),
    UPM_ratio_1 = UPM_ratio(1, targets_speed, x_speed),
    NNS_CDF_1   = NNS_CDF(x_speed, 1, plot = FALSE),
    times = 10
  )
  
  print(bench_after)
  saveRDS(bench_after, "NNS_after_timings.rds")
  
  cat("\n============================================================\n")
  cat(" BEFORE VS AFTER TIMING SUMMARY\n")
  cat("============================================================\n")
  
  before_summary <- summary(before_timings)
  after_summary  <- summary(bench_after)
  
  timing_compare <- data.frame(
    expr = before_summary$expr,
    before_median_ms = before_summary$median / 1e6,
    after_median_ms  = after_summary$median / 1e6,
    speedup = before_summary$median / after_summary$median
  )
  
  print(timing_compare, row.names = FALSE)
  
  write.csv(timing_compare, "NNS_before_after_timing_summary.csv", row.names = FALSE)
  
  cat("\nSaved timing comparison: NNS_before_after_timing_summary.csv\n")
}

cat("\nDone.\n")