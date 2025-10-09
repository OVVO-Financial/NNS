#' NNS Seasonality Test
#'
#' Seasonality test based on the coefficient of variation for the variable and lagged component series.  A result of 1 signifies no seasonality present.
#'
#' @param variable a numeric vector.
#' @param modulo integer(s); NULL (default) Used to find the nearest multiple(s) in the reported seasonal period.
#' @param mod.only logical; \code{TRUE} (default) Limits the number of seasonal periods returned to the specified \code{modulo}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference.
#' @return Returns a matrix of all periods exhibiting less coefficient of variation than the variable with \code{"all.periods"}; and the single period exhibiting the least coefficient of variation versus the variable with \code{"best.period"}; as well as a vector of \code{"periods"} for easy call into \link{NNS.ARMA.optim}.  If no seasonality is detected, \code{NNS.seas} will return ("No Seasonality Detected").
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## To call strongest period based on coefficient of variation:
#' NNS.seas(x, plot = FALSE)$best.period
#'
#' ## Using modulos for logical seasonal inference:
#' NNS.seas(x, modulo = c(2,3,5,7), plot = FALSE)
#' }
#' @export



NNS.seas <- function(variable,
                     modulo = NULL,
                     mod.only = TRUE,
                     plot = TRUE) {
  
  # Coerce to vector if tbl or data.table
  if (any(class(variable) %in% c("tbl", "data.table"))) {
    variable <- as.vector(unlist(variable, use.names = FALSE))
  }
  
  # Validate input
  if (!is.numeric(variable)) stop("Variable must be numeric")
  if (anyNA(variable)) stop("You have some missing values, please address.")
  if (any(is.infinite(variable))) stop("Infinite values not allowed")
  
  n <- length(variable)
  if (n < 5L) {
    return(data.table::data.table(
      "Period" = 0L,
      "Coefficient.of.Variation" = 0,
      "Variable.Coefficient.of.Variation" = 0,
      key = "Coefficient.of.Variation"
    ))
  }
  
  # If no modulo provided, do not restrict to modulo 
  if (is.null(modulo)) mod.only <- FALSE
  
  # Precompute trimmed variants used throughout
  variable_1 <- variable[1L:(n - 1L)]
  n1 <- n - 1L
  variable_2 <- if (n1 >= 2L) variable_1[1L:(n1 - 1L)] else numeric(0L)
  n2 <- length(variable_2)
  
  half_n <- n %/% 2L
  
  # Global reference: CV(variable); if mean==0, fallback to inverse |ACF(1)|
  mean_var <- mean(variable)
  use_cv <- (mean_var != 0)
  if (use_cv) {
    var_cov <- abs(stats::sd(variable) / mean_var)
  } else {
    a1 <- stats::acf(variable, lag.max = 1L, plot = FALSE)$acf[2L]
    var_cov <- abs(a1)^-1
  }
  
  # Helper: CV for subseries with robust fallback
  cv_or_fallback <- function(x) {
    if (length(x) < 2L) return(var_cov)
    if (use_cv) {
      z <- abs(stats::sd(x) / mean(x))
    } else {
      a1 <- stats::acf(x, lag.max = 1L, plot = FALSE)$acf[2L]
      z <- abs(a1)^-1
    }
    if (!is.finite(z)) var_cov else z
  }
  
  # Preallocate
  out  <- numeric(half_n); out1 <- numeric(half_n); out2 <- numeric(half_n)
  inst <- integer(half_n); inst1 <- integer(half_n); inst2 <- integer(half_n)
  
  # MAIN LOOP: reverse-stepped subseries per candidate period
  for (i in 1L:half_n) {
    idx  <- seq.int(n,  1L, by = -i)
    idx1 <- seq.int(n1, 1L, by = -i)
    idx2 <- if (n2) seq.int(n2, 1L, by = -i) else integer(0L)
    
    t  <- cv_or_fallback(variable[idx])
    t1 <- cv_or_fallback(variable_1[idx1])
    t2 <- cv_or_fallback(variable_2[idx2])
    
    # Accept if subseries CV <= variable CV
    if (t  <= var_cov) { inst[i]  <- i; out[i]  <- t  }
    if (t1 <= var_cov) { inst1[i] <- i; out1[i] <- t1 }
    if (t2 <= var_cov) { inst2[i] <- i; out2[i] <- t2 }
  }

  
  # Require all three staggered subseries to pass; average their CVs for stability
  keep <- (inst > 0L) & (inst1 > 0L) & (inst2 > 0L)
  
  if (any(keep)) {
    cv_mean <- rowMeans(cbind(out[keep], out1[keep], out2[keep]))
    periods <- inst[keep]
    M <- data.table::data.table(
      "Period" = periods,
      "Coefficient.of.Variation" = cv_mean,
      "Variable.Coefficient.of.Variation" = rep(var_cov, length(periods)),
      key = "Coefficient.of.Variation"
    )
  } else {
    # Reference: when nothing qualifies, return a single row at var_cov
    M <- data.table::data.table(
      "Period" = 1L,
      "Coefficient.of.Variation" = var_cov,
      "Variable.Coefficient.of.Variation" = var_cov,
      key = "Coefficient.of.Variation"
    )
  }
  
  # --- Modulo constraint handling ---
  if (!is.null(modulo)) {
    a <- M[["Period"]]
    plus  <- a + (modulo - a %% modulo)
    minus <- a - (a %% modulo)
    periods <- unique(c(minus, plus))
    periods <- periods[!is.na(periods) & periods > 0L]
    
    if (mod.only) {
      # Keep modulo-constrained periods AND append any missing multiples at var_cov
      keep_dt <- M[Period %in% periods]
      add_periods <- setdiff(periods, keep_dt[["Period"]])
      if (length(add_periods)) {
        add_dt <- data.table::data.table(
          "Period" = add_periods,
          "Coefficient.of.Variation" = rep(var_cov, length(add_periods)),
          "Variable.Coefficient.of.Variation" = rep(var_cov, length(add_periods))
        )
        keep_dt <- data.table::rbindlist(list(keep_dt, add_dt), use.names = TRUE)
      }
      # If still empty, fall back to a single row
      if (nrow(keep_dt) == 0L) {
        keep_dt <- data.table::data.table(
          "Period" = 1L,
          "Coefficient.of.Variation" = var_cov,
          "Variable.Coefficient.of.Variation" = var_cov
        )
      }
      data.table::setkey(keep_dt, "Coefficient.of.Variation")
      M <- keep_dt
    } else {
      # Ensure period 1 is present if it wasn't already
      if (!1L %in% M[["Period"]]) periods <- unique(c(periods, 1L))
      # Append missing modulo-suggested periods at var_cov
      add_periods <- setdiff(periods, M[["Period"]])
      if (length(add_periods) > 0L) {
        add_dt <- data.table::data.table(
          "Period" = add_periods,
          "Coefficient.of.Variation" = rep(var_cov, length(add_periods)),
          "Variable.Coefficient.of.Variation" = rep(var_cov, length(add_periods))
        )
        M <- data.table::rbindlist(list(M, add_dt), use.names = TRUE)
        data.table::setkey(M, "Coefficient.of.Variation")
      }
    }
  }
  
  # Honor original constraint: only periods < n/2
  M <- M[Period < n/2]
  
  # Plot (diagnostic)
  if (isTRUE(plot)) {
    overall_cv <- if (use_cv) abs(stats::sd(variable) / mean_var) else var_cov
    # Color mapping; guard division by zero
    if (is.finite(overall_cv) && overall_cv > 0) {
      predictive_strength <- pmax(pmin(1 - (M[["Coefficient.of.Variation"]] / overall_cv), 1), -1)
      color_ramp <- grDevices::colorRampPalette(c("blue", "green"))(100)
      point_colors <- color_ramp[cut(predictive_strength, breaks = 100, labels = FALSE)]
    } else {
      point_colors <- "black"
    }
    
    plot(M[["Period"]], M[["Coefficient.of.Variation"]],
         xlab = "Period",
         ylab = "Component Series CV",
         main = "Seasonality Detection via Predictive Power\n(Lower CV = Tighter Distribution = More Predictable)",
         ylim = c(0, 2 * if (is.finite(overall_cv)) overall_cv else 1),
         col = point_colors, pch = 19)
    
    # Highlight best period
    points(M[["Period"]][1L], M[["Coefficient.of.Variation"]][1L], pch = 19, col = "red", cex = 1.5)
    
    # Reference CV line and centered label
    abline(h = overall_cv, col = "red", lty = 2)
    xmid <- mean(par("usr")[1:2])              # center of current x-axis
    graphics::text(xmid, overall_cv,
                   labels = "Overall Series CV\n(Predictive Power Threshold)",
                   adj = c(0.5, 0.5), col = "red", xpd = NA)
  }
  
  # Return results
  list(
    "all.periods" = M,
    "best.period" = M[["Period"]][1L],
    "periods"     = as.integer(M[["Period"]])
  )
}



