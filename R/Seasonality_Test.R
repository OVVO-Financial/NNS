#' NNS Seasonality Test
#'
#' Seasonality test based on the coefficient of variation for the variable and lagged component series.  A linear trend is removed first (\code{detrend = TRUE}), keeping the mean of the variable, so that a trend does not mask the period.  For each lag \code{m}, the variable is split into its \code{m} component series \code{variable[seq(j, n, m)]}, \code{j = 1, ..., m}, and lag \code{m} is reported when the length-weighted mean of the component coefficients of variation is below the coefficient of variation of the variable.  Dispersion is the Gini mean difference, the degree-1 lower plus upper partial moment about each observation in turn, which is unbiased for any distribution at any component length.  The coefficient of variation is used when the (detrended) variable stays on one side of zero; otherwise dispersions are compared directly.  Lags are scanned up to \code{max.period}, by default \code{n / 5}, so every component series has at least 5 observations.  Periods added through \code{modulo} are kept up to \code{n / 2}.  On a series without seasonality about half of the lags pass by chance, so treat the result as candidate periods and validate them, for example with \link{NNS.ARMA.optim}, or set \code{p.value = TRUE}.  A result of 1 signifies no seasonality present.
#'
#' @param variable a numeric vector.
#' @param modulo integer(s); NULL (default) Adds the multiples of \code{modulo} on either side of each detected period (period 7 with \code{modulo = 5} gives 5 and 10).
#' @param mod.only logical; \code{TRUE} (default) Limits the number of seasonal periods returned to the specified \code{modulo}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference.
#' @param detrend logical; \code{TRUE} (default) removes the variable's linear trend, keeping its mean, before the lags are scored; an exactly linear variable is then constant and has no seasonal period.
#' @param max.period integer; \code{NULL} (default) scans lags up to \code{n / 5}.  A smaller value caps the scan, whose cost grows with \code{n} times the number of lags.
#' @param p.value logical; \code{FALSE} (default).  If \code{TRUE}, each period's dispersion ratio is compared with the best ratio over all lags in \code{nperm} shuffled copies of the (detrended) variable, which remove the original temporal ordering.  \code{all.periods} gains a \code{p.value} column, and \code{periods} and \code{best.period} keep only the periods with \code{p.value <= 0.05} (\code{1} if none); each copy contributes its best lag, accounting for scanning multiple lags under an exchangeable null. A serially dependent nonseasonal series need not satisfy that null.  Uses the random number generator.
#' @param dispersion \code{"auto"} (default) retains relative CV scoring on one-sided data and absolute Gini dispersion otherwise. \code{"absolute"} always compares Gini dispersions, by centering the series first; use this translation-invariant option for interval-scale data such as temperature.
#' @param nperm integer; 100 (default) number of shuffled copies when \code{p.value = TRUE}.
#' @return Returns a matrix of all periods exhibiting less coefficient of variation than the variable with \code{"all.periods"}; and the single period exhibiting the least coefficient of variation versus the variable with \code{"best.period"}; as well as a vector of \code{"periods"} for easy call into \link{NNS.ARMA.optim}.  If no seasonality is detected, \code{all.periods} holds a single row with \code{Period = 1} and a coefficient of variation equal to the variable's (\code{Period = 0} for fewer than 5 observations).  The coefficients of variation are those of the detrended variable when \code{detrend = TRUE}.
#' @note The historical column names \code{Coefficient.of.Variation} and \code{Variable.Coefficient.of.Variation} are retained in absolute mode; their values then contain absolute Gini dispersions.
#'
#' @note \code{best.period} is the strongest tested predictive lag, not an estimator guaranteed to return the smallest fundamental period; its multiples may also score well. The minimum five observations per component is a reliability default; \code{modulo} can supply longer domain-informed candidates. Shuffling removes serial dependence as well as seasonality, so the optional permutation comparison is not calibrated against every dependent nonseasonal process.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
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
#'
#' ## Keep only periods that beat shuffled copies of the series:
#' NNS.seas(x, plot = FALSE, p.value = TRUE)$best.period
#' }
#' @export



NNS.seas <- function(variable,
                     modulo = NULL,
                     mod.only = TRUE,
                     plot = TRUE,
                     detrend = TRUE,
                     max.period = NULL,
                     p.value = FALSE,
                     nperm = 100L,
                     dispersion = c("auto", "absolute")) {
  dispersion <- match.arg(dispersion)
  # API per NNS manual / Rd (arguments & defaults)  :contentReference[oaicite:3]{index=3}
  # Coerce tbl/data.table to numeric vector (repo reference)  :contentReference[oaicite:4]{index=4}
  if (any(class(variable) %in% c("tbl", "data.table"))) {
    variable <- as.vector(unlist(variable, use.names = FALSE))
  }
  if (!is.numeric(variable)) stop("Variable must be numeric")
  if (anyNA(variable)) stop("You have some missing values, please address.")
  if (any(is.infinite(variable))) stop("Infinite values not allowed")
  if (!is.null(max.period) && (!is.numeric(max.period) || length(max.period) != 1L ||
                               !is.finite(max.period) || max.period < 1))
    stop("max.period must be NULL or a single number >= 1.")
  if (isTRUE(p.value) && (!is.numeric(nperm) || length(nperm) != 1L || !is.finite(nperm) || nperm < 1))
    stop("nperm must be a single number >= 1.")

  x <- as.numeric(variable)
  if (dispersion == "absolute") x <- x - mean(x)
  n <- length(x)
  # Remove the linear trend, keeping the mean level: a trend inflates every
  # component's dispersion alike and can hide the period (#97).
  if (isTRUE(detrend) && n >= 5L) {
    t <- seq_len(n)
    b <- fast_lm(t, x)$coef[2]
    if (is.finite(b) && b != 0) {
      x <- x - b * (t - mean(t))
      # A second fit removes the slope's rounding, which grows with n, leaving the data's own.
      r <- x - fast_lm(t, x)$coef[2] * (t - mean(t))
      # An exactly linear series leaves only that rounding, which must not read as a period.
      if (max(abs(r - mean(r))) <= 8 * .Machine$double.eps * max(abs(variable))) x <- rep(mean(x), n)
    }
  }
  mp <- if (is.null(max.period)) NULL else as.integer(max.period)

  ans <- NNS_seas_cpp(
    variable           = x,
    modulo             = if (is.null(modulo)) NULL else as.integer(modulo),
    mod_only           = isTRUE(mod.only),
    max_period         = mp
  )

  # Optional significance: compare each period's dispersion ratio with the best
  # ratio over all candidate periods in shuffled copies of the (detrended)
  # series under an exchangeable null. Taking the best candidate of each copy
  # accounts for the full scan under that null, not arbitrary serial dependence. Periods that modulo added beyond the scan are scored on every
  # copy as well, so the null covers every candidate the observed table tests.
  if (isTRUE(p.value) && n >= 5L) {
    M <- ans$all.periods
    ratio <- function(A) A[["Coefficient.of.Variation"]] / A[["Variable.Coefficient.of.Variation"]]
    scan_max <- min(n %/% 2L, n %/% 5L)
    if (!is.null(mp)) scan_max <- max(1L, min(scan_max, mp))
    extra <- M$Period[M$Period > scan_max]
    null_best <- vapply(seq_len(as.integer(nperm)), function(b) {
      A <- if (length(extra)) {
        NNS_seas_cpp(variable = sample(x), modulo = as.integer(extra), mod_only = FALSE, max_period = mp)$all.periods
      } else {
        NNS_seas_cpp(variable = sample(x), modulo = NULL, mod_only = TRUE, max_period = mp)$all.periods
      }
      if (A$Period[1L] <= 1L) 1 else ratio(A)[1L]
    }, numeric(1L))
    obs <- ratio(M)
    pv <- vapply(obs, function(r) (1 + sum(null_best <= r)) / (1 + length(null_best)), numeric(1L))
    pv[M$Period <= 1L] <- NA_real_
    M$p.value <- pv
    sig <- M$Period[!is.na(pv) & pv <= 0.05]
    ans$all.periods <- M
    ans$periods <- if (length(sig)) as.integer(sig) else 1L
    ans$best.period <- ans$periods[1L]
  }
  
  # Plot (diagnostic)
  if (isTRUE(plot)) {
    M <- ans$all.periods
    
    # nothing to plot
    if (is.null(M) || nrow(M) == 0L) return(ans)
    
    overall_cv <- as.numeric(M$`Variable.Coefficient.of.Variation`[1L])
    metric_label <- if (all(x > 0) || all(x < 0)) "relative Gini dispersion"
                    else "Gini dispersion"
    
    # Predictive strength in [0,1] with guards
    n <- nrow(M)
    if (is.finite(overall_cv) && overall_cv > 0) {
      strength <- 1 - (M[["Coefficient.of.Variation"]] / overall_cv)
      strength <- pmin(pmax(strength, 0), 1)
      steel_pal <- grDevices::colorRampPalette(
        c("steelblue1", "steelblue2", "steelblue3", "steelblue4"))
      palette_ <- steel_pal(100L)
      idx <- pmax(1L, as.integer(round(strength * 99L + 1L)))
      point_colors <- palette_[idx]
    } else {
      point_colors <- rep("steelblue3", n)
    }
    
    # y-limits: nonnegative and wide enough to show the reference line if finite
    ymax <- if (is.finite(overall_cv) && overall_cv > 0) 2 * overall_cv else max(M[["Coefficient.of.Variation"]], 1)
    ylim <- c(0, ymax)
    
    plot(M[["Period"]], M[["Coefficient.of.Variation"]],
         xlab = "Period",
         ylab = paste("Component", metric_label),
         main = "Seasonality Candidates via Predictive Concentration\n(Lower dispersion = tighter component distributions)",
         ylim = ylim,
         col = point_colors, pch = 19)
    
    # highlight best period (table is keyed ascending by CV)
    points(M[["Period"]][1L], M[["Coefficient.of.Variation"]][1L],
           pch = 19, col = "red", cex = 1.5)
    
    # Reference CV line and centered label (only when finite)
    if (is.finite(overall_cv)) {
      abline(h = overall_cv, col = "red", lty = 2)
      usr <- graphics::par("usr")
      xmid <- mean(usr[1:2])
      graphics::text(xmid, overall_cv,
                     labels = paste0("Overall ", metric_label, "\n(Concentration threshold)"),
                     adj = c(0.5, 0.5), col = "red", xpd = NA)
    }
  }
  
  
  # Return results
  ans
}



