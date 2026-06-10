#' LPM VaR
#'
#' Generates a value at risk (VaR) quantile based on the Lower Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for left-tail VaR.
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is below.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, left-tail
#' LPM.VaR(0.05, 0, x)
#' }
#' @export

LPM.VaR <- function(percentile, degree, x) {
  x <- .NNS_prepare_VaR_x(x)
  percentile <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  if (degree == 0) {
    return(stats::quantile(x, percentile, na.rm = TRUE))
  }
  
  if (.NNS_is_supported_integer_degree(degree)) {
    return(.NNS_LPM_VaR_integer(percentile, as.integer(degree), x))
  }
  
  .NNS_LPM_VaR_optimize(percentile, degree, x)
}



#' UPM VaR
#'
#' Generates an upside value at risk (VaR) quantile based on the Upper Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for right-tail VaR.
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is above.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, right-tail
#' UPM.VaR(0.05, 0, x)
#' @export

UPM.VaR <- function(percentile, degree, x) {
  x <- .NNS_prepare_VaR_x(x)
  percentile <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  if (degree == 0) {
    return(stats::quantile(x, 1 - percentile, na.rm = TRUE))
  }
  
  if (.NNS_is_supported_integer_degree(degree)) {
    return(.NNS_UPM_VaR_integer(percentile, as.integer(degree), x))
  }
  
  .NNS_UPM_VaR_optimize(percentile, degree, x)
}



# ==============================================================================
# Internal helpers
# ==============================================================================

.NNS_prepare_VaR_x <- function(x) {
  if (inherits(x, c("tbl", "data.table"))) {
    x <- as.numeric(unlist(x))
  }
  
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  
  if (!length(x)) {
    stop("x must contain at least one non-NA numeric value.")
  }
  
  x
}



.NNS_is_supported_integer_degree <- function(degree) {
  length(degree) == 1L &&
    is.finite(degree) &&
    degree == as.integer(degree) &&
    degree >= 1 &&
    degree <= 4
}



# ------------------------------------------------------------------------------
# General integer-degree VaR inversion for degrees 1:4
# ------------------------------------------------------------------------------
#
# For integer degree d:
#
#   LPM_d(t) = sum((t - x_i)^d for x_i <= t)
#   UPM_d(t) = sum((x_i - t)^d for x_i > t)
#
# Within an interval between sorted unique observations, the below/above sets are
# fixed. Therefore LPM_d(t) and UPM_d(t) are degree-d polynomials in t.
#
# VaR solves:
#
#   LPM_d(t) / (LPM_d(t) + UPM_d(t)) = p
#
# equivalently:
#
#   (1 - p) * LPM_d(t) - p * UPM_d(t) = 0
#
# This implementation:
#   1. Sorts x once.
#   2. Builds prefix power sums P_0, P_1, ..., P_d.
#   3. Locates the correct order-statistic interval for each percentile.
#   4. Solves the exact degree-d polynomial on that interval using uniroot().
#
# This replaces the old behavior:
#   Vectorize(percentile) -> optimize() -> repeated full LPM.ratio / UPM.ratio scans.
# ------------------------------------------------------------------------------

.NNS_LPM_VaR_integer <- function(percentile, degree, x) {
  p <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  n <- length(x)
  x_sorted <- sort(x)
  x_min <- x_sorted[1L]
  x_max <- x_sorted[n]
  
  if (n == 1L || x_min == x_max) {
    return(rep(x_min, length(p)))
  }
  
  prep <- .NNS_prepare_integer_VaR_backend(x_sorted, degree)
  
  ratio_break <- .NNS_LPM_ratio_at_breaks(prep, degree)
  ratio_break <- pmin(pmax(ratio_break, 0), 1)
  
  # Protect findInterval from tiny floating-point non-monotonicity.
  ratio_break <- cummax(ratio_break)
  ratio_break[1L] <- 0
  ratio_break[length(ratio_break)] <- 1
  
  out <- numeric(length(p))
  
  left_tail <- p <= 0
  right_tail <- p >= 1
  middle <- !(left_tail | right_tail)
  
  out[left_tail] <- x_min
  out[right_tail] <- x_max
  
  if (any(middle)) {
    p_mid <- p[middle]
    
    interval <- findInterval(
      p_mid,
      ratio_break,
      rightmost.closed = TRUE
    )
    
    interval <- pmax(interval, 1L)
    interval <- pmin(interval, length(prep$unique_x) - 1L)
    
    out_mid <- numeric(length(p_mid))
    
    for (i in seq_along(p_mid)) {
      out_mid[i] <- .NNS_solve_LPM_integer_interval(
        percentile = p_mid[i],
        degree = degree,
        interval = interval[i],
        prep = prep
      )
    }
    
    out[middle] <- out_mid
  }
  
  out
}



.NNS_UPM_VaR_integer <- function(percentile, degree, x) {
  percentile <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  # UPM.ratio(t) = p is equivalent to LPM.ratio(t) = 1 - p.
  .NNS_LPM_VaR_integer(1 - percentile, degree, x)
}



.NNS_prepare_integer_VaR_backend <- function(x_sorted, degree) {
  n <- length(x_sorted)
  
  x_rle <- rle(x_sorted)
  unique_x <- x_rle$values
  k_break <- cumsum(x_rle$lengths)
  
  prefix_power <- matrix(
    0,
    nrow = degree + 1L,
    ncol = length(unique_x)
  )
  
  total_power <- numeric(degree + 1L)
  
  # Power 0 is count.
  prefix_power[1L, ] <- k_break
  total_power[1L] <- n
  
  if (degree >= 1L) {
    for (j in seq_len(degree)) {
      x_power <- x_sorted^j
      prefix_power[j + 1L, ] <- cumsum(x_power)[k_break]
      total_power[j + 1L] <- sum(x_power)
    }
  }
  
  list(
    n = n,
    unique_x = unique_x,
    k_break = k_break,
    prefix_power = prefix_power,
    total_power = total_power
  )
}



.NNS_LPM_raw_integer <- function(t, degree, prefix_power) {
  value <- 0
  
  for (j in 0:degree) {
    value <- value +
      choose(degree, j) *
      (-1)^j *
      t^(degree - j) *
      prefix_power[j + 1L]
  }
  
  value
}



.NNS_UPM_raw_integer <- function(t, degree, prefix_power, total_power) {
  suffix_power <- total_power - prefix_power
  value <- 0
  
  for (j in 0:degree) {
    value <- value +
      choose(degree, j) *
      (-1)^(degree - j) *
      t^(degree - j) *
      suffix_power[j + 1L]
  }
  
  value
}



.NNS_LPM_ratio_at_breaks <- function(prep, degree) {
  t <- prep$unique_x
  prefix_power <- prep$prefix_power
  total_power <- prep$total_power
  
  lpm <- numeric(length(t))
  upm <- numeric(length(t))
  
  for (j in 0:degree) {
    lpm <- lpm +
      choose(degree, j) *
      (-1)^j *
      t^(degree - j) *
      prefix_power[j + 1L, ]
    
    suffix_power <- total_power[j + 1L] - prefix_power[j + 1L, ]
    
    upm <- upm +
      choose(degree, j) *
      (-1)^(degree - j) *
      t^(degree - j) *
      suffix_power
  }
  
  ratio <- lpm / (lpm + upm)
  ratio[is.nan(ratio)] <- 0
  ratio
}



.NNS_pm_root_value_integer <- function(t, percentile, degree, prefix_power, total_power) {
  lpm <- .NNS_LPM_raw_integer(t, degree, prefix_power)
  upm <- .NNS_UPM_raw_integer(t, degree, prefix_power, total_power)
  
  (1 - percentile) * lpm - percentile * upm
}



.NNS_solve_LPM_integer_interval <- function(percentile, degree, interval, prep) {
  lower <- prep$unique_x[interval]
  upper <- prep$unique_x[interval + 1L]
  
  prefix_power <- prep$prefix_power[, interval]
  total_power <- prep$total_power
  
  if (lower == upper) {
    return(lower)
  }
  
  f <- function(t) {
    .NNS_pm_root_value_integer(
      t = t,
      percentile = percentile,
      degree = degree,
      prefix_power = prefix_power,
      total_power = total_power
    )
  }
  
  f_lower <- f(lower)
  f_upper <- f(upper)
  
  tol <- .Machine$double.eps^0.5
  
  if (is.finite(f_lower) && abs(f_lower) <= tol) {
    return(lower)
  }
  
  if (is.finite(f_upper) && abs(f_upper) <= tol) {
    return(upper)
  }
  
  if (
    is.finite(f_lower) &&
    is.finite(f_upper) &&
    f_lower * f_upper <= 0
  ) {
    root <- stats::uniroot(
      f,
      interval = c(lower, upper),
      tol = tol
    )$root
    
    return(pmin(pmax(root, lower), upper))
  }
  
  # Numerical safety fallback.
  # This still uses the cheap prefix-polynomial objective, not full LPM.ratio scans.
  root <- stats::optimize(
    function(t) abs(f(t)),
    interval = c(lower, upper)
  )$minimum
  
  pmin(pmax(root, lower), upper)
}



# ------------------------------------------------------------------------------
# Old-method fallback for unsupported degrees
# ------------------------------------------------------------------------------

.NNS_LPM_VaR_optimize <- function(percentile, degree, x) {
  p <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  x_min <- min(x)
  x_max <- max(x)
  
  if (x_min == x_max) {
    return(rep(x_min, length(p)))
  }
  
  vapply(
    p,
    function(pp) {
      func <- function(b) {
        abs(
          as.numeric(.Call("_NNS_LPM_ratio_RCPP", degree, b, x)) - pp
        )
      }
      
      stats::optimize(func, c(x_min, x_max))$minimum
    },
    numeric(1)
  )
}



.NNS_UPM_VaR_optimize <- function(percentile, degree, x) {
  p <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  x_min <- min(x)
  x_max <- max(x)
  
  if (x_min == x_max) {
    return(rep(x_min, length(p)))
  }
  
  vapply(
    p,
    function(pp) {
      func <- function(b) {
        abs(
          as.numeric(.Call("_NNS_UPM_ratio_RCPP", degree, b, x)) - pp
        )
      }
      
      stats::optimize(func, c(x_min, x_max))$minimum
    },
    numeric(1)
  )
}