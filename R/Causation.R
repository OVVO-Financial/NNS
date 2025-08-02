#' NNS Causation
#'
#' Returns the causality from observational data between two variables.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param factor.2.dummy logical; \code{FALSE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.  Includes dependent variable \code{y}.
#' @param tau options: ("cs", "ts", integer); 0 (default) Number of lagged observations to consider (for time series data).  Otherwise, set \code{(tau = "cs")} for cross-sectional data.  \code{(tau = "ts")} automatically selects the lag of the time series data, while \code{(tau = [integer])} specifies a time series lag.
#' @param plot logical; \code{FALSE} (default) Plots the raw variables, tau normalized, and cross-normalized variables.
#' @param p.value logical; \code{FALSE} (default) If \code{TRUE}, runs a permutation test to compute empirical p-values for the signed causation from x -> y.
#' @param nperm integer; number of permutations to use when \code{p.value = TRUE}. Default 100.
#' @param permute one of "both", "y", or "x"; which variable(s) to shuffle when constructing the null distribution.
#' @param seed optional integer seed for reproducibility of the permutation test.
#' @param conf.int numeric; 0.95 (default) confidence level for the partial-moment based interval computed on the permutation null distribution.
#'
#' @return If \code{p.value=FALSE} returns the original causation vector of length 3 (directional given/received and net), named either "C(x--->y)" or "C(y--->x)" in the third slot.  If \code{p.value=TRUE} returns a list with components:
#'  * \code{causation}: the original causation vector as above.
#'  * \code{p.value}: a list with empirical two-sided and one-sided p-values (x_causes_y, y_causes_x), the null distribution, the observed signed statistic, and metadata (permute, nperm).
#' If \code{p.value=TRUE} for a matrix, the function returns a list with components:
#'   * \code{causality}: the causality matrix.
#'   * \code{lower_CI}: matrix of lower confidence bounds (partial-moment based).
#'   * \code{upper_CI}: matrix of upper confidence bounds (partial-moment based).
#'   * \code{p.value}: matrix of empirical two-sided p-values.

#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#'
#' \dontrun{
#' ## x causes y...
#' set.seed(123)
#' x <- rnorm(1000) ; y <- x ^ 2
#' NNS.caus(x, y, tau = "cs")
#'
#' ## Causal matrix without per factor causation
#' NNS.caus(iris, tau = 0)
#'
#' ## Causal matrix with per factor causation
#' NNS.caus(iris, factor.2.dummy = TRUE, tau = 0)
#' }
#' @export


NNS.caus <- function(x, y = NULL,
                     factor.2.dummy = FALSE,
                     tau = 0,
                     plot = FALSE,
                     p.value = FALSE,
                     nperm = 100L,
                     permute = c("y", "x", "both"),
                     seed = NULL,
                     conf.int = 0.95){
  permute <- match.arg(permute)
  if(!is.null(seed)) set.seed(seed)
  
  # Base causation (delegates to core)
  cp <- NNS.caus_core(x = x, y = y,
                      factor.2.dummy = factor.2.dummy,
                      tau = tau,
                      plot = plot,
                      p.value = p.value,
                      nperm = nperm,
                      permute = permute,
                      seed = seed,
                      conf.int = conf.int)
  
  if (is.null(y)) return(cp)
  if (!isTRUE(p.value)) return(cp)
  
  # Compute observed signed statistic x -> y
  T_obs <- signed_from_cp(cp)
  
  # Build null distribution via permutations
  null_vals <- numeric(nperm)
  for(b in seq_len(nperm)){
    if(permute == "y"){
      y_perm <- sample(y, length(y), replace = FALSE)
      cp_perm <- NNS.caus_core(x = x, y = y_perm,
                               factor.2.dummy = factor.2.dummy,
                               tau = tau,
                               plot = FALSE)
    } else if(permute == "x"){
      x_perm <- sample(x, length(x), replace = FALSE)
      cp_perm <- NNS.caus_core(x = x_perm, y = y,
                               factor.2.dummy = factor.2.dummy,
                               tau = tau,
                               plot = FALSE)
    } else if(permute == "both"){
      x_perm <- sample(x, length(x), replace = FALSE)
      y_perm <- sample(y, length(y), replace = FALSE)
      cp_perm <- NNS.caus_core(x = x_perm, y = y_perm,
                               factor.2.dummy = factor.2.dummy,
                               tau = tau,
                               plot = FALSE)
    }
    null_vals[b] <- signed_from_cp(cp_perm)
  }
  
  # Empirical p-values (with +1 correction)
  p_two_sided <- (1 + sum(abs(null_vals) >= abs(T_obs))) / (1 + nperm)
  p_x_causes_y <- (1 + sum(null_vals >= T_obs)) / (1 + nperm)
  p_y_causes_x <- (1 + sum(null_vals <= T_obs)) / (1 + nperm)
  
  result <- list(
    causation = cp,
    p.value = list(
      two.sided = p_two_sided,
      x_causes_y = p_x_causes_y,
      y_causes_x = p_y_causes_x,
      null_distribution = null_vals,
      observed_signed = T_obs,
      permute = permute,
      nperm = nperm,
      lower_CI = LPM.VaR((1 - conf.int)/2, 0, null_vals),
      upper_CI = UPM.VaR((1 - conf.int)/2, 0, null_vals)
    )
  )
  # ensure no NAs in the permutation outputs
  result$p.value$two.sided <- ifelse(is.na(result$p.value$two.sided), 0, result$p.value$two.sided)
  result$p.value$x_causes_y <- ifelse(is.na(result$p.value$x_causes_y), 0, result$p.value$x_causes_y)
  result$p.value$y_causes_x <- ifelse(is.na(result$p.value$y_causes_x), 0, result$p.value$y_causes_x)
  result$p.value$null_distribution[is.na(result$p.value$null_distribution)] <- 0
  result$p.value$lower_CI <- ifelse(is.na(result$p.value$lower_CI), 0, result$p.value$lower_CI)
  result$p.value$upper_CI <- ifelse(is.na(result$p.value$upper_CI), 0, result$p.value$upper_CI)

  return(result)
}
