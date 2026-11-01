#' NNS Causation
#'
#' Returns the directional causation between two variables from observational data, as an index in [0, 1] for each direction (0.5 = no net direction).
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param factor.2.dummy logical; \code{FALSE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.  Includes dependent variable \code{y}.
#' @param tau options: ("cs", "ts", integer); 0 (default) Lag between cause and effect.  \code{(tau = [integer])} compares the lagged pairs \code{(x[t - tau], y[t])} for x ---> y and \code{(y[t - tau], x[t])} for y ---> x.  \code{(tau = 0)} or \code{(tau = "cs")} compares the contemporaneous pairs (cross-sectional data).  \code{(tau = "ts")} uses one lag for both directions, the mean of the first \link{NNS.seas} period at most \code{sqrt(length(x))} of each series.
#' @param plot logical; \code{FALSE} (default) Plots both lagged pairs on the copula scale with their indices.
#' @param p.value logical; \code{FALSE} (default) If \code{TRUE}, runs a permutation test of the signed statistic \code{C(x--->y) - C(y--->x)}.
#' @param nperm integer; number of permutations to use when \code{p.value = TRUE}. Default 100.
#' @param permute one of "both", "y", or "x"; which variable(s) to shuffle when constructing the null distribution.
#' @param seed optional integer seed for reproducibility of the permutation test.
#' @param conf.int numeric; 0.95 (default) confidence level for the partial-moment based interval computed on the permutation null distribution.
#'
#' @return If \code{p.value = FALSE} returns a named vector \code{c("C(x--->y)", "C(y--->x)")} of directional causation indices in [0, 1] that sum to 1: above 0.5 is evidence that the first direction leads, 0.5 is no net direction.  If \code{p.value = TRUE} returns a list with components:
#'  * \code{causation}: the vector above.
#'  * \code{p.value}: a list with empirical two-sided and one-sided p-values (x_causes_y, y_causes_x) of the signed statistic \code{C(x--->y) - C(y--->x)}, its null distribution and observed value, the partial-moment interval of the null, and metadata (permute, nperm).
#' For a matrix \code{x}, entry \code{[i, j]} is \code{C(column i ---> column j)}, with a 0 diagonal.  If \code{p.value = TRUE} for a matrix, the function returns a list with components:
#'   * \code{causality}: the causation matrix.
#'   * \code{lower_CI}, \code{upper_CI}: matrices with the partial-moment interval of each index under the permutation null.
#'   * \code{p.value_matrix}: matrix of empirical two-sided p-values of \code{C(i ---> j) - C(j ---> i)}.
#'
#' @details
#' Following book chapter 14, each direction is evaluated on explicitly lagged pairs mapped to the copula scale with per-column ranks (the degree-0 lower partial moment), so the result is invariant to shifts, rescaling and any monotone transform of either variable.  For x ---> y the raw score is
#' 
#' \deqn{\tilde{C}_{x \to y} = \frac{1}{2}\left[P_{x \to y} + \max(\rho_{x \to y} - \rho_{y \to x}, 0)\right],}
#' 
#' where \eqn{P_{x \to y}} is the share of the degree-1 deviation (lower plus upper partial moment about the mean) of \code{y[t]} resolved by conditioning on the position of \code{x[t - tau]}, and \eqn{\rho_{x \to y}} is \code{NNS.dep(x[t - tau], y[t], asym = TRUE)$Dependence}, the directional dependence built from the lagged co-partial moments of the pair.  The index is
#' 
#' \deqn{C_{x \to y} = \frac{\tilde{C}_{x \to y} + \kappa}{\tilde{C}_{x \to y} + \tilde{C}_{y \to x} + 2\kappa}, \quad C_{y \to x} = 1 - C_{x \to y},}
#' 
#' with a neutral pseudo-count \eqn{\kappa = 1/\sqrt{n - \tau}}, the order of the raw scores' sampling variation, so that weak or absent directional evidence stays near 0.5.  An invertible relationship (e.g. \code{y = 2 * x + 1}) has no cross-sectional direction and returns 0.5; the lag is what identifies direction in that case.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#'
#' \dontrun{
#' ## x causes y (y is a function of x, x is not a function of y)
#' set.seed(123)
#' x <- rnorm(1000) ; y <- x ^ 2
#' NNS.caus(x, y, tau = "cs")
#'
#' ## Lagged causation: y follows x by one period
#' x <- rnorm(500) ; y <- c(0, 0.8 * x[-500]) + 0.3 * rnorm(500)
#' NNS.caus(x, y, tau = 1)
#'
#' ## Causal matrix: entry [i, j] is C(column i ---> column j)
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
