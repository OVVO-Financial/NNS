#' NNS Causation
#'
#' Returns the directional lagged-association diagnostic between two variables from observational data, as an index in [0, 1] for each direction (0.5 = no net direction).
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param factor.2.dummy logical; \code{FALSE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.  Includes dependent variable \code{y}.
#' @param tau options: ("cs", "ts", integer); 0 (default) Lag between cause and effect.  \code{(tau = [integer])} compares the lagged pairs \code{(x[t - tau], y[t])} for x ---> y and \code{(y[t - tau], x[t])} for y ---> x.  \code{(tau = 0)} or \code{(tau = "cs")} compares the contemporaneous pairs (cross-sectional data).  \code{(tau = "ts")} uses one lag for both directions, the mean of the first \link{NNS.seas} period at most \code{sqrt(length(x))} of each series.
#' @param plot logical; \code{FALSE} (default) Plots both lagged pairs on the copula scale with their indices.
#' @param p.value logical; \code{FALSE} (default) If \code{TRUE}, runs a permutation test of the signed statistic \code{C(x--->y) - C(y--->x)}.
#' @param nperm integer; number of permutations to use when \code{p.value = TRUE}. Default 100.
#' @param permute one of "auto" (default), "y", "x", "both" or "shift"; how the null distribution is built.  \code{"shift"} circularly shifts \code{y} by random distinct offsets at least \code{max(tau + 1, n / 10)} steps from alignment, so the within-series ordering is retained, up to the wraparound boundary, while cross-series alignment changes. This does not condition on the effect's own past.  When at most \code{nperm} offsets are admissible, each is used once and the reported \code{nperm} is their number.  \code{"y"}, \code{"x"} and \code{"both"} shuffle the values, an independence null for cross-sections.  \code{"auto"} uses \code{"shift"} when the lag is positive (including \code{tau = "ts"}) and \code{"y"} otherwise.  The lag of \code{tau = "ts"} is resolved once on the observed series and kept for every permutation.  For a matrix \code{x} each choice applies to every pair, with column \code{i} as \code{x} and column \code{j} as \code{y}; \code{"auto"} shuffles every column for cross-sections.
#' @param seed optional integer seed for reproducibility of the permutation test.
#' @param conf.int numeric; 0.95 (default) confidence level for the partial-moment based interval computed on the permutation null distribution.
#'
#' @return If \code{p.value = FALSE} returns a named vector \code{c("C(x--->y)", "C(y--->x)")} of directional association indices in [0, 1] that sum to 1: above 0.5 means the first lagged direction has the larger score, and 0.5 is no net direction under this diagnostic.  If \code{p.value = TRUE} returns a list with components:
#'  * \code{causation}: the vector above.
#'  * \code{p.value}: a list with empirical two-sided and one-sided p-values (x_causes_y, y_causes_x) of the signed statistic \code{C(x--->y) - C(y--->x)}, its null distribution and observed value, the partial-moment interval of the null, and metadata (permute, nperm).
#' For a matrix \code{x}, entry \code{[i, j]} is \code{C(column i ---> column j)}, with a 0 diagonal.  If \code{p.value = TRUE} for a matrix, the function returns a list with components:
#'   * \code{causality}: the causation matrix.
#'   * \code{lower_CI}, \code{upper_CI}: matrices with the partial-moment interval of each index under the permutation null.
#'   * \code{p.value_matrix}: matrix of empirical two-sided p-values of \code{C(i ---> j) - C(j ---> i)}.
#'
#' @details
#' Following book chapter 14, each direction is evaluated on explicitly lagged pairs mapped to the copula scale with per-column ranks (the degree-0 lower partial moment), so the result is invariant to shifts, rescaling and any increasing transform of either variable (reversing a variable's order changes the copula, and so the index).  For x ---> y the raw score is
#' 
#' \deqn{\tilde{C}_{x \to y} = \frac{1}{2}\left[P_{x \to y} + \max(\rho_{x \to y} - \rho_{y \to x}, 0)\right],}
#' 
#' where \eqn{P_{x \to y}} is the share of the deviation of \code{y[t]} resolved by conditioning on the position of \code{x[t - tau]}, the deviation being the degree-2 lower plus upper partial moment about the mean with the \eqn{m/(m-1)} factor, used as a finite-sample normalization under the permutation reference, and \eqn{\rho_{x \to y}} is the debiased directional dependence of \code{y[t]} on \code{x[t - tau]} from \link{NNS.dep}'s cells (the lagged co-partial moments of the pair), averaged over the levels of the partition on \code{x[t - tau]}.  \link{NNS.dep} reports the strongest level; the average avoids that maximization step. The neutral value is 0.5, but finite-sample behavior depends on the data and the null construction.  The index is
#' 
#' \deqn{C_{x \to y} = \frac{\tilde{C}_{x \to y} + \kappa}{\tilde{C}_{x \to y} + \tilde{C}_{y \to x} + 2\kappa}, \quad C_{y \to x} = 1 - C_{x \to y},}
#' 
#' with a neutral pseudo-count \eqn{\kappa = 1/\sqrt{n - \tau}}, a regularization that shrinks weak directional evidence toward 0.5. It does not make the index's spread under independence constant: the spread grows with persistence. For independent pairs with n = 300 and tau = 20, the standard deviation of the index was about 0.08 for iid series, 0.14 for AR(1) series with coefficient 0.9 and 0.17 for random walks, and more than half of the persistent pairs fell more than 0.1 from 0.5. Distance from 0.5 alone is therefore not a significance criterion.  An invertible relationship (e.g. \code{y = 2 * x + 1}) has no cross-sectional direction and returns 0.5; a lag can produce unequal directional association in that case.  The index measures directional lagged association. It does not adjust for the response's own history, confounding variables, common drivers, or interventions. A common driver that reaches \code{y} later than \code{x} can favor x ---> y without a causal effect. Circular-shift inference requires a null model invariant under the allowed shifts; stationarity and preserved autocorrelation alone do not establish this. Restricted offsets and boundary wraparound make this a finite resampling diagnostic, not an exact conditional-causality test. Shuffled-value inference instead requires exchangeability. The reported interval is an interval of the permutation null, not an uncertainty interval for an identified causal parameter.
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
                     permute = c("auto", "y", "x", "both", "shift"),
                     seed = NULL,
                     conf.int = 0.95){
  permute <- match.arg(permute)
  if(!is.null(seed)) set.seed(seed)

  # Unwrap single-column data frames here as well as in the core: the
  # permutation loop below calls sample() on x and y, which on a data frame
  # would shuffle its one column as a single element and leave the rows in
  # place, so every null value would equal the observed statistic.
  if(is.data.frame(x) && ncol(x) == 1) x <- x[[1L]]
  if(is.data.frame(y) && ncol(y) == 1) y <- y[[1L]]
  
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

  # One numeric lag for the observed statistic and every permutation (#96):
  # tau = "ts" is resolved once on the observed series, so a permuted series
  # without a usable period can't fail, and the null compares like with like.
  xv <- as.numeric(x); yv <- as.numeric(y)
  tau_num <- if (identical(tau, "ts")) .nns_caus_ts_lag(xv, yv) else if (identical(tau, "cs")) 0L else as.integer(tau)
  # Time series keep their own dynamics under a circular shift; cross-sections
  # are shuffled.
  if (permute == "auto") permute <- if (tau_num > 0L) "shift" else "y"
  if (permute == "shift") {
    offsets <- .nns_caus_shift_draw(length(yv), tau_num, nperm)
    nperm <- length(offsets)
  }
  
  # Build null distribution via permutations
  null_vals <- numeric(nperm)
  for(b in seq_len(nperm)){
    x_perm <- xv; y_perm <- yv
    if(permute == "y") y_perm <- sample(yv, length(yv), replace = FALSE)
    if(permute == "x") x_perm <- sample(xv, length(xv), replace = FALSE)
    if(permute == "both"){
      x_perm <- sample(xv, length(xv), replace = FALSE)
      y_perm <- sample(yv, length(yv), replace = FALSE)
    }
    if(permute == "shift") y_perm <- .nns_caus_shift(yv, offsets[b])
    cp_perm <- NNS.caus_core(x = x_perm, y = y_perm,
                             factor.2.dummy = FALSE,
                             tau = tau_num,
                             plot = FALSE)
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
