#' NNS Monte Carlo Sampling
#'
#' Monte Carlo sampling from the maximum entropy bootstrap routine \link{NNS.meboot}, ensuring the replicates are sampled across the full [-1,1] range of \code{rho} (the target dependence of the detrended residuals, see \link{NNS.meboot}).
#'
#' @param x vector of data.
#' @param reps numeric; number of replicates to generate for each \code{rho} value, \code{30} default.  The defaults give 201 \code{rho} values and 6030 replicates.
#' @param lower_rho numeric \code{[-1,1]}; \code{-1} default will set the \code{from} argument in \code{seq(from, to, by)}.
#' @param upper_rho numeric \code{[-1,1]}; \code{1} default will set the \code{to} argument in \code{seq(from, to, by)}.
#' @param by numeric; \code{.01} default will set the \code{by} argument in \code{seq(-1, 1, step)}.
#' @param exp positive finite numeric; \code{1} default will exponentially weight maximum rho value if \code{exp > 1}.  Shrinks values towards \code{upper_rho}, keeping them within \code{[lower_rho, upper_rho]}.
#' @param type options("spearman", "pearson", "NNScor", "NNSdep"); \code{type = "spearman"}(default) dependence metric desired.
#' @param drift logical; \code{drift = TRUE} (default) preserves the drift of the original series.
#' @param target_drift numerical; \code{target_drift = NULL} (default) Specifies the desired drift when \code{drift = TRUE}, i.e. a risk-free rate of return.
#' @param target_drift_scale numerical; instead of calculating a \code{target_drift}, provide a scalar to the existing drift when \code{drift = TRUE}.
#' @param xmin numeric; the lower limit for the left tail.
#' @param xmax numeric; the upper limit for the right tail.
#' @param ... possible additional arguments to be passed to \link{NNS.meboot}.
#'
#' @return
#' \itemize{
#'   \item{ensemble} average observation over all replicates as a vector.
#'   \item{replicates} maximum entropy bootstrap replicates as a list for each \code{rho}.
#' }
#'
#' @references Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank Correlations in Maximum Entropy Bootstrap and Improved Monte Carlo Simulations.  \doi{10.2139/ssrn.3621614}
#'
#' @examples
#' \dontrun{
#' # To generate a set of MC sampled time-series to AirPassengers
#' MC_samples <- NNS.MC(AirPassengers, reps = 10, lower_rho = -1, upper_rho = 1, by = .5, xmin = 0)
#' }
#' @export


NNS.MC <- function(x,
                   reps = 30,
                   lower_rho = -1,
                   upper_rho = 1,
                   by = .01,
                   exp = 1,
                   type = "spearman",
                   drift = TRUE,
                   target_drift = NULL,
                   target_drift_scale = NULL,
                   xmin = NULL,
                   xmax = NULL, ...){


  scalar_finite <- function(z) is.numeric(z) && length(z) == 1L && is.finite(z)
  if (!scalar_finite(exp) || exp <= 0) stop("exp must be a single finite number > 0.")
  if (!scalar_finite(lower_rho) || !scalar_finite(upper_rho) ||
      lower_rho < -1 || upper_rho > 1 || lower_rho > upper_rho)
    stop("lower_rho and upper_rho must satisfy -1 <= lower_rho <= upper_rho <= 1.")
  if (!scalar_finite(by) || by <= 0) stop("by must be a single finite number > 0.")

  rhos <- seq(lower_rho, upper_rho, by)
  l <- length(rhos)
  
  # Warp each side of 0 (or of the bound nearest 0) within its own span, so the grid stays inside the bounds.
  if(exp != 1){
    pivot <- min(max(0, lower_rho), upper_rho)
    low <- rhos <= pivot
    if(pivot > lower_rho){
      w <- ((pivot - rhos[low]) / (pivot - lower_rho))^exp
      rhos[low] <- (1 - w) * pivot + w * lower_rho
    }
    if(upper_rho > pivot){
      w <- ((rhos[!low] - pivot) / (upper_rho - pivot))^(1/exp)
      rhos[!low] <- (1 - w) * pivot + w * upper_rho
    }
  }
  
  exp_rhos <- rev(rhos)

  if(is.null(target_drift)){
    if(!is.null(target_drift_scale)){
      replicates <- NNS.meboot(x = x, reps = reps, rho = exp_rhos, type = type, drift = TRUE,
                               target_drift_scale = target_drift_scale, 
                               xmin = xmin, xmax = xmax, ...)["replicates",]
    } else {
      replicates <- NNS.meboot(x = x, reps = reps, rho = exp_rhos, type = type, drift = drift,
                               xmin = xmin, xmax = xmax, ...)["replicates",]
    } 
  } else {
    replicates <- NNS.meboot(x = x, reps = reps, rho = exp_rhos, type = type, drift = TRUE,
                             target_drift = target_drift,
                             xmin = xmin, xmax = xmax, ...)["replicates",]
  }


  ensemble <- rowMeans(do.call(cbind, replicates))

  names(replicates) <- paste0("rho = ", exp_rhos)
  
  return(.NNS.out(list("ensemble" = ensemble, "replicates" = replicates)))
}