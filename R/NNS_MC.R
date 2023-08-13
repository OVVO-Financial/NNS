#' NNS Monte Carlo Sampling
#'
#' Monte Carlo sampling from the maximum entropy bootstrap routine \link{NNS.meboot}, ensuring the replicates are sampled from the full [-1,1] correlation space.
#'
#' @param x vector of data.
#' @param reps numeric; number of replicates to generate, \code{30} default.
#' @param lower_rho numeric \code{[-1,1]}; \code{.01} default will set the \code{from} argument in \code{seq(from, to, by)}.
#' @param upper_rho numeric \code{[-1,1]}; \code{.01} default will set the \code{to} argument in \code{seq(from, to, by)}.
#' @param by numeric; \code{.01} default will set the \code{by} argument in \code{seq(-1, 1, step)}.
#' @param exp numeric; \code{1} default will exponentially weight maximum rho value if \code{exp > 1}.  Shrinks values towards \code{upper_rho}.
#' @param type options("spearman", "pearson", "NNScor", "NNSdep"); \code{type = "spearman"}(default) dependence metric desired.
#' @param drift logical; \code{TRUE} default preserves the drift of the original series.
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
#' @references Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank Correlations in Maximum Entropy Bootstrap and Improved Monte Carlo Simulations
#' \href{https://www.ssrn.com/abstract=3621614}{https://www.ssrn.com/abstract=3621614}
#'
#' @examples
#' \dontrun{
#' # To generate a set of MC sampled time-series to AirPassengers
#' MC_samples <- NNS.MC(AirPassengers, xmin = 0)
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
                   xmin = NULL,
                   xmax = NULL, ...){


  rhos <- seq(lower_rho, upper_rho, by)
  l <- length(rhos)
  
  neg_rhos <- abs(rhos[rhos<0])
  pos_rhos <- rhos[rhos>0]
  
  exp_rhos <- rev(c((neg_rhos^exp)*-1, pos_rhos^(1/exp)))
  
  
  samples <- suppressWarnings(NNS.meboot(x = x, reps = reps, rho = exp_rhos, type = type, drift = drift,
                            xmin = xmin, xmax = xmax, ...))
  
  replicates <- samples["replicates",]
  
  rm(samples)

  ensemble <- Rfast::rowmeans(do.call(cbind, replicates))

  names(replicates) <- paste0("rho = ", exp_rhos)
  
  return(list("ensemble" = ensemble, "replicates" = replicates))
}