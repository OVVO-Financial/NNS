#' NNS Stochastic Superiority
#'
#' Computes stochastic superiority between two numeric vectors as the empirical
#' probability that an observation from \code{x} exceeds an observation from
#' \code{y}, with optional tie adjustment and optional confidence intervals via
#' maximum entropy bootstrap.
#'
#' \code{NNS.SS} returns:
#' \deqn{P(X > Y),}
#' the tie probability
#' \deqn{P(X = Y),}
#' and the tie-adjusted stochastic superiority measure
#' \deqn{P^* = P(X > Y) + \frac{1}{2} P(X = Y).}
#'
#' When \code{confidence.interval = TRUE}, confidence bounds for \code{P^*}
#' are computed from \code{\link{NNS.meboot}} bootstrap replicates using
#' \code{\link{LPM.VaR}} and \code{\link{UPM.VaR}} with \code{degree = 0}.
#'
#' @usage
#' NNS.SS(
#'   x,
#'   y,
#'   confidence.interval = FALSE,
#'   reps = 999,
#'   ci = 0.95,
#'   rho = 1
#' )
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param confidence.interval logical; \code{FALSE} (default) returns only the
#' empirical stochastic superiority measures. Set to \code{TRUE} to compute
#' bootstrap confidence intervals for \code{p_star}.
#' @param reps numeric; number of maximum entropy bootstrap replicates used when
#' \code{confidence.interval = TRUE}. Default is \code{999}.
#' @param ci numeric in \eqn{(0, 1)}; confidence level used for the bootstrap
#' interval when \code{confidence.interval = TRUE}. Default is \code{0.95}.
#' @param rho numeric; dependence target passed to \code{\link{NNS.meboot}}.
#' Default is \code{1}.
#'
#' @details
#' Missing values are removed from both \code{x} and \code{y} using
#' \code{stats::na.omit}. The empirical estimates are computed via a fast sorted
#' comparison routine rather than explicit pairwise expansion of all
#' \code{x}-\code{y} combinations.
#'
#' For continuous data, \code{p_tie} will typically be zero, so \code{p_star}
#' and \code{p_gt} will be identical up to numerical precision. For discrete
#' data, \code{p_star} provides the standard tie-adjusted superiority measure.
#'
#' When \code{confidence.interval = TRUE}, the interval is constructed from the
#' empirical bootstrap distribution of \code{p_star}, where
#' \eqn{\alpha = 1 - ci}. The lower bound is obtained from
#' \code{\link{LPM.VaR}} evaluated at \eqn{\alpha / 2}, and the upper bound is
#' obtained from \code{\link{UPM.VaR}} evaluated at \eqn{\alpha / 2}, both with
#' \code{degree = 0}.
#'
#' @return
#' If \code{confidence.interval = FALSE}, returns a list containing:
#' \describe{
#'   \item{\code{p_gt}}{empirical probability that \code{x > y}.}
#'   \item{\code{p_tie}}{empirical probability that \code{x = y}.}
#'   \item{\code{p_star}}{tie-adjusted stochastic superiority probability.}
#' }
#'
#' If \code{confidence.interval = TRUE}, returns a list containing:
#' \describe{
#'   \item{\code{p_gt}}{empirical probability that \code{x > y}.}
#'   \item{\code{p_tie}}{empirical probability that \code{x = y}.}
#'   \item{\code{p_star}}{tie-adjusted stochastic superiority probability.}
#'   \item{\code{lower}}{lower confidence bound for \code{p_star}.}
#'   \item{\code{upper}}{upper confidence bound for \code{p_star}.}
#'   \item{\code{ci}}{confidence level used.}
#'   \item{\code{reps}}{number of bootstrap replicates used.}
#'   \item{\code{boot_vals}}{bootstrap replicate values of \code{p_star}.}
#' }
#'
#' @note
#' This function measures stochastic superiority as a pairwise exceedance
#' probability. This is distinct from first-, second-, or third-degree
#' stochastic dominance; see \code{\link{NNS.FSD}}, \code{\link{NNS.SSD}}, and
#' \code{\link{NNS.TSD}} for dominance testing.
#'
#' @author
#' Fred Viole, OVVO Financial Systems
#'
#' @references
#' \itemize{
#'   \item Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank
#'   Correlations in Maximum Entropy Bootstrap and Improved Monte Carlo
#'   Simulations. \doi{10.2139/ssrn.3621614}
#'   \item Viole, F. and Nawrocki, D. (2013)
#'   \emph{Nonlinear Nonparametric Statistics: Using Partial Moments}.
#'   ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(200, mean = 0.4, sd = 1)
#' y <- rnorm(200, mean = 0.0, sd = 1)
#'
#' # Empirical stochastic superiority
#' NNS.SS(x, y)
#'
#' # With confidence intervals
#' NNS.SS(x, y, confidence.interval = TRUE, reps = 999, ci = 0.95)
#'
#' # Discrete example with ties
#' x <- sample(1:5, 100, replace = TRUE)
#' y <- sample(1:5, 100, replace = TRUE)
#' NNS.SS(x, y)
#' }
#'
#' @export


NNS.SS <- function(x,
                   y,
                   confidence.interval = FALSE,
                   reps = 999,
                   ci = 0.95,
                   rho = 1) {
  
  x <- as.numeric(stats::na.omit(x))
  y <- as.numeric(stats::na.omit(y))
  
  if (length(x) == 0L || length(y) == 0L) {
    stop("x and y must both contain at least one non-missing value.")
  }
  
  if (!is.logical(confidence.interval) || length(confidence.interval) != 1L || is.na(confidence.interval)) {
    stop("confidence.interval must be a single TRUE/FALSE value.")
  }
  
  if (!confidence.interval) {
    return(stoch_superiority_cpp(x = x, y = y))
  }
  
  if (!is.numeric(reps) || length(reps) != 1L || reps < 2) {
    stop("reps must be a single number >= 2.")
  }
  
  if (!is.numeric(ci) || length(ci) != 1L || ci <= 0 || ci >= 1) {
    stop("ci must be a single number in (0, 1).")
  }
  
  empirical <- stoch_superiority_cpp(x = x, y = y)
  
  # NNS maximum-entropy bootstrap replicates
  x_boots <- NNS::NNS.meboot(x = x, reps = reps, rho = rho)["replicates", ]$replicates
  y_boots <- NNS::NNS.meboot(x = y, reps = reps, rho = rho)["replicates", ]$replicates
  
  boot_vals <- vapply(seq_len(reps), function(i) {
    stoch_superiority_cpp(
      x = x_boots[, i],
      y = y_boots[, i]
    )[["p_star"]]
  }, numeric(1))
  
  alpha <- (1 - ci) / 2
  
  list(
    p_gt      = empirical$p_gt,
    p_tie     = empirical$p_tie,
    p_star    = empirical$p_star,
    lower     = as.numeric(NNS::LPM.VaR(alpha, degree = 0, x = boot_vals)),
    upper     = as.numeric(NNS::UPM.VaR(alpha, degree = 0, x = boot_vals)),
    ci        = ci,
    reps      = reps,
    boot_vals = boot_vals
  )
}