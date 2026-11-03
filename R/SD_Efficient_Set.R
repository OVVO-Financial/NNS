#' NNS SD Efficient Set
#'
#' Determines the set of stochastic dominant variables for various degrees.
#'
#' @param x a numeric matrix or data frame.
#' @param degree numeric options: (1, 2, 3); Degree of stochastic dominance test from (1, 2 or 3).
#' @param type \code{"discrete"} (default) uses empirical FSD. At degree 1, \code{"continuous"} retains the degree-one area-ratio (Omega/expectile) ordering, which need not imply FSD or preference for every increasing utility. SSD and TSD are unaffected by \code{type}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @return Returns set of stochastic dominant variable names.
#' @note The curves are compared on the merged sample points and, for TSD and continuous area order, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous area order) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \doi{10.2139/ssrn.3002675}
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y<-rnorm(100) ; z<-rnorm(100)
#' A <- cbind(x, y, z)
#' NNS.SD.efficient.set(A, 1)
#' }
#' @export



NNS.SD.efficient.set <- function(x, degree, type = "discrete", status = TRUE) {
  .Call(`_NNS_NNS_SD_efficient_set_parallel_cpp`,
        as.matrix(x), as.integer(degree), as.character(type), as.logical(status))
}
