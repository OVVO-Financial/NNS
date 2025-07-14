#' LPM VaR
#'
#' Generates a value at risk (VaR) quantile based on the Lower Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for left-tail VaR (vectorized).
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is below.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
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
  if (inherits(x, c("tbl","data.table"))) x <- as.numeric(unlist(x))
  percentile <- pmin(pmax(percentile, 0), 1)
  
  if (degree == 0) {
    return(quantile(x, percentile, na.rm = TRUE))
  }
  
  # here we call the C++ ratio directly:
  func <- function(b) {
    abs(
       as.numeric(.Call("_NNS_LPM_ratio_RCPP", degree, b, x)) - percentile
    )
  }
  
  if (min(x) != max(x)) {
    optimize(func, c(min(x), max(x)))$minimum
  } else {
    min(x)
  }
}
LPM.VaR <- Vectorize(LPM.VaR, vectorize.args = "percentile")



#' UPM VaR
#'
#' Generates an upside value at risk (VaR) quantile based on the Upper Partial Moment ratio
#' @param percentile numeric [0, 1]; The percentile for right-tail VaR (vectorized).
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is above.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, right-tail
#' UPM.VaR(0.05, 0, x)
#' @export

UPM.VaR <- function(percentile, degree, x) {
  if (inherits(x, c("tbl","data.table"))) x <- as.numeric(unlist(x))
  percentile <- pmin(pmax(percentile, 0), 1)
  
  if (degree == 0) {
    return(quantile(x, 1 - percentile, na.rm = TRUE))
  }

  func <- function(b) {
    abs(
      as.numeric(.Call("_NNS_UPM_ratio_RCPP", degree, b, x)) - percentile
    )
  }
  
  if (min(x) != max(x)) {
    optimize(func, c(min(x), max(x)))$minimum
  } else {
    min(x)
  }
}
UPM.VaR <- Vectorize(UPM.VaR, vectorize.args = "percentile")

