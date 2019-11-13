#' LPM VaR
#'
#' Generates a value at risk (VaR) based on the Lower Partial Moment ratio.
#' @param percentile numeric [0, 1]; The percentile for left-tail VaR.
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is above.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 95th percentile VaR (left-tail)
#' LPM.VaR(0.95, 0, x)
#' @export

LPM.VaR <- function(percentile, degree, x){

    sort_x <- sort(x)
    index <- findInterval(1 - percentile, LPM.ratio(0, sort_x, x))
    return(mean(sort_x[index:(index + 1)]))

}

#' UPM VaR
#'
#' Generates an upside value at risk (VaR) based on the Upper Partial Moment ratio
#' @param percentile numeric [0, 1]; The percentile for right-tail VaR.
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is below.
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 95th percentile VaR (right-tail)
#' UPM.VaR(0.95, 0, x)
#' @export

UPM.VaR <- function(percentile, degree, x){

    sort_x <- sort(x,decreasing = TRUE)
    index <- findInterval(1 - percentile, UPM.ratio(0, sort_x, x))
    return(mean(sort_x[(index - 1):index]))

}
