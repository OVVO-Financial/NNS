#' LPM VaR
#'
#' Generates a value at risk (VaR) quantile based on the Lower Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for left-tail VaR (vectorized).
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is below.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5% quantile, left-tail
#' LPM.VaR(0.05, 0, x)
#' @export

LPM.VaR <- function(percentile, degree, x){

    if(any(class(x)=="tbl")) x <- as.vector(unlist(x))

    percentile <- pmax(pmin(percentile, 1), 0)

    l <- length(x)

    if(degree == 0){
        td <- tdigest::tdigest(x, compression = max(100, log(l,10)*100))
        q <- tryCatch(tdigest::tquantile(td, percentile),
                      error = quantile(x, percentile, na.rm = TRUE))
        return(q)
    } else {
        func <- function(b){
            abs(LPM.ratio(degree, b, x) - percentile)
        }
        return(optimize(func, c(min(x),max(x)))$minimum)
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
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5% quantile, right-tail
#' UPM.VaR(0.05, 0, x)
#' @export

UPM.VaR <- function(percentile, degree, x){

    if(any(class(x)=="tbl")) x <- as.vector(unlist(x))

    percentile <- pmax(pmin(percentile, 1), 0)

    l <- length(x)
    if(degree==0){
        td <- tdigest::tdigest(x, compression = max(100, log(l,10)*100))
        q <- tryCatch(tdigest::tquantile(td, 1 - percentile),
                      error = quantile(x, 1 - percentile, na.rm = TRUE))
        return(q)
    } else {
        func <- function(b){
            abs(LPM.ratio(degree, b, x) - (1 - percentile))
        }
        return(optimize(func, c(min(x),max(x)))$minimum)
    }

}

UPM.VaR <- Vectorize(UPM.VaR, vectorize.args = "percentile")

