#' LPM VaR
#'
#' Generates a VaR based on the Lower Partial Moment ratio
#' @param percentile numeric [0,1]; The percentile for left-tail VaR.
#' @param degree integer; \code{(degree=0)} for discrete distributions, \code{(degree=1)} for continuous distributions.
#' @param x a numeric vector.
#' @param extend options: ("yes",NULL); \code{NULL} (default) Sets the \code{"extendInt"} argument from \link{uniroot}.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is above.
#' @note If endpoint error is generated, set \code{(extend="yes")}.
#' @keywords VaR
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' ## For 95th percentile VaR (left-tail)
#' LPM.VaR(0.95,0,x)
#' @export

LPM.VaR <- function(percentile,degree,x,extend=NULL){

  f<- function(tgt) (sum((tgt - (x[x <= tgt]))^degree)/length(x))/((sum((tgt - (x[x <= tgt]))^degree)/length(x))+ sum(((x[x > tgt]) - tgt)^degree)/length(x))- (1-percentile)

  return(uniroot(f,lower=min(x),upper = max(x),extendInt = extend)$root)

}

#' UPM VaR
#'
#' Generates an upside VaR based on the Upper Partial Moment ratio
#' @param percentile numeric [0,1]; The percentile for right-tail VaR.
#' @param degree integer; \code{(degree=0)} for discrete distributions, \code{(degree=1)} for continuous distributions.
#' @param x a numeric vector.
#' @param extend options: ("yes",NULL); \code{NULL} (default) Sets the \code{"extendInt"} argument from \link{uniroot}.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is below.
#' @keywords VaR
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' ## For 95th percentile VaR (right-tail)
#' UPM.VaR(0.95,0,x)
#' @export

UPM.VaR <- function(percentile,degree,x,extend=NULL){


  f<- function(tgt) (sum(((x[x > tgt]) - tgt)^degree)/length(x))/((sum((tgt - (x[x <= tgt]))^degree)/length(x))+ sum(((x[x > tgt]) - tgt)^degree)/length(x))- (1-percentile)

  return(uniroot(f,lower=min(x),upper = max(x),extendInt = extend)$root)

}
