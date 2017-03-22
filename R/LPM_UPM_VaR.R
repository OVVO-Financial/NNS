#' LPM VaR
#'
#' Generates a VaR based on the Lower Partial Moment ratio
#' @param percentile numeric [0,1]; The percentile for VaR.
#' @param degree integer; \code{(degree=0)} for discrete distributions, \code{(degree=1)} for continuous distributions.
#' @param x a numeric vector.
#' @keywords VaR
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' LPM.VaR(0.95,0,x)
#' @export

LPM.VaR <- function(percentile,degree,x){

  f<- function(tgt) (sum((tgt - (x[x <= tgt]))^degree)/length(x))/((sum((tgt - (x[x <= tgt]))^degree)/length(x))+ sum(((x[x > tgt]) - tgt)^degree)/length(x))- (1-percentile)

  return(uniroot(f,lower=min(x),upper = max(x))$root)

}

#' UPM VaR
#'
#' Generates an upside VaR based on the Upper Partial Moment ratio
#' @param percentile numeric [0,1]; The percentile for VaR.
#' @param degree integer; \code{(degree=0)} for discrete distributions, \code{(degree=1)} for continuous distributions.
#' @param x a numeric vector.
#' @keywords VaR
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' UPM.VaR(0.95,0,x)
#' @export

UPM.VaR <- function(percentile,degree,x){


  f<- function(tgt) (sum(((x[x > tgt]) - tgt)^degree)/length(x))/((sum((tgt - (x[x <= tgt]))^degree)/length(x))+ sum(((x[x > tgt]) - tgt)^degree)/length(x))- (1-percentile)

  return(uniroot(f,lower=min(x),upper = max(x))$root)

}
