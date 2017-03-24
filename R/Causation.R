#' NNS Causation
#'
#' Returns the causality from observational data between two variables
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param tau integer; Number of lagged observations to consider.
#' @param plot logical; \code{FALSE} (default) Plots the raw variables, tau normalized, and cross-normalized variables.
#' @return Returns the directional causation and quantity of association.
#' @keywords causation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## x clearly causes y...
#' set.seed(123)
#' x<-rnorm(100); y<-x^2
#' NNS.caus(x,y,1)
#' @export

NNS.caus <- function(x,y,tau,plot=FALSE){
  Causation.x.given.y = Uni.caus(x,y,tau,plot = plot)
  Causation.y.given.x = Uni.caus(y,x,tau,plot = plot)

  if(abs(Causation.x.given.y)<abs(Causation.y.given.x)){
      return(c(Causation.x.given.y = Causation.x.given.y,
               Causation.y.given.x = Causation.y.given.x,
        "C(y--->x)" =  Causation.y.given.x-Causation.x.given.y))
  } else {
      return(c(Causation.x.given.y = Causation.x.given.y,
                Causation.y.given.x = Causation.y.given.x,
                "C(x--->y)" = Causation.x.given.y-Causation.y.given.x))
    }
}
