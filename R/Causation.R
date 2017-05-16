#' NNS Causation
#'
#' Returns the causality from observational data between two variables
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param tau integer; Number of lagged observations to consider.
#' @param plot logical; \code{FALSE} (default) Plots the raw variables, tau normalized, and cross-normalized variables.
#' @return Returns the directional causation (x ---> y) or (y ---> x) and net quantity of association.  For causal matrix, gross quantity of association is returned as (x[column] ---> y[row]).
#' @keywords causation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## x clearly causes y...
#' set.seed(123)
#' x<-rnorm(100); y<-x^2
#' NNS.caus(x,y,1)
#'
#' ## Causal matrix
#' \dontrun{
#' NNS.caus(data.matrix(iris),tau = 0)
#' }
#' @export

NNS.caus <- function(x,y,tau,plot=FALSE){

  if(!missing(y)){
  Causation.x.given.y = Uni.caus(x,y,tau,plot = plot)
  Causation.y.given.x = Uni.caus(y,x,tau,plot = plot)

  if(abs(Causation.x.given.y)<=abs(Causation.y.given.x)){
      return(c(Causation.x.given.y = Causation.x.given.y,
               Causation.y.given.x = Causation.y.given.x,
        "C(y--->x)" =  Causation.y.given.x-Causation.x.given.y))
  } else {
      return(c(Causation.x.given.y = Causation.x.given.y,
                Causation.y.given.x = Causation.y.given.x,
                "C(x--->y)" = Causation.x.given.y-Causation.y.given.x))
    }
  } else {

    causes = list()

    for(i in 1:ncol(x)){
      causes[[i]]=sapply(1:ncol(x), function(b) Uni.caus(x[,i],x[,b],plot = F,tau=tau))
    }

    causes=do.call(cbind,causes)
    diag(causes)=1

    colnames(causes)=colnames(x)
    rownames(causes)=colnames(x)

    return(causes)
}


}
