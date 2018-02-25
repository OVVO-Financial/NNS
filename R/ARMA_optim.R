#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of given `seasonal.periods` in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum.
#'
#' @param variable a numeric vector.
#' @param training.set numeric; \code{NULL} (defualt) Sets the number of variable observations
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor=c(12,24,36))}
#'
#' @return Returns a list containing a vector of optimal seasonal periods \code{$period} and the minimum SSE value \code{$SSE}.
#'
#' @note The number of combinations will grow prohibitively large, they should be kept to a minimum.
#'
#' \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
#' @keywords Autoregressive model
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#'
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data
#' \dontrun{
#' NNS.ARMA.optim(AirPassengers,training.set=132,seasonal.periods=seq(12,60,12))
#' }
#'
#' @export

NNS.ARMA.optim=function(variable,training.set,seasonal.factor){

  limit=length(variable)/2
  if(max(seasonal.factor)>=limit){stop(paste0("Please set maximum [seasonal.factor] to less than ",floor(limit)))}

  variable = as.numeric(variable)

  h=length(variable)-training.set
  test.set=tail(variable,h)

  nns.estimates=list()
  seasonal.combs=list()

  seasonals=seasonal.factor

for(i in 1:length(seasonals)){
  # Combinations of seasonal factors
  seasonal.combs[[i]]=combn(seasonals,i)

  nns.estimates.indiv=numeric()

  for(k in 1:ncol(seasonal.combs[[i]])){

    # Find the min SSE for a given seasonals sequence
    nns.estimates.indiv[[k]]=sum((NNS.ARMA(variable,training.set = training.set,h=h,seasonal.factor =seasonal.combs[[i]][,k],method='lin',plot=FALSE)-test.set)^2)
  }

  nns.estimates[[i]]=nns.estimates.indiv

}

min.estimate=sapply(nns.estimates,min)
min.estimate.index=which.min(min.estimate)
min.estimate.entry=which.min(nns.estimates[[min.estimate.index]])

return(list(periods=seasonal.combs[[min.estimate.index]][,min.estimate.entry],SSE=min(min.estimate)))

}
