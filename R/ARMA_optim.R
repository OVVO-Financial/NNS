#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param training.set numeric; \code{NULL} (defualt) Sets the number of variable observations
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor=c(12,24,36))}
#' @param method options: ("comb","seq"); \code{"comb"} Tries all combinations of \code{"seasonal.factor"} provided, while \code{"seq"} (default) tests each by adding element to previous iteration of \code{"seasonal.factor"}.
#'
#' @return Returns a list containing a vector of optimal seasonal periods \code{$period}, the minimum SSE value \code{$SSE}, and the \code{$method} identifying which \link{NNS.ARMA} method was used.
#'
#' @note The number of combinations will grow prohibitively large, they should be kept to a minimum when \code{(method="comb")}.
#'
#' \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
#' @keywords Autoregressive model
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#'
#' ## Nonlinear NNS.ARMA period optimization using 5 yearly lags on AirPassengers monthly data
#' \dontrun{
#' nns.optims <- NNS.ARMA.optim(AirPassengers,training.set=132,seasonal.factor=seq(12,60,12))
#'}
#'
#' ## Then use optimal parameters in NNS.ARMA
#' \dontrun{
#' NNS.ARMA(AirPassengers,training.set=132,seasonal.factor=nns.optims$periods,method=nns.optims$method)
#' }
#'
#' @export

NNS.ARMA.optim=function(variable,training.set,seasonal.factor,method='seq'){

  limit=length(variable)/2
  if(max(seasonal.factor)>=limit){stop(paste0("Please set maximum [seasonal.factor] to less than ",floor(limit)))}

  variable = as.numeric(variable)

  h=length(variable)-training.set
  test.set=tail(variable,h)

  nns.estimates=list()
  seasonal.combs=list()

for(i in 1:length(seasonal.factor)){
  nns.estimates.indiv=numeric()

  # Combinations of seasonal factors
  if(method=='comb'){

    seasonal.combs[[i]]=combn(seasonal.factor,i)

    for(k in 1:ncol(seasonal.combs[[i]])){

      # Find the min SSE for a given seasonals sequence
      nns.estimates.indiv[k]=sum((NNS.ARMA(variable,training.set = training.set,h=h,seasonal.factor =seasonal.combs[[i]][,k],method='lin',plot=FALSE)-test.set)^2)
    }

    nns.estimates[[i]]=nns.estimates.indiv

    nns.estimates.indiv=numeric()

    }
  else{
    nns.estimates.indiv[i]=sum((NNS.ARMA(variable,training.set = training.set,h=h,seasonal.factor =seasonal.factor[1:i],method='lin',plot=FALSE)-test.set)^2)
  }


}
  methods.SSE=numeric()
  if(method=='comb'){
      min.estimate=sapply(nns.estimates,min)
      min.estimate.index=which.min(min.estimate)
      min.estimate.entry=which.min(nns.estimates[[min.estimate.index]])

      nns.periods=seasonal.combs[[min.estimate.index]][,min.estimate.entry]
      benchmark.SSE=min(min.estimate)


      for(j in c("lin","both","nonlin")){
        methods.SSE[j]=sum((NNS.ARMA(variable,training.set = training.set,h=h,seasonal.factor =nns.periods,method=j,plot=FALSE)-test.set)^2)
      }

      nns.SSE=min(methods.SSE)
      nns.method=c("lin","both","nonlin")[which.min(methods.SSE)]


  } else{
      min.estimate=which.min(nns.estimates.indiv)
      nns.periods=seasonal.factor[1:min.estimate]
      benchmark.SSE=nns.estimates.indiv[min.estimate]

      for(j in c("lin","both","nonlin")){
        methods.SSE[j]=sum((NNS.ARMA(variable,training.set = training.set,h=h,seasonal.factor =nns.periods,method=j,plot=FALSE)-test.set)^2)
      }

      nns.SSE=min(methods.SSE)
      nns.method=c("lin","both","nonlin")[which.min(methods.SSE)]

  }

  return(list(periods=nns.periods,SSE=nns.SSE,method=nns.method))
}
