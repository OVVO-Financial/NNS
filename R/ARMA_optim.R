#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param training.set numeric; \code{NULL} (defualt) Sets the number of variable observations
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param method options: ("comb", "seq"); \code{"comb"} Tries all combinations of \code{"seasonal.factor"} provided, while \code{"seq"} (default) tests each by adding element to previous iteration of \code{"seasonal.factor"}.
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#'
#' @return Returns a list containing a vector of optimal seasonal periods \code{$period}, the minimum objective function value \code{$obj.fn}, and the \code{$method} identifying which \link{NNS.ARMA} method was used.
#'
#' @note The number of combinations will grow prohibitively large, they should be kept to a minimum when \code{(method = "comb")}.
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
#' nns.optims <- NNS.ARMA.optim(AirPassengers, training.set = 132, seasonal.factor = seq(12, 36, 12))
#'}
#'
#' ## Then use optimal parameters in NNS.ARMA
#' \dontrun{
#' NNS.ARMA(AirPassengers, training.set=132, seasonal.factor = nns.optims$periods,
#' method = nns.optims$method)
#' }
#'
#' @export

NNS.ARMA.optim=function(variable, training.set, seasonal.factor, method = "seq", negative.values = FALSE, obj.fn = expression(sum((predicted - actual)^2))){

  if(!is.null(training.set)){
    seasonal.factor = seasonal.factor[seasonal.factor<training.set/exp(1)]
  } else {
    seasonal.factor = seasonal.factor[seasonal.factor<length(variable)/exp(1)]
  }

  variable = as.numeric(variable)

  h = length(variable) - training.set
  actual = tail(variable, h)

  nns.estimates = list()
  seasonal.combs = list()
  best.path = list()

  for(i in 1 : length(seasonal.factor)){
    nns.estimates.indiv = numeric()

    # Combinations of seasonal factors
    if(method == 'comb'){

      seasonal.combs[[i]] = combn(seasonal.factor, i)

      if(i > 2){
        best.path = unlist(best.path[[i-1]])
        seasonal.combs[[i]] = seasonal.combs[[i]][,apply(seasonal.combs[[i]],2, function(z) sum(best.path%in%z))==length(best.path)]
      }

      for(k in 1 : ncol(seasonal.combs[[i]])){

        # Find the min SSE for a given seasonals sequence
        predicted = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = seasonal.combs[[i]][ , k], method = 'lin', plot = FALSE, negative.values = negative.values)

        nns.estimates.indiv[k] = eval(obj.fn)
        nns.estimates.indiv[is.na(nns.estimates.indiv)] = Inf
      }

      nns.estimates[[i]] = nns.estimates.indiv
      nns.estimates.indiv = numeric()
      best.path[[i]] = seasonal.combs[[i]][which.min(nns.estimates[[i]])]


      if(i > 1 && (min(nns.estimates[[i]]) > min(nns.estimates[[i-1]]))) break

    } else {
      predicted = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = seasonal.factor[1 : i], method = 'lin', plot = FALSE, negative.values = negative.values)

      nns.estimates.indiv = eval(obj.fn)
      nns.estimates.indiv[is.na(nns.estimates.indiv)] = Inf
    }


  }
  methods.SSE = numeric()
  if(method == 'comb'){
    min.estimate = sapply(nns.estimates,min)
    min.estimate.index = which.min(min.estimate)
    min.estimate.entry = which.min(nns.estimates[[min.estimate.index]])

    nns.periods = seasonal.combs[[min.estimate.index]][ , min.estimate.entry]
    benchmark.SSE = min(min.estimate)



    for(j in c("lin", "both", "nonlin")){
      predictions = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = j, plot = FALSE, negative.values = negative.values)

      methods.SSE[j] = eval(obj.fn)
      methods.SSE[is.na(methods.SSE)] = Inf
  #    if(j==2 && (methods.SSE[2] > methods.SSE[1])) break
    }

    nns.SSE = min(methods.SSE)
    nns.method = c("lin", "both", "nonlin")[which.min(methods.SSE)]

  } else{
    min.estimate = which.min(nns.estimates.indiv)
    nns.periods = seasonal.factor[1 : min.estimate]
    benchmark.SSE = nns.estimates.indiv[min.estimate]

    for(j in c("lin", "both", "nonlin")){predictions = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = j, plot = FALSE, negative.values = negative.values)

      methods.SSE[j] = eval(obj.fn)
      methods.SSE[is.na(methods.SSE)] = Inf
  #    if(j==2 && (methods.SSE[2] > methods.SSE[1])) break
    }

    nns.SSE = min(methods.SSE)
    nns.method = c("lin", "both", "nonlin")[which.min(methods.SSE)]

  }

  return(list(periods = nns.periods,
              obj.fn = nns.SSE,
              method = nns.method))
}
