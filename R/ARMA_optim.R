#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param training.set numeric; \code{NULL} (defualt) Sets the number of variable observations
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param method options: ("comb", "seq"); \code{"comb"} (default) Tries all combinations of \code{"seasonal.factor"} provided, while \code{"seq"}  tests each by adding element to previous iteration of \code{"seasonal.factor"}.
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param depth integer; 1 (default) Sets the level from which further combinations are generated containing only members from prior level's best \code{seasonal.factors}.
#' @param print.trace logical; \code{TRUE} (defualt) Prints current iteration information.  Suggested as backup in case of error, best parameters to that point still known and copyable!
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

NNS.ARMA.optim=function(variable, training.set,
                        seasonal.factor,
                        method = "comb",
                        negative.values = FALSE,
                        obj.fn = expression(sum((predicted - actual)^2)),
                        depth = 1,
                        print.trace = TRUE){

  if(!is.null(training.set)){
    seasonal.factor = seasonal.factor[seasonal.factor<training.set/exp(1) & seasonal.factor>0]
  } else {
    seasonal.factor = seasonal.factor[seasonal.factor<length(variable)/exp(1) & seasonal.factor>0]
  }

  variable = as.numeric(variable)

  h = length(variable) - training.set
  actual = tail(variable, h)

  nns.estimates = list()
  seasonal.combs = list()
  best.path = list()
  seasonal.methods = list()
  estimate.methods = list()

for(j in c('lin','nonlin','both')){

  for(i in 1 : length(seasonal.factor)){
    nns.estimates.indiv = numeric()

    # Combinations of seasonal factors
    if(method == 'comb'){

      seasonal.combs[[i]] = combn(seasonal.factor, i)

      if(i > depth){
        best.path[[i]] = unlist(best.path[[i-1]])
        seasonal.combs[[i]] = seasonal.combs[[i]][,apply(seasonal.combs[[i]],2, function(z) sum(best.path[[i]]%in%z))==length(best.path[[i]])]
      }

      if(is.null(ncol(seasonal.combs[[i]]))){ break }

      for(k in 1 : ncol(seasonal.combs[[i]])){

          # Find the min SSE for a given seasonals sequence
          predicted = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = seasonal.combs[[i]][ , k], method = j, plot = FALSE, negative.values = negative.values)

          nns.estimates.indiv[k] = eval(obj.fn)
          nns.estimates.indiv[is.na(nns.estimates.indiv)] = Inf
      } # for k in ncol(seasonal.combs)


      nns.estimates[[i]] = nns.estimates.indiv
      nns.estimates.indiv = numeric()

      best.path[[i]] = seasonal.combs[[i]][,which.min(nns.estimates[[i]])]
      nns.estimates[[i]] = min(nns.estimates[[i]])

      if(print.trace){
        print(paste("CURRENT ", j, " PATH MEMBER = ", paste("c(", paste(unlist(best.path[[i]]), collapse = ", ")),")"))
        print(paste0("CURRENT ", j, " PREDICTION = ", nns.estimates[[i]]))
      }


      if(i > 1 && (nns.estimates[[i]] > nns.estimates[[i-1]])) break

    } else {
        predicted = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = seasonal.factor[1 : i], method = 'lin', plot = FALSE, negative.values = negative.values)

        nns.estimates.indiv = eval(obj.fn)
        nns.estimates.indiv[is.na(nns.estimates.indiv)] = Inf
    }


  } # for i in 1:length(seasonal factor)


  # Need 3 estimates with their seasonal.combs...
  seasonal.methods[[which(c("lin",'nonlin','both')==j)]] = best.path[[i-1]]
  estimate.methods[[which(c("lin",'nonlin','both')==j)]] = nns.estimates[[i-1]]

  if(print.trace){
    print(paste0("BEST ", j, " PATH MEMBER = ", paste("c(", paste(unlist(best.path[[i-1]]), collapse = ", ")),")"))
    print(paste0("BEST ", j, " PREDICTION = ", nns.estimates[[i-1]]))
  }
} # for j in c('lin','nonlin','both')



  if(method == 'comb'){
    nns.periods = seasonal.methods[[which.min(unlist(estimate.methods))]]
    nns.method = c("lin","nonlin","both")[which.min(unlist(estimate.methods))]
    nns.SSE = min(unlist(estimate.methods))
  } else {
    methods.SSE = numeric()
    min.estimate = which.min(nns.estimates.indiv)
    nns.periods = seasonal.factor[1 : min.estimate]
    benchmark.SSE = nns.estimates.indiv[min.estimate]

    for(j in c("lin", "both", "nonlin")){predictions = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = j, plot = FALSE, negative.values = negative.values)

      methods.SSE[j] = eval(obj.fn)
      methods.SSE[is.na(methods.SSE)] = Inf

    }

    nns.SSE = min(methods.SSE)
    nns.method = c("lin", "both", "nonlin")[which.min(methods.SSE)]

  }

  return(list(periods = nns.periods,
              obj.fn = nns.SSE,
              method = nns.method))
}
