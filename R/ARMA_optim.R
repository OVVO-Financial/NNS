#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param training.set numeric; \code{NULL} (defualt) Sets the number of variable observations
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param linear.approximation logical; \code{TRUE} (default) Uses the best linear output from \code{NNS.reg} to generate a nonlinear and mixture regression for comparison.  \code{FALSE} is a more exhaustive search over the objective space.
#' @param depth integer; \code{depth = 1} (default) Sets the level from which further combinations are generated containing only members from prior level's best \code{seasonal.factors}.
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
#' ## Nonlinear NNS.ARMA period optimization using 2 yearly lags on AirPassengers monthly data
#' \dontrun{
#' nns.optims <- NNS.ARMA.optim(AirPassengers[1:132], training.set = 120,
#' seasonal.factor = seq(12, 24, 6))
#'}
#'
#' ## Then use optimal parameters in NNS.ARMA to predict 12 periods in-sample
#' \dontrun{
#' NNS.ARMA(AirPassengers, h=12, training.set=132,
#' seasonal.factor = nns.optims$periods, method = nns.optims$method)
#' }
#'
#' @export

NNS.ARMA.optim=function(variable, training.set,
                        seasonal.factor,
                        negative.values = FALSE,
                        obj.fn = expression(sum((predicted - actual)^2)),
                        objective = "min",
                        linear.approximation = TRUE,
                        depth = 1,
                        print.trace = TRUE){

    training.set = floor(training.set)

    if(!is.null(training.set)){
        seasonal.factor = seasonal.factor[seasonal.factor<training.set/3 & seasonal.factor>0]
  } else {
        seasonal.factor = seasonal.factor[seasonal.factor<length(variable)/3 & seasonal.factor>0]
  }


  variable = as.numeric(variable)

  h = length(variable) - training.set
  actual = tail(variable, h)

  nns.estimates = list()
  seasonal.combs = list()

  overall.seasonals = list()
  overall.estimates = list()

  previous.estimates = list()
  previous.seasonals = list()

for(j in c('lin','nonlin','both')){
    current.seasonals = list()
    current.estimate = numeric()
    seasonal.combs = list()

    for(i in 1 : length(seasonal.factor)){
        nns.estimates.indiv = numeric()

            seasonal.combs[[i]] = combn(seasonal.factor, i)


                if(i == 1){
                    if(linear.approximation  && j!="lin"){
                        seasonal.combs[[1]] = matrix(unlist(overall.seasonals[[1]]),ncol=1)
                        current.seasonals[[1]] = unlist(overall.seasonals[[1]])
                    } else {current.seasonals[[i]] = as.numeric(unlist(seasonal.combs[[1]]))}
                } else {
                    if(linear.approximation  && j!="lin"){
                        next
                    } else {
                        current.seasonals[[i]] = as.numeric(unlist(current.seasonals[[i-1]])) }
                }

            if(i > depth){
                seasonal.combs[[i]] = seasonal.combs[[i]][,apply(seasonal.combs[[i]],2, function(z) sum(current.seasonals[[i]]%in%z))==length(current.seasonals[[i]])]
            }



            if(is.null(ncol(seasonal.combs[[i]]))){ break }
            if(dim(seasonal.combs[[i]])[2]==0){ break }

            for(k in 1 : ncol(seasonal.combs[[i]])){
                # Find the min (obj.fn) for a given seasonals sequence
                predicted = NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = seasonal.combs[[i]][ , k], method = j, plot = FALSE, negative.values = negative.values)

                nns.estimates.indiv[k] = eval(obj.fn)
                if(objective=='min'){
                    nns.estimates.indiv[is.na(nns.estimates.indiv)] = Inf
                } else {
                    nns.estimates.indiv[is.na(nns.estimates.indiv)] = -Inf
                }
            } # for k in ncol(seasonal.combs)


            nns.estimates[[i]] = nns.estimates.indiv
            nns.estimates.indiv = numeric()

            if(objective=='min'){
                current.seasonals[[i]] = seasonal.combs[[i]][,which.min(nns.estimates[[i]])]
                current.estimate[i] = min(nns.estimates[[i]])
                if(i > 1 && current.estimate[i] > current.estimate[i-1]){
                  current.seasonals = current.seasonals[-length(current.estimate)]
                  current.estimate = current.estimate[-length(current.estimate)]
                  break
                }
            } else {
                current.seasonals[[i]] = seasonal.combs[[i]][,which.max(nns.estimates[[i]])]
                current.estimate[i] = max(nns.estimates[[i]])
                if(i > 1 && current.estimate[i] <= current.estimate[i-1]){
                  current.seasonals = current.seasonals[-length(current.estimate)]
                  current.estimate = current.estimate[-length(current.estimate)]
                  break
                }
            }


            if(print.trace){
                if(i == 1){
                  print(j)
                  print("COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:")
                }
                print(paste("CURRENT method = ", paste0("'",j,"'"), ", seasonal.factor = ", paste("c(", paste(unlist(current.seasonals[[i]]), collapse = ", ")),")"))
                print(paste0("CURRENT ", j, " OBJECTIVE FUNCTION = ", current.estimate[i]))
            }


### BREAKING PROCEDURE FOR IDENTICAL PERIODS ACROSS METHODS
        if(which(c("lin",'nonlin','both')==j) > 1 ){
            if(sum(as.numeric(unlist(current.seasonals[[i]]))%in%as.numeric(unlist(previous.seasonals[[which(c("lin",'nonlin','both')==j)-1]][i])))==length(as.numeric(unlist(current.seasonals[[i]])))){

              if(objective=='min'){
                if(current.estimate[i] >= previous.estimates[[which(c("lin",'nonlin','both')==j)-1]][i]) break
                } else{
                  if(current.estimate[i] <= previous.estimates[[which(c("lin",'nonlin','both')==j)-1]][i])  break
                }
            }
        }

    } # for i in 1:length(seasonal factor)

    previous.seasonals[[which(c("lin",'nonlin','both')==j)]] = current.seasonals
    previous.estimates[[which(c("lin",'nonlin','both')==j)]] = current.estimate

    overall.seasonals[[which(c("lin",'nonlin','both')==j)]] = current.seasonals[length(current.estimate)]
    overall.estimates[[which(c("lin",'nonlin','both')==j)]] = current.estimate[length(current.estimate)]

    if(print.trace){
        if(i > 1){
            print(paste0("BEST method = ", paste0("'",j,"'"),  ", seasonal.factor = ", paste("c(", paste(unlist(current.seasonals[length(current.estimate)]), collapse = ", ")),")"))
            print(paste0("BEST ", j, " OBJECTIVE FUNCTION = ", current.estimate[length(current.estimate)]))
        }
      else {
            print(paste0("BEST method = ", paste0("'",j,"'"), " PATH MEMBER = ", paste("c(", paste(unlist(current.seasonals), collapse = ", ")),")"))
            print(paste0("BEST ", j, " OBJECTIVE FUNCTION = ", current.estimate[1]))
        }
    }



} # for j in c('lin','nonlin','both')

        if(objective=='min'){
            nns.periods = unlist(overall.seasonals[[which.min(unlist(overall.estimates))]])
            nns.method = c("lin","nonlin","both")[which.min(unlist(overall.estimates))]
            nns.SSE = min(unlist(overall.estimates))
        } else {
            nns.periods = unlist(overall.seasonals[[which.max(unlist(overall.estimates))]])
            nns.method = c("lin","nonlin","both")[which.max(unlist(overall.estimates))]
            nns.SSE = max(unlist(overall.estimates))
        }

return(list(periods = nns.periods,
            obj.fn = nns.SSE,
            method = nns.method))
}
