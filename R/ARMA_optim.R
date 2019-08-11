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
#' @param print.trace logical; \code{TRUE} (defualt) Prints current iteration information.  Suggested as backup in case of error, best parameters to that point still known and copyable!
#' @param ncores integer; value specifying the number of cores to be used in the parallelized procedure. If NULL (default), the number of cores to be used is equal to half the number of cores of the machine.
#' @param subcores integer; value specifying the number of cores to be used in the parallelized procedure in the subroutine \link{NNS.ARMA}.  If NULL (default), the number of cores to be used is equal to half the number of cores of the machine - 1.
#'
#' @return Returns a list containing:
#' \itemize{
#' \item{\code{$period}} a vector of optimal seasonal periods
#' \item{\code{$weights}} the optimal weights of each seasonal period between an equal weight or NULL weighting
#' \item{\code{$obj.fn}} the minimum objective function value
#' \item{\code{$method}} the method identifying which \link{NNS.ARMA} method was used.
#' \item{\code{$ensemble}} a logical indicator representing whether the ensemble method using a naive \link{NNS.ARMA} with \code{best.periods = length(seasonal.factor)} generated a better \code{obj.fn} result.
#'}
#' @note The number of combinations will grow prohibitively large, they should be kept as small as possible.
#'
#' \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#'
#' ## Nonlinear NNS.ARMA period optimization using 2 yearly lags on AirPassengers monthly data
#' \dontrun{
#' nns.optims <- NNS.ARMA.optim(AirPassengers[1:132], training.set = 120,
#' seasonal.factor = seq(12, 24, 6))
#'
#' ## Then use optimal parameters in NNS.ARMA to predict 12 periods in-sample
#' NNS.ARMA(AirPassengers, h = 12, training.set = 132,
#' seasonal.factor = nns.optims$periods, method = nns.optims$method)
#'
#' ## If {$ensemble == TRUE} then take the mean of both forecasts, using {best.periods = length(seasonal.factor)}
#' in the {NNS.ARMA} function.
#'
#' ## Using optimized paramaters from {NNS.ARMA.optim}
#' nns.estimates.1 <- NNS.ARMA(AirPassengers, h = 12, training.set = 132,
#' seasonal.factor = nns.optims$periods, method = nns.optims$method)
#'
#' ## Using all {seasonal.factor} considered in {NNS.ARMA.optim}
#' nns.estimates.2 <- NNS.ARMA(AirPassengers, h = 12, training.set = 132,
#' seasonal.factor = FALSE, best.periods = length(seq(12,24,6)), method = nns.optims$method)
#'
#' ## Combined estimates
#' nns.combined.estimate <- rowMeans(cbind(nns.estimates.1, nns.estimates.2))
#' }
#'
#' @export

NNS.ARMA.optim <- function(variable, training.set,
                        seasonal.factor,
                        negative.values = FALSE,
                        obj.fn = expression(sum((predicted - actual)^2)),
                        objective = "min",
                        linear.approximation = TRUE,
                        print.trace = TRUE,
                        ncores = NULL, subcores = NULL){



  if (is.null(ncores)) {
      cores <- detectCores()
      num_cores <- as.integer(cores / 2)
  } else {
      cores <- detectCores()
      num_cores <- ncores
  }

  if (is.null(subcores)) {
    subcores <- as.integer(cores / 2) - 1
  }

  if((num_cores+subcores)>cores){ stop(paste0("Please ensure total number of cores [ncores + subcores] is less than ", cores))}

  if(is.null(training.set)){training.set <- 0}
  training.set <- floor(training.set)

  variable <- as.numeric(variable)

  h <- length(variable) - training.set
  actual <- tail(variable, h)

  if(is.null(training.set)){
      l <- length(variable)
  } else {
      l <- training.set
  }

  denominator <- min(5,max(2, as.integer(l/100)))

  seasonal.factor <- seasonal.factor[seasonal.factor<=(l/denominator)]

  if(length(seasonal.factor)==0){stop(paste0('Please ensure "seasonal.factor" contains elements less than ', l/denominator, ", otherwise use cross-validation of seasonal factors as demonstrated in the vignette >>> Getting Started with NNS: Forecasting"))}

  nns.estimates <- list()
  seasonal.combs <- list()

  overall.seasonals <- list()
  overall.estimates <- list()

  previous.estimates <- list()
  previous.seasonals <- list()

  for(j in c('lin','nonlin','both')){
      current.seasonals <- list()
      current.estimate <- numeric()
      seasonal.combs <- list()

      for(i in 1 : length(seasonal.factor)){
          nns.estimates.indiv <- list()

          if(i == 1){
              seasonal.combs[[i]] <- t(seasonal.factor)
          } else {
              remaining.index <- !(seasonal.factor%in%current.seasonals[[i-1]])
              if(sum(remaining.index)==0){ break }
              seasonal.combs[[i]] <- rbind(replicate(length(seasonal.factor[remaining.index]), current.seasonals[[i-1]]), as.integer(seasonal.factor[remaining.index]))
          }

      if(i == 1){
          if(linear.approximation  && j!="lin"){
              seasonal.combs[[1]] <- matrix(unlist(overall.seasonals[[1]]),ncol=1)
              current.seasonals[[1]] <- unlist(overall.seasonals[[1]])
          } else {
              current.seasonals[[i]] <- as.integer(unlist(seasonal.combs[[1]]))
          }
      } else {
          if(linear.approximation  && j!="lin"){
              next
          } else {
              current.seasonals[[i]] <- as.integer(unlist(current.seasonals[[i-1]]))
          }
      }



      if(is.null(ncol(seasonal.combs[[i]]))){ break }
      if(dim(seasonal.combs[[i]])[2]==0){ break }



      if(j!="lin" && linear.approximation){

          # Find the min (obj.fn) for a given seasonals sequence
          actual <- tail(variable, h)

          predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = unlist(overall.seasonals[[1]]), method = j, plot = FALSE, negative.values = negative.values, ncores = subcores)

          nns.estimates.indiv <- eval(obj.fn)

      } else {

          if(num_cores>1){
              cl <- makeCluster(num_cores)
              registerDoParallel(cl)
          } else { cl <- NULL }

          nns.estimates.indiv <- foreach(k = 1 : ncol(seasonal.combs[[i]]),.packages = "NNS")%dopar%{
          actual <- tail(variable, h)

          predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor =  seasonal.combs[[i]][ , k], method = j, plot = FALSE, negative.values = negative.values, ncores = subcores)

          eval(obj.fn)

        }

        if(!is.null(cl)){
            stopCluster(cl)
            registerDoSEQ()
        }

      }

      gc()
      nns.estimates.indiv <- unlist(nns.estimates.indiv)

      if(objective=='min'){
          nns.estimates.indiv[is.na(nns.estimates.indiv)] <- Inf
      } else {
          nns.estimates.indiv[is.na(nns.estimates.indiv)] <- -Inf
      }

      nns.estimates[[i]] <- nns.estimates.indiv
      nns.estimates.indiv <- numeric()

      if(objective=='min'){
          current.seasonals[[i]] <- seasonal.combs[[i]][,which.min(nns.estimates[[i]])]
          current.estimate[i] <- min(nns.estimates[[i]])
        if(i > 1 && current.estimate[i] > current.estimate[i-1]){
          current.seasonals <- current.seasonals[-length(current.estimate)]
          current.estimate <- current.estimate[-length(current.estimate)]
          break
        }
      } else {
          current.seasonals[[i]] <- seasonal.combs[[i]][,which.max(nns.estimates[[i]])]
          current.estimate[i] <- max(nns.estimates[[i]])
        if(i > 1 && current.estimate[i] < current.estimate[i-1]){
          current.seasonals <- current.seasonals[-length(current.estimate)]
          current.estimate <- current.estimate[-length(current.estimate)]
          break
        }
      }


      if(print.trace){
        if(i == 1){
          print(paste0("CURRNET METHOD: ",j))
          print("COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:")
        }
        print(paste("NNS.ARMA(... method = ", paste0("'",j,"'"), ", seasonal.factor = ", paste("c(", paste(unlist(current.seasonals[[i]]), collapse = ", ")),") ...)"))
        print(paste0("CURRENT ", j, " OBJECTIVE FUNCTION = ", current.estimate[i]))
      }


      ### BREAKING PROCEDURE FOR IDENTICAL PERIODS ACROSS METHODS
      if(which(c("lin",'nonlin','both')==j) > 1 ){
          if(sum(as.numeric(unlist(current.seasonals[[i]]))%in%as.numeric(unlist(previous.seasonals[[which(c("lin",'nonlin','both')==j)-1]][i])))==length(as.numeric(unlist(current.seasonals[[i]])))){

              if(objective=='min'){
                  if(current.estimate[i] >= previous.estimates[[which(c("lin",'nonlin','both')==j)-1]][i]) break
              } else {
                  if(current.estimate[i] <= previous.estimates[[which(c("lin",'nonlin','both')==j)-1]][i]) break
              }
          }
      }

      if(j!='lin' && linear.approximation){ break }

    } # for i in 1:length(seasonal factor)

    previous.seasonals[[which(c("lin",'nonlin','both')==j)]] <- current.seasonals
    previous.estimates[[which(c("lin",'nonlin','both')==j)]] <- current.estimate

    overall.seasonals[[which(c("lin",'nonlin','both')==j)]] <- current.seasonals[length(current.estimate)]
    overall.estimates[[which(c("lin",'nonlin','both')==j)]] = current.estimate[length(current.estimate)]

    if(print.trace){
        if(i > 1){
            print(paste0("BEST method = ", paste0("'",j,"'"),  ", seasonal.factor = ", paste("c(", paste(unlist(current.seasonals[length(current.estimate)]), collapse = ", ")),")"))
            print(paste0("BEST ", j, " OBJECTIVE FUNCTION = ", current.estimate[length(current.estimate)]))
        } else {
            print(paste0("BEST method = ", paste0("'",j,"'"), " PATH MEMBER = ", paste("c(", paste(unlist(current.seasonals), collapse = ", ")),")"))
            print(paste0("BEST ", j, " OBJECTIVE FUNCTION = ", current.estimate[1]))
        }
    }

    gc()



  } # for j in c('lin','nonlin','both')


  if(objective=='min'){
      nns.periods <- unlist(overall.seasonals[[which.min(unlist(overall.estimates))]])
      nns.method <- c("lin","nonlin","both")[which.min(unlist(overall.estimates))]
      nns.SSE <- min(unlist(overall.estimates))

    if(length(nns.periods)>1){
        predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, ncores = subcores, weights = rep((1/length(nns.periods)),length(nns.periods)))

        weight.SSE <- eval(obj.fn)

        if(weight.SSE<nns.SSE){
            nns.weights <- rep((1/length(nns.periods)),length(nns.periods))

            predicted <- rowMeans(predicted, NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = FALSE, best.periods = length(seasonal.factor), method = nns.method, weights = nns.weights, plot = FALSE, negative.values = negative.values, ncores = subcores))

            ensemble.SSE <- eval(obj.fn)<nns.SSE
        } else {
            nns.weights <- NULL

            predicted <- rowMeans(predicted, NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = FALSE, best.periods = length(seasonal.factor), method = nns.method, weights = nns.weights, plot = FALSE, negative.values = negative.values, ncores = subcores))

            ensemble.SSE <- eval(obj.fn)<nns.SSE
        }
    } else {
        nns.weights <- NULL

        predicted <- rowMeans(predicted, NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = FALSE, best.periods = length(seasonal.factor), method = nns.method, weights = nns.weights, plot = FALSE, negative.values = negative.values, ncores = subcores))

        ensemble.SSE <- eval(obj.fn)<nns.SSE
    }

  } else {
      nns.periods <- unlist(overall.seasonals[[which.max(unlist(overall.estimates))]])
      nns.method <- c("lin","nonlin","both")[which.max(unlist(overall.estimates))]
      nns.SSE <- max(unlist(overall.estimates))

      if(length(nns.periods)>1){
          predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, ncores = subcores, weights = rep((1/length(nns.periods)),length(nns.periods)))

          weight.SSE <- eval(obj.fn)

          if(weight.SSE>nns.SSE){
              nns.weights <- rep((1/length(nns.periods)),length(nns.periods))

              predicted <- rowMeans(predicted, NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = FALSE, best.periods = length(seasonal.factor), method = nns.method, plot = FALSE, negative.values = negative.values, ncores = subcores))

              ensemble.SSE <- eval(obj.fn)>weight.SSE

          } else {
              nns.weights <- NULL

              predicted <- rowMeans(predicted, NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = FALSE, best.periods = length(seasonal.factor), method = nns.method, plot = FALSE, negative.values = negative.values, ncores = subcores))

              ensemble.SSE <- eval(obj.fn)>nns.SSE
          }
      } else {
          nns.weights <- NULL

          predicted <- rowMeans(predicted, NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = FALSE, best.periods = length(seasonal.factor), method = nns.method, plot = FALSE, negative.values = negative.values, ncores = subcores))

          ensemble.SSE <- eval(obj.fn)>nns.SSE
      }

  }



  return(list(periods = nns.periods,
              weights = nns.weights,
              obj.fn = nns.SSE,
              method = nns.method,
              ensemble = ensemble.SSE))
}
