#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param training.set integer; Sets the number of variable observations as the training set.  See \code{Note} below for recommended uses.
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param linear.approximation logical; \code{TRUE} (default) Uses the best linear output from \code{NNS.reg} to generate a nonlinear and mixture regression for comparison.  \code{FALSE} is a more exhaustive search over the objective space.
#' @param lin.only logical; \code{FALSE} (default) Restricts the optimization to linear methods only.
#' @param print.trace logical; \code{TRUE} (defualt) Prints current iteration information.  Suggested as backup in case of error, best parameters to that point still known and copyable!
#' @param ncores integer; value specifying the number of cores to be used in the parallelized procedure. If NULL (default), the number of cores to be used is equal to half the number of cores of the machine.
#'
#' @return Returns a list containing:
#' \itemize{
#' \item{\code{$period}} a vector of optimal seasonal periods
#' \item{\code{$weights}} the optimal weights of each seasonal period between an equal weight or NULL weighting
#' \item{\code{$obj.fn}} the objective function value
#' \item{\code{$method}} the method identifying which \link{NNS.ARMA} method was used.
#' \item{\code{$bias.shift}} a numerical result of the overall bias of the optimum objective function result.  To be added to the final result when using the \link{NNS.ARMA} with the derived parameters.
#'}
#' @note
#' \itemize{
#' \item{} Typically, \code{(training.set = length(variable) - 2 * length(forecast horizon))} is used for optimization.  Smaller samples would use \code{(training.set = length(variable) - length(forecast horizon))} in order to preserve information.
#'
#' \item{} The number of combinations will grow prohibitively large, they should be kept as small as possible.  \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
#' \item{} If variable cannot logically assume negative values, then the \code{$bias.shift} must be limited to 0 via a \code{pmax(0,...)} call.
#'}
#'
#'
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#'
#' ## Nonlinear NNS.ARMA period optimization using 2 yearly lags on AirPassengers monthly data
#' \dontrun{
#' nns.optims <- NNS.ARMA.optim(AirPassengers[1:132], training.set = 120,
#' seasonal.factor = seq(12, 24, 6))
#'
#' ## Then use optimal parameters in NNS.ARMA to predict 12 periods in-sample.
#' ## Note the {$bias.shift} usage in the {NNS.ARMA} function:
#' nns.estimates <- NNS.ARMA(AirPassengers, h = 12, training.set = 132,
#' seasonal.factor = nns.optims$periods, method = nns.optims$method) + nns.optims$bias.shift
#'
#' ## If variable cannot logically assume negative values
#' nns.estimates <- pmax(0, nns.estimates)
#' }
#'
#' @export

NNS.ARMA.optim <- function(variable, training.set,
                        seasonal.factor,
                        negative.values = FALSE,
                        obj.fn = expression(sum((predicted - actual)^2)),
                        objective = "min",
                        linear.approximation = TRUE,
                        lin.only = FALSE,
                        print.trace = TRUE,
                        ncores = NULL){

  if(any(class(variable)=="tbl")) variable <- as.vector(unlist(variable))

  if(length(training.set)==length(variable)){ stop("Please provide a 'training.set' value (integer) less than the length of the variable.")}

  if(is.null(obj.fn)){ stop("Please provide an objective function")}
  objective <- tolower(objective)

  if (is.null(ncores)) {
      num_cores <- as.integer(parallel::detectCores()) - 1
  } else {
      num_cores <- ncores
  }

  if(is.null(training.set)){stop("Please use the length of the variable less the forecast period as the training.set value.")}

  training.set <- floor(training.set)

  variable <- as.numeric(variable)

  h <- length(variable) - training.set
  actual <- tail(variable, h)

  if(is.null(training.set)) l <- length(variable) else l <- training.set

  denominator <- min(5, max(2, ifelse((l/100)%%1 < .5, floor(l/100), ceiling(l/100))))

  seasonal.factor <- seasonal.factor[seasonal.factor <= (l/denominator)]
  seasonal.factor <- unique(seasonal.factor)

  if(length(seasonal.factor)==0) seasonal.factor <- 1

  oldw <- getOption("warn")
  options(warn = -1)

  nns.estimates <- list()
  seasonal.combs <- list()

  overall.seasonals <- list()
  overall.estimates <- list()

  previous.estimates <- list()
  previous.seasonals <- list()

  if(lin.only) methods <- "lin" else methods <- c('lin','nonlin','both')

  for(j in methods){
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

          predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = unlist(overall.seasonals[[1]]), method = j, plot = FALSE, negative.values = negative.values, ncores = 1)

          nns.estimates.indiv <- eval(obj.fn)

      } else {

          if(num_cores>1){
            cl <- parallel::makeCluster(num_cores)
            doParallel::registerDoParallel(cl)
          }

          nns.estimates.indiv <- foreach(k = 1 : ncol(seasonal.combs[[i]]),.packages = c("NNS", "data.table"))%dopar%{
          actual <- tail(variable, h)

          predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor =  seasonal.combs[[i]][ , k], method = j, plot = FALSE, negative.values = negative.values, ncores = 1)

          nns.estimates.indiv <- eval(obj.fn)

        }

        if(num_cores>1){
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
    overall.estimates[[which(c("lin",'nonlin','both')==j)]] <- current.estimate[length(current.estimate)]


    if(print.trace){
        if(i > 1){
            print(paste0("BEST method = ", paste0("'",j,"'"),  ", seasonal.factor = ", paste("c(", paste(unlist(current.seasonals[length(current.estimate)]), collapse = ", "))," )"))
            print(paste0("BEST ", j, " OBJECTIVE FUNCTION = ", current.estimate[length(current.estimate)]))
        } else {
            print(paste0("BEST method = ", paste0("'",j,"'"), " PATH MEMBER = ", paste("c(", paste(unlist(current.seasonals), collapse = ", "))," )"))
            print(paste0("BEST ", j, " OBJECTIVE FUNCTION = ", current.estimate[1]))
        }
    }

    gc()

    if(lin.only) predicted <- current.estimate
    if(j!='lin' && lin.only){ break }

  } # for j in c('lin','nonlin','both')



  if(objective=='min'){
      nns.periods <- unlist(overall.seasonals[[which.min(unlist(overall.estimates))]])
      nns.method <- c("lin","nonlin","both")[which.min(unlist(overall.estimates))]
      nns.SSE <- min(unlist(overall.estimates))

    if(length(nns.periods)>1){
        predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, ncores = 1, weights = rep((1/length(nns.periods)),length(nns.periods)))

        weight.SSE <- eval(obj.fn)

        if(weight.SSE<nns.SSE){
            nns.weights <- rep((1/length(nns.periods)),length(nns.periods))

            bias <- mode(predicted - actual)
            predicted <- predicted+bias
            bias.SSE <- eval(obj.fn)

            if(bias.SSE>weight.SSE){bias <- 0}

        } else {
            nns.weights <- NULL

            bias <- mode(predicted - actual)
            predicted <- predicted+bias
            bias.SSE <- eval(obj.fn)

            if(bias.SSE>nns.SSE){bias <- 0}
        }
    } else {
        nns.weights <- NULL

        bias <- mode(predicted - actual)
        predicted <- predicted+bias
        bias.SSE <- eval(obj.fn)

        if(bias.SSE>nns.SSE){bias <- 0}
    }

  } else {
      nns.periods <- unlist(overall.seasonals[[which.max(unlist(overall.estimates))]])
      nns.method <- c("lin","nonlin","both")[which.max(unlist(overall.estimates))]
      nns.SSE <- max(unlist(overall.estimates))

      if(length(nns.periods)>1){
          predicted <- NNS.ARMA(variable, training.set = training.set, h = h, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, ncores = 1, weights = rep((1/length(nns.periods)),length(nns.periods)))

          weight.SSE <- eval(obj.fn)

          if(weight.SSE>nns.SSE){
              nns.weights <- rep((1/length(nns.periods)),length(nns.periods))

              bias <- mode(predicted - actual)
              predicted <- predicted+bias
              bias.SSE <- eval(obj.fn)

              if(bias.SSE<weight.SSE){bias <- 0}

          } else {
              nns.weights <- NULL

              bias <- mode(predicted - actual)
              predicted <- predicted+bias
              bias.SSE <- eval(obj.fn)
              if(objective=="min"){
                if(bias.SSE<=nns.SSE){bias <- 0}
              } else {
                if(bias.SSE>=nns.SSE){bias <- 0}
              }
          }
      } else {
          nns.weights <- NULL

          bias <- mode(predicted - actual)
          predicted <- predicted+bias
          bias.SSE <- eval(obj.fn)
          if(objective=="min"){
              if(bias.SSE<=nns.SSE){bias <- 0}
          } else {
              if(bias.SSE>=nns.SSE){bias <- 0}
          }


      }

  }

  options(warn = oldw)

  return(list(periods = nns.periods,
              weights = nns.weights,
              obj.fn = nns.SSE,
              method = nns.method,
              bias.shift = bias))
}
