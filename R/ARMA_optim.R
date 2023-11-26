#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  Minimum sum of squared errors (forecast-actual) is used to determine optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param h integer; \code{NULL} (default) Number of periods to forecast out of sample.  If \code{NULL}, \code{h = length(variable) - training.set}.
#' @param training.set integer; \code{NULL} (default) Sets the number of variable observations as the training set.  See \code{Note} below for recommended uses.
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.  It will automatically select \code{(negative.values = TRUE)} if the minimum value of the \code{variable} is negative.
#' @param obj.fn expression;
#' \code{expression(cor(predicted, actual, method = "spearman") / sum((predicted - actual)^2))} (default) Rank correlation / sum of squared errors is the default objective function.  Any \code{expression(...)} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"max"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param linear.approximation logical; \code{TRUE} (default) Uses the best linear output from \code{NNS.reg} to generate a nonlinear and mixture regression for comparison.  \code{FALSE} is a more exhaustive search over the objective space.
#' @param lin.only logical; \code{FALSE} (default) Restricts the optimization to linear methods only.
#' @param pred.int numeric [0, 1]; \code{NULL} Returns the associated prediction intervals for the final estimate.  Constructed using the maximum entropy bootstrap \link{meboot} on the final estimates.
#' @param print.trace logical; \code{TRUE} (defualt) Prints current iteration information.  Suggested as backup in case of error, best parameters to that point still known and copyable!
#'
#' @return Returns a list containing:
#' \itemize{
#' \item{\code{$period}} a vector of optimal seasonal periods
#' \item{\code{$weights}} the optimal weights of each seasonal period between an equal weight or NULL weighting
#' \item{\code{$obj.fn}} the objective function value
#' \item{\code{$method}} the method identifying which \link{NNS.ARMA} method was used.
#' \item{\code{$shrink}} whether to use the \code{shrink} parameter in \link{NNS.ARMA}.
#' \item{\code{$nns.regress}} whether to smooth the variable via \link{NNS.reg} before forecasting.
#' \item{\code{$bias.shift}} a numerical result of the overall bias of the optimum objective function result.  To be added to the final result when using the \link{NNS.ARMA} with the derived parameters.
#' \item{\code{$errors}} a vector of model errors from internal calibration.
#' \item{\code{$results}} a vector of length \code{h}.
#' \item{\code{$lower.pred.int}} a vector of lower prediction intervals per forecast point.
#' \item{\code{$upper.pred.int}} a vector of upper prediction intervals per forecast point.
#'}
#' @note
#' \itemize{
#' \item{} Typically, \code{(training.set = length(variable) - 2 * length(forecast horizon))} is used for optimization.  Smaller samples would use \code{(training.set = length(variable) - length(forecast horizon))} in order to preserve information.
#'
#' \item{} The number of combinations will grow prohibitively large, they should be kept as small as possible.  \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
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
#' ## To use optimal parameters in NNS.ARMA to predict 12 periods in-sample.
#' ## Note the {$bias.shift} usage in the {NNS.ARMA} function:
#' nns.estimates <- NNS.ARMA(AirPassengers, h = 12, training.set = 132,
#' seasonal.factor = nns.optims$periods, method = nns.optims$method) + nns.optims$bias.shift
#'
#' ## If variable cannot logically assume negative values
#' nns.estimates <- pmax(0, nns.estimates)
#'
#'
#' ## To predict out of sample using best parameters:
#' NNS.ARMA.optim(AirPassengers[1:132], h = 12, seasonal.factor = seq(12, 24, 6))
#' 
#' ## Incorporate any objective function from external packages (such as \code{Metrics::mape})
#' NNS.ARMA.optim(AirPassengers[1:132], h = 12, seasonal.factor = seq(12, 24, 6),
#' obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")
#' }
#'
#' @export

NNS.ARMA.optim <- function(variable,
                           h = NULL,
                           training.set = NULL,
                           seasonal.factor,
                           negative.values = FALSE,
                           obj.fn = expression(cor(predicted, actual, method = "spearman") / sum((predicted - actual)^2)),
                           objective = "max",
                           linear.approximation = TRUE,
                           lin.only = FALSE,
                           pred.int = NULL,
                           print.trace = TRUE){
  
  if(any(class(variable)%in%c("tbl","data.table"))) variable <- as.vector(unlist(variable))
  
  if(sum(is.na(variable)) > 0) stop("You have some missing values, please address.")
  
  n <- length(variable)
  
  if(is.null(obj.fn)){ stop("Please provide an objective function")}
  objective <- tolower(objective)
 
  if(is.null(training.set) && is.null(h)) stop("Please use the length of the variable less the desired forecast period as the training.set value, or provide a value for h.")
  
  if(min(variable)<0) negative.values <- TRUE
  
  variable <- as.numeric(variable)
  
  if(!is.null(h) && h > 0) h_oos <- h_is <- h else {
    h <- 1
    h_is <- length(variable) - training.set
    h_oos <- NULL
  }
  
  if(h_is < .1 * n) h_eval <- 2*h_is else h_eval <- h_is
  
  actual <- tail(variable, h_eval)
  
  if(is.null(training.set)){
    l <- n - h_eval 
    training.set <- l
  } else l <- training.set
  
  if(l <= .5 * n) stop("Please provide a 'training.set' value (integer) less than (2*h) or a smaller (h).")
  if(training.set == n) stop("Please provide a 'training.set' value (integer) less than the length of the variable.")
  
  denominator <- min(4, max(3, ifelse((l/100)%%1 < .5, floor(l/100), ceiling(l/100))))
  
  seasonal.factor <- seasonal.factor[seasonal.factor <= (l/denominator)]
  seasonal.factor <- unique(seasonal.factor)
  
  if(length(seasonal.factor)==0) stop(paste0('Please ensure "seasonal.factor" contains elements less than ', l/denominator, ", otherwise use cross-validation of seasonal factors as demonstrated in the vignette >>> Getting Started with NNS: Forecasting"))
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  seasonal.combs <- nns.estimates <- vector(mode = "list")
 
  previous.seasonals <- previous.estimates <- overall.estimates <- overall.seasonals <- vector(mode = "list")

  
  if(lin.only) methods <- "lin" else methods <- c("lin", "nonlin", "both")

  for(j in methods){
    seasonal.combs <- current.seasonals <- vector(mode = "list")
    current.estimate <- numeric()

    
    for(i in 1 : length(seasonal.factor)){
      if(i == 1){
        seasonal.combs[[i]] <- t(seasonal.factor)
      } else {
        remaining.index <- !(seasonal.factor%in%current.seasonals[[i-1]])
        if(sum(remaining.index)==0){ break }
        seasonal.combs[[i]] <- rbind(replicate(length(seasonal.factor[remaining.index]), current.seasonals[[i-1]]), as.integer(seasonal.factor[remaining.index]))
      }
      
      if(i == 1){
        if(linear.approximation  && j!="lin"){
          seasonal.combs[[1]] <- matrix(unlist(overall.seasonals[[1]]), ncol=1)
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
      
      if(is.null(ncol(seasonal.combs[[i]])) || dim(seasonal.combs[[i]])[2]==0) break 
      
      if(j!="lin" && linear.approximation){
        # Find the min (obj.fn) for a given seasonals sequence
        actual <- tail(variable, h_eval)
        
        predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = unlist(overall.seasonals[[1]]), method = j, plot = FALSE, negative.values = negative.values)
        
        nns.estimates.indiv <- eval(obj.fn)
      } else {
        nns.estimates.indiv <- lapply(1 : ncol(seasonal.combs[[i]]), function(k) {
          actual <- tail(variable, h_eval)
          
          predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor =  seasonal.combs[[i]][ , k], method = j, plot = FALSE, negative.values = negative.values)
          
          return(eval(obj.fn))
        })
      }

      nns.estimates.indiv <- unlist(nns.estimates.indiv)
     
      if(objective=='min') nns.estimates.indiv[is.na(nns.estimates.indiv)] <- Inf else nns.estimates.indiv[is.na(nns.estimates.indiv)] <- -Inf
      
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
    
    
    if(lin.only) predicted <- current.estimate
    if(j!='lin' && lin.only) break
    
  } # for j in c("lin", "nonlin", "both")
  

  if(objective == "min"){
    nns.periods <- unlist(overall.seasonals[[which.min(unlist(overall.estimates))]])
    nns.method <- c("lin","nonlin","both")[which.min(unlist(overall.estimates))]
    nns.SSE <- min(unlist(overall.estimates))
    
    if(length(nns.periods)>1){
      predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = rep((1/length(nns.periods)),length(nns.periods)))
      
      weight.SSE <- eval(obj.fn)
      
      if(weight.SSE < nns.SSE){
        nns.weights <- rep((1/length(nns.periods)),length(nns.periods))
        
        errors <- predicted - actual
        bias <- gravity(na.omit(errors))
        if(is.na(bias)) bias <- 0
        predicted <- predicted - bias
        bias.SSE <- eval(obj.fn)
        
        if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE > weight.SSE) bias <- 0
      } else {
        nns.weights <- NULL
        
        errors <- predicted - actual
        bias <- gravity(na.omit(errors))
        if(is.na(bias)) bias <- 0
        predicted <- predicted - bias
        bias.SSE <- eval(obj.fn)
        
        if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE >= nns.SSE) bias <- 0
      }
    } else {
      nns.weights <- NULL
      
      errors <- predicted - actual
      bias <- gravity(na.omit(errors))
      if(is.na(bias)) bias <- 0
      predicted <- predicted - bias
      bias.SSE <- eval(obj.fn)
      
      if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE >= nns.SSE) bias <- 0
    }
    
  } else {
    nns.periods <- unlist(overall.seasonals[[which.max(unlist(overall.estimates))]])
    nns.method <- c("lin","nonlin","both")[which.max(unlist(overall.estimates))]
    nns.SSE <- max(unlist(overall.estimates))
    
    if(length(nns.periods) > 1){
      predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = rep((1/length(nns.periods)),length(nns.periods)))
      
      weight.SSE <- eval(obj.fn)
      
      if(weight.SSE > nns.SSE){
        nns.weights <- rep((1/length(nns.periods)),length(nns.periods))
        
        errors <- predicted - actual
        bias <- gravity(na.omit(errors))
        if(is.na(bias)) bias <- 0
        predicted <- predicted - bias
        bias.SSE <- eval(obj.fn)
        
        if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE <= weight.SSE) bias <- 0
        
      } else {
        nns.weights <- NULL
        
        errors <- predicted - actual
        bias <- gravity(na.omit(errors))
        if(is.na(bias)) bias <- 0
        predicted <- predicted - bias
        bias.SSE <- eval(obj.fn)
        
        if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE <= nns.SSE) bias <- 0
      }
    } else {
      nns.weights <- NULL
      predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values)
      
      errors <- predicted - actual
      bias <- gravity(na.omit(errors))
      if(is.na(bias)) bias <- 0
      predicted <- predicted - bias
      bias.SSE <- eval(obj.fn)
      if(objective=="min"){
        if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE >= nns.SSE) bias <- 0
      } else {
        if(is.na(bias.SSE)) bias <- 0 else if(bias.SSE <= nns.SSE) bias <- 0
      }
    }
  }
  
  predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = TRUE)
  
  if(objective == "min"){
    if(eval(obj.fn) < nns.SSE) nns.shrink = TRUE else nns.shrink = FALSE
  }
  
  if(objective == "max"){
    if(eval(obj.fn) > nns.SSE) nns.shrink = TRUE else nns.shrink = FALSE
  }
  
  
  regressed_variable <- NNS.reg(1:length(variable), variable, plot = FALSE)$Fitted.xy$y.hat
  
  predicted <- NNS.ARMA(regressed_variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = TRUE)
  
  nns.regress <- FALSE
  
  if(objective == "min"){
    if(eval(obj.fn) < nns.SSE){
      variable <- regressed_variable
      nns.regress <- TRUE
    }
  }
  
  if(objective == "max"){
    if(eval(obj.fn) > nns.SSE){
      variable <- regressed_variable
      nns.regress <- TRUE
    }
  }
  
  
  options(warn = oldw)
  
  
  if(is.null(h_oos)){
    model.results <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = nns.shrink) - bias
  } else {
    model.results <- NNS.ARMA(variable, h = h_oos, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = nns.shrink) - bias
  }
  
  if(!is.null(pred.int)){
      lower_PIs <- model.results - abs(LPM.VaR((1-pred.int)/2, 0, errors)) - abs(bias)
      upper_PIs <- model.results + abs(UPM.VaR((1-pred.int)/2, 0, errors)) + abs(bias)
  } else {
      upper_PIs <- lower_PIs <- NULL
  } 
  
  if(!negative.values){
      model.results <- pmax(0, model.results)
      lower_PIs <- pmax(0, lower_PIs)
      upper_PIs <- pmax(0, upper_PIs)
  }
  

  
  
  return(list(periods = nns.periods,
              weights = nns.weights,
              obj.fn = nns.SSE,
              method = nns.method,
              shrink = nns.shrink,
              nns.regress = nns.regress,
              bias.shift = -bias,
              errors = errors,
              results = model.results,
              lower.pred.int = lower_PIs,
              upper.pred.int = upper_PIs))
}