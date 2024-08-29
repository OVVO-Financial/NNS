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
#' @param pred.int numeric [0, 1]; 0.95 (default) Returns the associated prediction intervals for the final estimate.  Constructed using the maximum entropy bootstrap \link{NNS.meboot} on the final estimates.
#' @param print.trace logical; \code{TRUE} (default) Prints current iteration information.  Suggested as backup in case of error, best parameters to that point still known and copyable!
#' @param plot logical; \code{FALSE} (default)
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
#' \item{} Typically, \code{(training.set = 0.8 * length(variable)} is used for optimization.  Smaller samples could use \code{(training.set = 0.9 * length(variable))} (or larger) in order to preserve information.
#'
#' \item{} The number of combinations will grow prohibitively large, they should be kept as small as possible.  \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
#'}
#'
#'
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' 
#' @examples
#'
#' ## Nonlinear NNS.ARMA period optimization using 2 yearly lags on AirPassengers monthly data
#' \dontrun{
#' nns.optims <- NNS.ARMA.optim(AirPassengers[1:132], training.set = 120,
#' seasonal.factor = seq(12, 24, 6))
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
                           obj.fn =  expression( mean((predicted - actual)^2) / (NNS::Co.LPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual)) + NNS::Co.UPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual)) )  ),
                           objective = "min",
                           linear.approximation = TRUE,
                           pred.int = 0.95,
                           print.trace = TRUE,
                           plot = FALSE){
  
  if(any(class(variable)%in%c("tbl","data.table"))) variable <- as.vector(unlist(variable))
  
  if(sum(is.na(variable)) > 0) stop("You have some missing values, please address.")
  
  n <- length(variable)
  
  if(is.null(obj.fn)){ stop("Please provide an objective function")}
  objective <- tolower(objective)
  
  if(is.null(training.set) && is.null(h)) stop("Please use the length of the variable less the desired forecast period as the [training.set] value, or provide a value for [h].")
  
  variable <- as.numeric(variable)
  OV <- variable
  
  if(min(variable) < 0) negative.values <- TRUE
  
  if(!is.null(h) && h > 0) h_oos <- h_is <- h else {
    h <- NULL
    h_oos <- NULL
  }
  
  if(is.null(training.set)) training.set <- .8 * n
  
  h_eval <- h_is <- n - training.set
  
  actual <- tail(variable, h_eval)
  
  if(training.set <= .5 * n) stop("Please provide a larger [training.set] value (integer) or a smaller [h].")
  if(training.set == n) stop("Please provide a [training.set] value (integer) less than the length of the variable.")
  
  denominator <- min(4, max(3, ifelse((training.set/100)%%1 < .5, floor(training.set/100), ceiling(training.set/100))))
  
  seasonal.factor <- seasonal.factor[seasonal.factor <= (training.set/denominator)]
  seasonal.factor <- unique(seasonal.factor)
  
  if(length(seasonal.factor)==0) stop(paste0('Please ensure [seasonal.factor] contains elements less than ', training.set/denominator, ", otherwise use cross-validation of seasonal factors as demonstrated in the vignette >>> Getting Started with NNS: Forecasting"))
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  seasonal.combs <- nns.estimates <- vector(mode = "list")
  
  previous.seasonals <- previous.estimates <- overall.estimates <- overall.seasonals <- vector(mode = "list")
  
  
  methods <- c("lin", "nonlin", "both")
  
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
      
      if(j=="lin"){

        nns.estimates.indiv <- lapply(1 : ncol(seasonal.combs[[i]]), function(k) {
          actual <- tail(variable, h_eval)
          message("Testing seasonal.factor ", paste(unlist(seasonal.combs[[i]][ , k]), ","), "\r", appendLF = FALSE)
          predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor =  seasonal.combs[[i]][ , k], method = "lin", plot = FALSE, negative.values = negative.values)
          
          return(eval(obj.fn))
        })
      }
      
      if(j=="nonlin" && linear.approximation){
        # Find the min (obj.fn) for a given seasonals sequence
        actual <- tail(variable, h_eval)
        
        predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = unlist(overall.seasonals[[1]]), method = j, plot = FALSE, negative.values = negative.values)
        nonlin.predicted <- predicted
        
        nns.estimates.indiv <- eval(obj.fn)
      }
      
      if(j=="both" && linear.approximation){
        # Find the min (obj.fn) for a given seasonals sequence
        actual <- tail(variable, h_eval)
       
        lin.predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = unlist(overall.seasonals[[1]]), method = "lin", plot = FALSE, negative.values = negative.values)
        predicted <- both.predicted <- (lin.predicted + nonlin.predicted) / 2

        nns.estimates.indiv <- eval(obj.fn)
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
      if(which(c("lin","nonlin","both")==j) > 1 ){
        if(sum(as.numeric(unlist(current.seasonals[[i]]))%in%as.numeric(unlist(previous.seasonals[[which(c("lin","nonlin","both")==j)-1]][i])))==length(as.numeric(unlist(current.seasonals[[i]])))){
          
          if(objective=='min'){
            if(current.estimate[i] >= previous.estimates[[which(c("lin","nonlin","both")==j)-1]][i]) break
          } else {
            if(current.estimate[i] <= previous.estimates[[which(c("lin","nonlin","both")==j)-1]][i]) break
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
  
  final.predicted <- predicted
  
  predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = TRUE)
  
  if(objective == "min"){
    if(eval(obj.fn) < nns.SSE){
      nns.shrink = TRUE
      final.predicted <- predicted
    }  else nns.shrink = FALSE
  }
  
  if(objective == "max"){
    if(eval(obj.fn) > nns.SSE){
      nns.shrink = TRUE
      final.predicted <- predicted
    }  else nns.shrink = FALSE
  }
  
  
  regressed_variable <- NNS.reg(1:length(variable), variable, plot = FALSE)$Fitted.xy$y.hat
  
  predicted <- NNS.ARMA(regressed_variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = TRUE)
  
  nns.regress <- FALSE
  
  if(objective == "min"){
    if(eval(obj.fn) < nns.SSE){
      variable <- regressed_variable
      nns.regress <- TRUE
      final.predicted <- predicted
    }
  }
  
  if(objective == "max"){
    if(eval(obj.fn) > nns.SSE){
      variable <- regressed_variable
      nns.regress <- TRUE
      final.predicted <- predicted
    }
  }
  
  lower_PIs_is <- final.predicted - abs(LPM.VaR((1-pred.int)/2, 0, errors)) - abs(bias)
  upper_PIs_is <- final.predicted + abs(UPM.VaR((1-pred.int)/2, 0, errors)) + abs(bias)
  
  options(warn = oldw)
  
  
  if(is.null(h_oos)){
    if(is.null(h)) h <- h_eval
    model.results <- NNS.ARMA(OV, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = nns.shrink) - bias
  } else {
    if(is.null(h)) h <- h_oos
    model.results <- NNS.ARMA(OV, h = h_oos, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = nns.shrink) - bias
  }
  

  lower_PIs <- model.results - abs(LPM.VaR((1-pred.int)/2, 0, errors)) - abs(bias)
  upper_PIs <- model.results + abs(UPM.VaR((1-pred.int)/2, 0, errors)) + abs(bias)
  
  
  if(!negative.values){
    model.results <- pmax(0, model.results)
    lower_PIs <- pmax(0, lower_PIs)
    upper_PIs <- pmax(0, upper_PIs)
  }

  if(plot){
    if(is.null(h_oos)) xlim <- c(1, max((training.set + h))) else xlim <- c(1, max((n + h)))
    
    plot(OV, type = 'l', lwd = 2, main = "NNS.ARMA Forecast", col = 'steelblue',
         xlim = xlim,
         ylab =  "Variable",
         ylim = c(min(model.results, variable,  unlist(lower_PIs), unlist(upper_PIs) ), 
                  max(model.results, variable,  unlist(lower_PIs), unlist(upper_PIs) )) )
    
    lfp <- length(final.predicted)
    
    starting.point <- min(training.set, min(n - lfp))
    
    lines((starting.point + 1) : (starting.point + lfp), final.predicted, col = "red", lwd = 2, lty = 2)
    
    polygon(c((starting.point + 1) : (starting.point + lfp), rev((starting.point + 1) : (starting.point + lfp))),
            c(lower_PIs_is, rev(upper_PIs_is)),
            col = rgb(70/255, 130/255, 180/255, alpha = 0.5),
            border = NA)
    
    lines(OV, lwd = 2, col = "steelblue")
    lines((starting.point + 1) : (starting.point + lfp), final.predicted, col = "red", lwd = 2, lty = 2)
    
    legend("topleft", legend = c("Variable", "Internal Validation"), 
           col = c("steelblue", "red"), lty = c(1, 2), bty = "n", lwd = 2)
    
    if(!is.null(h_oos)){
      lines((n + 1) : (n + h), model.results, col = "red", lwd = 2)
      
      polygon(c((n + 1) : (n + h), rev((n + 1) : (n + h))),
              c(lower_PIs, rev(upper_PIs)),
              col = rgb(1, 192/255, 203/255, alpha = 0.5),
              border = NA)
      
      legend("topleft", legend = c("Variable", "Internal Validation", "Forecast"), 
             col = c("steelblue", "red", "red"), lty = c(1, 2, 1), bty = "n", lwd = 2)
    } 
      
    
    
    
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