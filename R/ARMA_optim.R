#' NNS ARMA Optimizer
#'
#' Wrapper function for optimizing any combination of a given \code{seasonal.factor} vector in \link{NNS.ARMA}.  The objective function \code{obj.fn} on a validation window (the last \code{length(variable) - training.set} observations) is used to determine the optimum across all \link{NNS.ARMA} methods.
#'
#' @param variable a numeric vector.
#' @param h integer; \code{NULL} (default) Number of periods to forecast out of sample.  If \code{NULL}, \code{h = length(variable) - training.set}.
#' @param training.set integer; \code{NULL} (default) Sets the number of variable observations as the training set.  See \code{Note} below for recommended uses.
#' @param seasonal.factor integers; Multiple frequency integers considered for \link{NNS.ARMA} model, i.e. \code{(seasonal.factor = c(12, 24, 36))}.
#' @param lin.only logical; \code{FALSE} (default) For fast optimization of the linear regression method.  More robust than \code{lin.only = TRUE}.
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.  It will automatically select \code{(negative.values = TRUE)} if the training data \code{variable[1 : training.set]} has a negative value, and for the out-of-sample forecast if \code{variable} does.
#' @param obj.fn expression;
#' The default is the mean squared error divided by the degree-1 co-movement of \code{predicted} and \code{actual}, \code{mean((predicted - actual)^2) / (Co.LPM(1, predicted, actual, mean(predicted), mean(actual)) + Co.UPM(1, predicted, actual, mean(predicted), mean(actual)))} (\code{Inf} when the co-movement is 0), to be minimized.  Any \code{expression(...)} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param linear.approximation logical; \code{TRUE} (default) Uses the best linear output from \code{NNS.reg} to generate a nonlinear and mixture regression for comparison.  \code{FALSE} is a more exhaustive search over the objective space.
#' @param pred.int numeric [0, 1]; 0.95 (default) Returns the associated prediction intervals for the final estimate.  Constructed from the returned model's internal validation residuals (actual - forecast), adding their \link{LPM.VaR} and \link{UPM.VaR} at \code{(1 - pred.int)/2} to the final estimates.
#' @param print.trace logical; \code{TRUE} (default) Prints current iteration information.  Suggested as backup in case of error, best parameters to that point still known and copyable!
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param plot logical; \code{FALSE} (default)
#' @param bias.correction \code{"split"} (default) estimates a gravity bias from the first half of the validation errors and accepts it only if the objective improves on the later half. At least four validation observations are required. \code{"none"} disables correction; \code{"legacy"} reproduces same-window estimation and acceptance. Model selection still uses the validation window, so these diagnostics are not independent out-of-sample performance estimates.
#'
#' @return Returns a list containing:
#' \itemize{
#' \item{\code{$period}} a vector of optimal seasonal periods
#' \item{\code{$weights}} the optimal weights of each seasonal period between an equal weight or NULL weighting
#' \item{\code{$obj.fn}} the \code{obj.fn} value of the returned model's internal validation forecasts, after \code{$bias.shift}.
#' \item{\code{$method}} the method identifying which \link{NNS.ARMA} method was used.
#' \item{\code{$shrink}} whether to use the \code{shrink} parameter in \link{NNS.ARMA}.
#' \item{\code{$nns.regress}} whether to smooth the variable via \link{NNS.reg} before forecasting.
#' \item{\code{$bias.shift}} the accepted additive correction (0 when rejected or disabled). Add it to the forecast from \link{NNS.ARMA} with the returned model parameters, then apply the nonnegative support clamp if required.
#' \item{\code{$bias.validation}} correction policy, estimated/applied shifts, chronological estimate and acceptance rows, and raw/adjusted acceptance scores.
#' \item{\code{$raw.results}} forecasts before the additive bias correction.
#' \item{\code{$validation.forecast}} the final internal validation forecasts, including correction and support clamp.
#' \item{\code{$errors}} a vector of the returned model's internal validation errors (forecast - actual), before \code{$bias.shift}.
#' \item{\code{$results}} a vector of length \code{h}.
#' \item{\code{$lower.pred.int}} a vector of lower prediction intervals per forecast point.
#' \item{\code{$upper.pred.int}} a vector of upper prediction intervals per forecast point.
#'}
#' @note
#' \itemize{
#' \item{} Typically, \code{(training.set = 0.8 * length(variable))} is used for optimization.  Smaller samples could use \code{(training.set = 0.9 * length(variable))} (or larger) in order to preserve information.
#'
#' \item{} The number of combinations will grow prohibitively large, they should be kept as small as possible.  \code{seasonal.factor} containing an element too large will result in an error.  Please reduce the maximum \code{seasonal.factor}.
#'
#' \item{} Set \code{(ncores = 1)} if routine is used within a parallel architecture.
#'}
#'
#'
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
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
                           lin.only = FALSE,
                           negative.values = FALSE,
                           obj.fn =  expression( {
                             .denom <- NNS::Co.LPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual)) +
                                       NNS::Co.UPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual))
                             if (.denom == 0) Inf else mean((predicted - actual)^2) / .denom
                           } ),
                           objective = "min",
                           linear.approximation = TRUE,
                           ncores = NULL,
                           pred.int = 0.95,
                           print.trace = TRUE,
                           plot = FALSE,
                           bias.correction = c("split", "none", "legacy")){
  
  if(any(class(variable)%in%c("tbl","data.table"))) variable <- as.vector(unlist(variable))
  
  if(anyNA(variable)) stop("You have some missing values, please address.")
  
  n <- length(variable)
  
  bias.correction <- match.arg(bias.correction)
  if(is.null(obj.fn)){ stop("Please provide an objective function")}
  objective <- tolower(objective)
  
  if(is.null(training.set) && is.null(h)) stop("Please use the length of the variable less the desired forecast period as the [training.set] value, or provide a value for [h].")
  
  variable <- as.numeric(variable)
  OV <- variable
  
  if(!is.null(h) && h > 0) h_oos <- h else {
    h <- NULL
    h_oos <- NULL
  }
  
  if(is.null(training.set)) training.set <- floor(.8 * n)
  training.set <- as.integer(training.set)
  
  h_eval <- as.integer(n - training.set)
  
  actual <- tail(variable, h_eval)
  
  if(training.set <= .5 * n) stop("Please provide a larger [training.set] value (integer) or a smaller [h].")
  if(training.set == n) stop("Please provide a [training.set] value (integer) less than the length of the variable.")
  
  # Validation forecasts see only the training window; out-of-sample results are fitted on the whole series.
  if(min(variable[1 : training.set]) < 0) negative.values <- TRUE
  results.negative.values <- negative.values || (!is.null(h_oos) && min(variable) < 0)
  
  denominator <- min(4, max(3, ifelse((training.set/100)%%1 < .5, floor(training.set/100), ceiling(training.set/100))))
  
  seasonal.factor <- seasonal.factor[seasonal.factor <= (training.set/denominator)]
  seasonal.factor <- unique(seasonal.factor)
  
  if(length(seasonal.factor)==0) stop(paste0('Please ensure [seasonal.factor] contains elements less than ', training.set/denominator, ", otherwise use cross-validation of seasonal factors as demonstrated in the vignette >>> Getting Started with NNS: Forecasting"))
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  seasonal.combs <- nns.estimates <- vector(mode = "list")
  
  previous.seasonals <- previous.estimates <- overall.estimates <- overall.seasonals <- vector(mode = "list")
  
  
  methods <- c("lin", "nonlin", "both")
  if(lin.only) methods <- "lin"

  for(j in methods){
    seasonal.combs <- current.seasonals <- vector(mode = "list")
    current.estimate <- numeric()
    
    if (j == "lin") {
      # Determine the number of cores to use
      auto_cores <- is.null(ncores)
      num_cores <- if (auto_cores) {
        max(1L, as.integer(parallel::detectCores()) - 1L, na.rm = TRUE)
      } else {
        max(1L, as.integer(ncores), na.rm = TRUE)
      }

      # Workload gate: standing up a process cluster costs ~1s (spawn +
      # exporting data to workers). For small series / few seasonal candidates
      # that fixed cost can exceed the compute it parallelizes. The series-length
      # / candidate-count heuristic only gates the AUTOMATIC default
      # (ncores = NULL); an explicit ncores is an intentional request and is
      # honored whenever there is more than one seasonal candidate to spread
      # across cores.
      use_parallel <- num_cores > 1 &&
        length(seasonal.factor) >= 2L &&
        (!auto_cores || (length(variable) >= 500L && length(seasonal.factor) >= 8L))

      # Manage cluster creation
      cl <- NULL
      if (use_parallel) {
        forked <- FALSE
        cl <- tryCatch(
          { cl0 <- parallel::makeForkCluster(num_cores); forked <- TRUE; cl0 },
          error = function(e) parallel::makeCluster(num_cores)
        )
        # Fork workers inherit the already-loaded NNS namespace; only the PSOCK
        # fallback needs an explicit library() call on each worker.
        if (!forked) parallel::clusterEvalQ(cl, library(NNS))
      }
    }
    
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
      
      if (j == "lin") {
        # Parallel or sequential computation based on the workload gate
        nns.estimates.indiv <- if (use_parallel) {
          parallel::clusterExport(
            cl,
            varlist = c("variable", "h_eval", "training.set", "seasonal.combs", "i", "obj.fn", "negative.values", "NNS.ARMA", "print.trace"),
            envir = environment()
          )
          parallel::parLapply(cl, 1:ncol(seasonal.combs[[i]]), function(k) {
            actual <- tail(variable, h_eval)
            predicted <- NNS.ARMA(
              variable,
              training.set = training.set,
              h = h_eval,
              seasonal.factor = seasonal.combs[[i]][, k],
              method = "lin",
              plot = FALSE,
              negative.values = negative.values
            )
            eval(obj.fn)
          })
        } else {
          lapply(1:ncol(seasonal.combs[[i]]), function(k) {
            actual <- tail(variable, h_eval)
            predicted <- NNS.ARMA(
              variable,
              training.set = training.set,
              h = h_eval,
              seasonal.factor = seasonal.combs[[i]][, k],
              method = "lin",
              plot = FALSE,
              negative.values = negative.values
            )
            eval(obj.fn)
          })
        }
        
        # Ensure output is unlisted
        nns.estimates.indiv <- unlist(nns.estimates.indiv)
      }
      
      if(j != "lin" && !linear.approximation){
        # Exhaustive search: score every candidate seasonal set with the method itself.
        nns.estimates.indiv <- unlist(lapply(1:ncol(seasonal.combs[[i]]), function(k) {
          actual <- tail(variable, h_eval)
          predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval,
                                seasonal.factor = seasonal.combs[[i]][, k], method = j,
                                plot = FALSE, negative.values = negative.values)
          eval(obj.fn)
        }))
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
        
        # Score the returned model: method = "both" averages lin and nonlin inside
        # the recursion, unlike the mean of two separate forecasts.
        predicted <- both.predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = unlist(overall.seasonals[[1]]), method = "both", plot = FALSE, negative.values = negative.values)
        
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
    
    if (j == "lin") {
      # Clean up cluster
      if (!is.null(cl)) {
        parallel::stopCluster(cl)
        invisible(gc(verbose = FALSE))  # Clean up memory
      }
    }
  
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
    if(lin.only) nns.method <- "lin" else nns.method <- c("lin","nonlin","both")[which.min(unlist(overall.estimates))]
    nns.SSE <- min(unlist(overall.estimates))
    predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = NULL)
    
    if(length(nns.periods) > 1){
      # Candidate: equal weights across the selected periods.
      predicted.null <- predicted
      predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = rep(1/length(nns.periods), length(nns.periods)))
      weight.SSE <- eval(obj.fn)
      if(is.na(weight.SSE)) weight.SSE <- Inf
      predicted <- predicted.null
      
      if(weight.SSE < nns.SSE){
        nns.weights <- rep((1/length(nns.periods)),length(nns.periods))
        nns.SSE <- weight.SSE
        predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights)
      } else nns.weights <- NULL
    } else nns.weights <- NULL
  } else {
    nns.periods <- unlist(overall.seasonals[[which.max(unlist(overall.estimates))]])
    if(lin.only) nns.method <- "lin" else nns.method <- c("lin","nonlin","both")[which.max(unlist(overall.estimates))]
    nns.SSE <- max(unlist(overall.estimates))
    predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = NULL)
    
    if(length(nns.periods) > 1){
      # Candidate: equal weights across the selected periods.
      predicted.null <- predicted
      predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = rep(1/length(nns.periods), length(nns.periods)))
      weight.SSE <- eval(obj.fn)
      if(is.na(weight.SSE)) weight.SSE <- -Inf
      predicted <- predicted.null
      
      if(weight.SSE > nns.SSE){
        nns.weights <- rep((1/length(nns.periods)),length(nns.periods))
        nns.SSE <- weight.SSE
        predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights)
      } else nns.weights <- NULL
    } else nns.weights <- NULL
  }

  final.predicted <- predicted
  best.SSE <- nns.SSE
 
  predicted <- NNS.ARMA(variable, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = TRUE)
  
  if(objective == "min"){
    if(eval(obj.fn) < nns.SSE){
      nns.shrink = TRUE
      final.predicted <- predicted
      best.SSE <- eval(obj.fn)
    }  else nns.shrink = FALSE
  }
  
  if(objective == "max"){
    if(eval(obj.fn) > nns.SSE){
      nns.shrink = TRUE
      final.predicted <- predicted
      best.SSE <- eval(obj.fn)
    }  else nns.shrink = FALSE
  }
  
  
  # Smooth the training window only, so the validation score never sees the validation data.
  regressed_training <- NNS.reg(1:training.set, variable[1:training.set], plot = FALSE, smooth = TRUE)$Fitted.xy$y.hat
  
  predicted <- NNS.ARMA(regressed_training, training.set = training.set, h = h_eval, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = negative.values, weights = nns.weights, shrink = TRUE)
  
  nns.regress <- FALSE
  
  if(objective == "min"){
    if(eval(obj.fn) < best.SSE){
      nns.regress <- TRUE
      final.predicted <- predicted
    }
  }
  
  if(objective == "max"){
    if(eval(obj.fn) > best.SSE){
      nns.regress <- TRUE
      final.predicted <- predicted
    }
  }

  # The bias parameter has its own chronological estimate/acceptance split.
  # The selected base model still used the internal validation window.
  raw.validation <- final.predicted
  errors <- raw.validation - actual
  bias.validation <- .nns_arma_bias_choice(raw.validation, actual, obj.fn, objective,
                                           negative.values, bias.correction)
  bias <- -bias.validation$applied.shift
  final.predicted <- raw.validation - bias
  if(!negative.values) final.predicted <- pmax(0, final.predicted)
  predicted <- final.predicted
  final.SSE <- eval(obj.fn)

  residuals <- actual - final.predicted
  lower_PIs_is <- final.predicted + LPM.VaR((1-pred.int)/2, 0, residuals)
  upper_PIs_is <- final.predicted + UPM.VaR((1-pred.int)/2, 0, residuals)
  
  options(warn = oldw)
  
  
  # nns.regress was tested on the smoothed series with shrink = TRUE.
  if(nns.regress) nns.shrink <- TRUE
  
  if(is.null(h_oos)){
    if(is.null(h)) h <- h_eval
    raw.results <- raw.validation
    model.results <- final.predicted
  } else {
    if(is.null(h)) h <- h_oos
    final.variable <- if(nns.regress) NNS.reg(1:n, OV, plot = FALSE, smooth = TRUE)$Fitted.xy$y.hat else OV
    raw.results <- NNS.ARMA(final.variable, h = h_oos, seasonal.factor = nns.periods, method = nns.method, plot = FALSE, negative.values = results.negative.values, weights = nns.weights, shrink = nns.shrink)
    model.results <- raw.results - bias
  }
  
  
  # The bands sit on the forecasts the function returns, so clamp those first.
  if(!results.negative.values) model.results <- pmax(0, model.results)
  lower_PIs <- model.results + LPM.VaR((1-pred.int)/2, 0, residuals)
  upper_PIs <- model.results + UPM.VaR((1-pred.int)/2, 0, residuals)

  if(!results.negative.values){
    lower_PIs <- pmax(0, lower_PIs)
    upper_PIs <- pmax(0, upper_PIs)
  }
  if(!negative.values){
    lower_PIs_is <- pmax(0, lower_PIs_is)
    upper_PIs_is <- pmax(0, upper_PIs_is)
  }
  
  if(plot){
    if(is.null(h_oos)) xlim <- c(1, max((training.set + h))) else xlim <- c(1, max((n + h)))
    
    plot(OV, type = 'l', lwd = 2, main = "NNS.ARMA Forecast", col = 'steelblue',
         xlim = xlim,
         ylab =  "Variable",
         ylim = c(min(model.results, variable,  unlist(lower_PIs), unlist(upper_PIs) ), 
                  max(model.results, variable,  unlist(lower_PIs), unlist(upper_PIs) )) )
    
    lfp <- length(final.predicted)
    
    starting.point <- as.integer(n - lfp)
    
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
  
  
  return(.NNS.out(list(periods = nns.periods,
              weights = nns.weights,
              obj.fn = final.SSE,
              method = nns.method,
              shrink = nns.shrink,
              nns.regress = nns.regress,
              bias.shift = -bias,
              bias.validation = bias.validation,
              raw.results = raw.results,
              validation.forecast = final.predicted,
              errors = errors,
              results = model.results,
              lower.pred.int = lower_PIs,
              upper.pred.int = upper_PIs)))
}
# Separate correction estimation from its later acceptance evaluation. This helper
# deliberately does not claim that the base model selection is out of sample.
.nns_arma_bias_choice <- function(predicted, actual, obj.fn, objective,
                                  negative.values, policy) {
  n <- length(actual)
  estimate <- test <- integer()
  candidate <- applied <- 0
  before <- after <- NA_real_
  status <- if (policy == "none") "disabled" else "insufficient_validation"
  if (policy == "legacy") estimate <- test <- seq_len(n)
  if (policy == "split" && n >= 4L) {
    estimate <- seq_len(n %/% 2L)
    test <- seq.int(length(estimate) + 1L, n)
  }
  if (length(estimate) && length(test)) {
    candidate <- -gravity(predicted[estimate] - actual[estimate])
    if (!is.finite(candidate)) candidate <- 0
    score <- function(p) eval(obj.fn,
      envir = list(predicted = p, actual = actual[test]), enclos = parent.frame())
    before <- score(predicted[test])
    adjusted <- predicted[test] + candidate
    if (!negative.values) adjusted <- pmax(0, adjusted)
    after <- score(adjusted)
    accepted <- isTRUE(if (objective == "min") after < before else after > before)
    if (accepted && is.finite(after)) applied <- candidate
    status <- if (applied != 0) "accepted" else "rejected"
  }
  list(policy = policy, estimated.shift = candidate, applied.shift = applied,
       estimate.rows = estimate, acceptance.rows = test, raw.score = before,
       adjusted.score = after, accepted = applied != 0, status = status)
}
