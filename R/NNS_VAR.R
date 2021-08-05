#' NNS VAR
#'
#' Nonparametric vector autoregressive model incorporating \link{NNS.ARMA} estimates of variables into \link{NNS.reg} for a multi-variate time-series forecast.
#'
#' @param variables a numeric matrix or data.frame of contemporaneous time-series to forecast.
#' @param h integer; 1 (default) Number of periods to forecast. \code{(h = 0)} will return just the interpolated and extrapolated values.
#' @param tau positive integer [ > 0]; 1 (default) Number of lagged observations to consider for the time-series data.  Vector for single lag for each respective variable or list for multiple lags per each variable.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all") method for reducing regressors via \link{NNS.stack}.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for dimension reduction in the lagged variable matrix.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param obj.fn expression;
#' \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.ARMA.optim}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param nowcast logical; \code{FALSE} (default) internal call for \link{NNS.nowcast}.
#'
#' @return Returns the following matrices of forecasted variables:
#' \itemize{
#'  \item{\code{"interpolated_and_extrapolated"}} Returns a \code{data.frame} of the linear interpolated and \link{NNS.ARMA} extrapolated values to replace \code{NA} values in the original \code{variables} argument.  This is required for working with variables containing different frequencies, e.g. where \code{NA} would be reported for intra-quarterly data when indexed with monthly periods.
#'  \item{\code{"relevant_variables"}} Returns the relevant variables from the dimension reduction step.
#'
#'  \item{\code{"univariate"}} Returns the univariate \link{NNS.ARMA} forecasts.
#'
#'  \item{\code{"multivariate"}} Returns the multi-variate \link{NNS.reg} forecasts.
#'
#'  \item{\code{"ensemble"}} Returns the ensemble of both \code{"univariate"} and \code{"multivariate"} forecasts.
#'  }
#'
#' @note
#' \itemize{
#' \item \code{dim.red.method = "cor"} is significantly faster than the other methods, but comes at the expense of ignoring possible nonlinear relationships between lagged variables.
#' \item Not recommended for factor variables, even after transformed to numeric.  \link{NNS.reg} is better suited for factor or binary regressor extrapolation.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric Vector Autoregression Using NNS"
#' \url{https://www.ssrn.com/abstract=3489550}
#'
#' Viole, F. (2020) "NOWCASTING with NNS"
#' \url{https://www.ssrn.com/abstract=3586658}
#'
#' Viole, F. (2019) "Forecasting Using NNS"
#' \url{https://www.ssrn.com/abstract=3382300}
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
#' \url{https://link.springer.com/article/10.1007/s10614-017-9713-5}
#'
#' Vinod, H. and Viole, F. (2018) "Clustering and Curve Fitting by Line Segments"
#' \url{https://www.preprints.org/manuscript/201801.0090/v1}
#'
#' @examples
#'
#'  \dontrun{
#'  ####################################################
#'  ### Standard Nonparametric Vector Autoregression ###
#'  ####################################################
#'
#'  set.seed(123)
#'  x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#'  A <- cbind(x = x, y = y, z = z)
#'
#'  ### Using lags 1:4 for each variable
#'  NNS.VAR(A, h = 12, tau = 4, status = TRUE)
#'
#'  ### Using lag 1 for variable 1, lag 3 for variable 2 and lag 3 for variable 3
#'  NNS.VAR(A, h = 12, tau = c(1,3,3), status = TRUE)
#'
#'  ### Using lags c(1,2,3) for variables 1 and 3, while using lags c(4,5,6) for variable 2
#'  NNS.VAR(A, h = 12, tau = list(c(1,2,3), c(4,5,6), c(1,2,3)), status = TRUE)
#'
#'  ### CONFIDENCE INTERVALS FOR PREDICTIONS
#'  # Store NNS.VAR output
#'  nns_estimate <- NNS.VAR(A, h = 12, tau = 4, status = TRUE)
#'
#'  # Create bootstrap replicates using NNS.meboot
#'  replicates <- NNS.meboot(nns_estimate$ensemble[,1])$replicates
#'
#'  # Apply UPM.VaR and LPM.VaR for desired confidence interval
#'  # Tail percentage used in first argument per {LPM.VaR} and {UPM.VaR} functions
#'  upper_CIs <- apply(replicates, 1, function(g) UPM.VaR(.025, 0, g))
#'  lower_CIs <- apply(replicates, 1, function(g) LPM.VaR(.025, 0, g))
#'
#'  # View results
#'  cbind(nns_estimate$ensemble[,1], lower_CIs, upper_CIs)
#'
#'
#'  #########################################
#'  ### NOWCASTING with Mixed Frequencies ###
#'  #########################################
#'
#'  library(Quandl)
#'  econ_variables <- Quandl(c("FRED/GDPC1", "FRED/UNRATE", "FRED/CPIAUCSL"),type = 'ts',
#'                           order = "asc", collapse = "monthly", start_date="2000-01-01")
#'
#'  ### Note the missing values that need to be imputed
#'  head(econ_variables)
#'  tail(econ_variables)
#'
#'
#'  NNS.VAR(econ_variables, h = 12, tau = 12, status = TRUE)
#'  }
#'
#' @export



NNS.VAR <- function(variables,
                    h,
                    tau = 1,
                    dim.red.method = "cor",
                    obj.fn = expression( sum((predicted - actual)^2) ),
                    objective = "min",
                    status = TRUE,
                    ncores = NULL,
                    nowcast = FALSE){

  oldw <- getOption("warn")
  options(warn = -1)

  if(nowcast) dates <- zoo::as.yearmon(zoo::as.yearmon(rownames(variables)[1]) + seq(0, (dim(variables)[1] + (h-1)))/12)

  if(any(class(variables)=="tbl")) variables <- as.data.frame(variables)

  dim.red.method <- tolower(dim.red.method)
  if(sum(dim.red.method%in%c("cor","nns.dep","nns.caus","all"))==0){ stop('Please ensure the dimension reduction method is set to one of "cor", "nns.dep", "nns.caus" or "all".')}

  nns_IVs <- list()


  if(is.null(colnames(variables))){
    var_names <- character()
    for(i in 1:ncol(variables)){
      var_names[i] <- paste0("x",i)
    }
    colnames(variables) <- var_names
  }

  colnames(variables) <- gsub(" - ", "...", colnames(variables))

  # Parallel process...
  if (is.null(ncores)) {
    num_cores <- as.integer(parallel::detectCores()) - 1
  } else {
    num_cores <- ncores
  }

  if(num_cores>1) doParallel::registerDoParallel(cores = num_cores)

  if(status) message("Currently generating univariate estimates...","\r", appendLF=TRUE)

  na_s <- numeric()


  nns_IVs <- foreach(i = 1:ncol(variables), .packages = c("NNS", "data.table"))%dopar%{
    # For Interpolation / Extrapolation of all missing values
    index <- seq_len(dim(variables)[1])
    last_point <- tail(index, 1)
    a <- cbind.data.frame("index" = index, variables)
    a <- a[, c(1,(i+1))]
    interpolation_start <- which(!is.na(a[,2]))[1]
    interpolation_point <- tail(which(!is.na(a[,2])), 1)
    a <- a[complete.cases(a),]

    if(dim(a)[1]<last_point){
      nns_IVs$interpolation <- NNS.reg(a[,1], a[,2], order = "max",
                                       point.est = index, plot = FALSE,
                                       ncores = 1, point.only = TRUE)$Point.est

      new_variable <- nns_IVs$interpolation
    } else {
      new_variable <- variables[,i]
      nns_IVs$interpolation <- new_variable
    }

    na_s[i] <- tail(index, 1) - interpolation_point
    if(na_s[i] > 0){
      periods <- NNS.seas(new_variable, modulo = min(tau[[min(i, length(tau))]]),
                          mod.only = FALSE, plot = FALSE)$periods

      ts <- interpolation_point - 2*(h + na_s[i])
      if(ts < 100) ts <- interpolation_point - (h + na_s[i])

      b <- NNS.ARMA.optim(new_variable, seasonal.factor = periods,
                          training.set = ts,
                          obj.fn = obj.fn,
                          objective = objective,
                          print.trace = status,
                          ncores = 1)

      nns_IVs$results <- NNS.ARMA(new_variable, h = (h + na_s[i]), seasonal.factor = b$periods, weights = b$weights,
                                  method = b$method, ncores = 1, plot = FALSE) + b$bias.shift



      na_s_extrapolation <- rowMeans(cbind(tail(nns_IVs$interpolation, na_s[i]), head(nns_IVs$results, na_s[i])))
      nns_IVs$interpolation <- c(nns_IVs$interpolation, na_s_extrapolation)

      nns_IVs$obj_fn <- b$obj.fn
    } else {
      nns_IVs$results <- nns_IVs$interpolation
    }

    return(list(nns_IVs, na.omit(na_s[i]), head(nns_IVs$results, na_s[i])))
  }


  univariate_extrapolation <- lapply(nns_IVs, `[[`, 3)
  na_s <- unlist(lapply(nns_IVs, `[[`, 2))
  nns_IVs <- lapply(nns_IVs, `[[`, 1)

  nns_IVs_interpolated_extrapolated <- data.frame(do.call(cbind, lapply(lapply(nns_IVs, `[[`, 1), function(x) head(x, dim(variables)[1]))))

  nns_IVs_results <- data.frame(do.call(cbind, lapply(lapply(nns_IVs, `[[`, 2), function(x) tail(x, h))))
  colnames(nns_IVs_results) <- colnames(variables)


  if(h == 0){
    colnames(nns_IVs_interpolated_extrapolated) <- as.character(colnames(variables))
    return(nns_IVs_interpolated_extrapolated)
  }


  new_values <- list()

  # Combine interpolated / extrapolated / forecasted IVs onto training data.frame
  for(i in 1:ncol(variables)){
    new_values[[i]] <- c(nns_IVs_interpolated_extrapolated[,i], nns_IVs_results[,i])
  }



  new_values <- data.frame(do.call(cbind, new_values))
  colnames(new_values) <- as.character(colnames(variables))

  nns_IVs_interpolated_extrapolated <- head(new_values, dim(variables)[1])

  # Now lag new forecasted data.frame
  lagged_new_values <- lag.mtx(new_values, tau = tau)

  # Keep original variables as training set
  lagged_new_values_train <- head(lagged_new_values, dim(lagged_new_values)[1] - h)


  if(status) message("Currently generating multi-variate estimates...", "\r", appendLF = TRUE)


  if(num_cores>1){
    if(status) message("Parallel process running, status unavailable... \n","\r",appendLF=FALSE)
    status <- FALSE
  }

  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }

  nns_DVs <- list()
  relevant_vars <- list()

  lists <- foreach(i = 1:ncol(variables), .packages = c("NNS", "data.table"), .combine = 'comb', .init = list(list(), list()),
                   .multicombine = TRUE)%dopar%{


                     if(status) message("Variable ", i, " of ", ncol(variables), appendLF = TRUE)

                     # Dimension reduction NNS.reg to reduce variables
                     cor_threshold <- NNS.stack(IVs.train = lagged_new_values_train[, -i],
                                                DV.train = lagged_new_values_train[, i],
                                                ts.test = 2*h, folds = 1,
                                                obj.fn = obj.fn,
                                                objective = objective,
                                                method = 2,
                                                dim.red.method = dim.red.method,
                                                order = NULL)

                     if(any(dim.red.method == "cor" | dim.red.method == "all")){
                       rel.1 <- abs(cor(cbind(lagged_new_values_train[, i], lagged_new_values_train[, -i]), method = "spearman"))
                     }

                     if(any(dim.red.method == "nns.dep" | dim.red.method == "all")){
                       rel.2 <- NNS.dep(cbind(lagged_new_values_train[, i], lagged_new_values_train[, -i]))$Dependence
                     }

                     if(any(dim.red.method == "nns.caus" | dim.red.method == "all")){
                       rel.3 <- NNS.caus(cbind(lagged_new_values_train[, i], lagged_new_values_train[, -i]))
                     }

                     if(dim.red.method == "cor"){
                       rel_vars <- rel.1[-1,1]
                     }

                     if(dim.red.method == "nns.dep"){
                       rel_vars <- rel.2[-1,1]
                     }

                     if(dim.red.method == "nns.caus"){
                       rel_vars <- rel.3[1,-1]
                     }

                     if(dim.red.method == "all"){
                       rel_vars <- ((rel.1+rel.2+rel.3)/3)[1, -1]
                     }

                     rel_vars <- names(rel_vars[rel_vars > cor_threshold$NNS.dim.red.threshold])
                     rel_vars <- rel_vars[rel_vars!=i]

                     relevant_vars <- rel_vars

                     if(any(length(rel_vars)==0 | is.null(relevant_vars))){
                       rel_vars <- names(lagged_new_values_train)
                     }

                     # NNS.stack() cross-validates the parameters of the multivariate NNS.reg() and dimension reduction NNS.reg()
                     if(length(rel_vars)>1){
                       DV_values <- NNS.stack(lagged_new_values_train[, rel_vars],
                                              lagged_new_values_train[, i],
                                              IVs.test =  tail(lagged_new_values[, rel_vars], h),
                                              obj.fn = obj.fn,
                                              objective = objective,
                                              ts.test = 2*h, folds = 1,
                                              status = status, ncores = num_cores,
                                              dim.red.method = dim.red.method,
                                              order = NULL, stack = TRUE)


                       nns_DVs <- DV_values$stack
                       nns_DVs[is.na(nns_DVs)] <- nns_IVs_results[is.na(nns_DVs),i]
                     } else {
                       nns_DVs <- nns_IVs_results[,i]
                     }

                     list(nns_DVs, relevant_vars)

                   }

  if(num_cores>1) registerDoSEQ()

  nns_DVs <- lists[[1]]
  relevant_vars <- lists[[2]]

  nns_DVs <- data.frame(do.call(cbind, nns_DVs))
  nns_DVs <- head(nns_DVs, h)

  RV <- lapply(relevant_vars, function(x) if(length(x)==0){NA} else {x})

  colnames(nns_DVs) <- colnames(variables)

  RV <- do.call(cbind, lapply(RV, `length<-`, max(lengths(RV))))
  colnames(RV) <- as.character(colnames(variables))

  uni <- numeric()
  multi <- numeric()

  for(i in 1:length(colnames(RV))){
    if(length(na.omit(RV[,i]) > 0)){
      given_var <- unlist(strsplit(colnames(RV)[i], split = "_tau"))[1]
      observed_var <- do.call(rbind,(strsplit(na.omit(RV[,i]), split = "_tau")))[,1]

      equal_tau <- sum(given_var==observed_var)
      unequal_tau <- sum(given_var!=observed_var)

      if(equal_tau > unequal_tau) uni[i] <-  .5 + .5*((equal_tau)/(equal_tau + unequal_tau)) else uni[i] <-  .5 - .5*((unequal_tau)/(equal_tau + unequal_tau))

      if(equal_tau == unequal_tau)uni[i] <- 0.5

      multi[i] <- 1 - uni[i]
    } else {
      uni[i] <- 0.5
      multi[i] <- 0.5
    }
  }


  forecasts <- data.frame(Reduce(`+`,list(t(t(nns_IVs_results)*uni) , t(t(nns_DVs)*multi))))
  colnames(forecasts) <- colnames(variables)

  rownames(nns_IVs_interpolated_extrapolated) <- rownames(variables)
  colnames(nns_IVs_interpolated_extrapolated) <- colnames(variables)

  colnames(nns_IVs_results) <- colnames(variables)
  rownames(nns_IVs_results) <- tail(dates, h)
  colnames(nns_DVs) <- colnames(variables)
  rownames(nns_DVs) <- tail(dates, h)
  colnames(forecasts) <- colnames(variables)
  rownames(forecasts) <- tail(dates, h)

  options(warn = oldw)

  if(sum(na_s)>0){
    return( list("interpolated_and_extrapolated" = nns_IVs_interpolated_extrapolated,
                 "relevant_variables" = data.frame(RV),
                 univariate = nns_IVs_results,
                 multivariate = nns_DVs,
                 ensemble = forecasts) )
  } else {
    return( list("relevant_variables" = data.frame(RV),
                 univariate = nns_IVs_results,
                 multivariate = nns_DVs,
                 ensemble = forecasts) )

  }
}
