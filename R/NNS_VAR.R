#' NNS VAR
#'
#' Nonparametric vector autoregressive model incorporating \link{NNS.ARMA} estimates of variables into \link{NNS.reg} for a multi-variate time-series forecast.
#'
#' @param variables a numeric matrix or data.frame of contemporaneous time-series to forecast.
#' @param h integer; 1 (default) Number of periods to forecast.
#' @param tau integer; 0 (default) Number of lagged observations to consider for the time-series data.
#' @param obj.fn expression;
#' \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param epochs integer; \code{100} (default) Total number of feature combinations to run.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.ARMA.optim}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns the following matrices of forecasted variables:
#' \itemize{
#'  \item{\code{"univariate"}} Returns the univariate \link{NNS.ARMA} forecasts.
#'
#'  \item{\code{"multivariate"}} Returns the multi-variate \link{NNS.reg} forecasts.
#'
#'  \item{\code{"ensemble"}} Returns the ensemble of both \code{"univariate"} and \code{"multivariate"} forecasts.
#'  }
#'
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
#'
#' Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric Vector Autoregression Using NNS"
#' \url{https://ssrn.com/abstract=3489550}
#'
#' Viole, F. (2019) "Forecasting Using NNS"
#' \url{https://ssrn.com/abstract=3382300}
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
#'  set.seed(123)
#'  x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#'  A <- cbind(x = x, y = y, z = z)
#'  NNS.VAR(A, h = 12, tau = 4, status = TRUE)
#'  }
#'
#' @export



NNS.VAR <- function(variables,
                    h,
                    tau = 0,
                    obj.fn = expression( sum((predicted - actual)^2) ),
                    objective = "min",
                    epochs = 100,
                    status = TRUE,
                    ncores = NULL){

  nns_IVs <- list()

  # Parallel process...
  if (is.null(ncores)) {
    cores <- detectCores()
    num_cores <- cores - 1
  } else {
    cores <- detectCores()
    num_cores <- ncores
  }

  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)

  if(status){
    message("Currently generating univariate estimates...","\r", appendLF=TRUE)
  }

  nns_IVs <- foreach(i = 1:ncol(variables), .packages = 'NNS')%dopar%{
    variable <- variables[, i]

    periods <- NNS.seas(variable, modulo = tau,
                        mod.only = FALSE, plot = FALSE)$periods

    b <- NNS.ARMA.optim(variable, seasonal.factor = periods,
                        training.set = length(variable) - 2*h,
                        obj.fn = obj.fn,
                        objective = objective,
                        print.trace = status,
                        ncores = 1)

    nns_IVs$results <- NNS.ARMA(variable, h = h, seasonal.factor = b$periods, weights = b$weights,
             method = b$method, ncores = 1, plot = FALSE) + b$bias.shift

    nns_IVs$obj_fn <- b$obj.fn

    return(nns_IVs)

  }

  stopCluster(cl)
  registerDoSEQ()

  nns_IVs_results <- do.call(cbind, lapply(nns_IVs, `[[`, 1))
  colnames(nns_IVs_results) <- colnames(variables)

  # Combine forecasted IVs onto training data.frame
  new_values <- rbind(variables, nns_IVs_results)

  # Now lag new forecasted data.frame
  lagged_new_values <- lag.mtx(new_values, tau = tau)

  # Keep original variables as training set
  lagged_new_values_train <- head(lagged_new_values, dim(lagged_new_values)[1] - h)

  # Select tau = 0 as test set DVs
  DVs <- which(grepl("tau.0", colnames(lagged_new_values)))

  nns_DVs <- list()
  DV_obj_fn <- list()

  if(status){
    message("Currently generating multi-variate estimates...", "\r", appendLF = TRUE)
  }

  for(i in DVs){
    index <- which(DVs%in%i)
    if(status){
      message("Variable ", index, " of ", length(DVs), appendLF = TRUE)
    }

# NNS.boost() is an ensemble method comparable to xgboost, and aids in dimension reduction
    nns_boost_est <- NNS.boost(lagged_new_values_train[, -i], lagged_new_values_train[, i],
                               IVs.test = tail(lagged_new_values_train[, -i], h),
                               obj.fn = obj.fn,
                               objective = objective,
                               ts.test = 2*h, folds = 1,
                               depth = "max",
                               learner.trials = epochs,
                               ncores = num_cores, type = NULL,
                               feature.importance = FALSE)

# NNS.stack() cross-validates the parameters of the multivariate NNS.reg() and dimension reduction NNS.reg()
    relevant_vars <- colnames(lagged_new_values)%in%names(nns_boost_est$feature.weights)

    DV_values <- NNS.stack(lagged_new_values_train[, relevant_vars],
                                  lagged_new_values_train[, i],
                                  IVs.test =  tail(lagged_new_values[, relevant_vars], h),
                                  obj.fn = obj.fn,
                                  objective = objective,
                                  order = "max",
                                  ts.test = 2*h, folds = 1,
                                  status = status, ncores = num_cores)

    nns_DVs[[index]] <- DV_values$stack

    DV_obj_fn[[index]] <- sum( (c(DV_values$OBJfn.reg, DV_values$OBJfn.dim.red) / (DV_values$OBJfn.reg + DV_values$OBJfn.dim.red)) * c(DV_values$OBJfn.reg, DV_values$OBJfn.dim.red))

  }

  nns_DVs <- do.call(cbind, nns_DVs)
  colnames(nns_DVs) <- colnames(variables)

  if(objective=="min"){
      IV_weights <- 1/unlist(lapply(nns_IVs, `[[`, 2))
      DV_weights <- 1/unlist(DV_obj_fn)

  } else {
      IV_weights <- unlist(lapply(nns_IVs, `[[`, 2))
      DV_weights <- unlist(DV_obj_fn)
  }

  denom <- (IV_weights + DV_weights)
  IV_weights <- IV_weights / denom
  DV_weights <- DV_weights / denom

  IV_weights <- rep(IV_weights, each = dim(nns_IVs_results)[1])
  DV_weights <- rep(DV_weights, each = dim(nns_DVs)[1])

  forecasts <- (IV_weights * nns_IVs_results + DV_weights * nns_DVs)
  colnames(forecasts) <- colnames(variables)

  return( list(univariate = nns_IVs_results, multivariate = nns_DVs, ensemble = forecasts) )

}
