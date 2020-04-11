#' NNS VAR
#'
#' Nonparametric vector autoregressive model incorporating \link{NNS.ARMA} estimates of variables into \link{NNS.reg} for a multi-variate time-series forecast.
#'
#' @param variables a numeric matrix or data.frame of contemporaneous time-series to forecast.
#' @param h integer; 1 (default) Number of periods to forecast.
#' @param tau positive integer [ > 0]; 1 (default) Number of lagged observations to consider for the time-series data.  Vector for single lag for each respective variable or list for multiple lags per each variable.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all") method for reducing regressors via \link{NNS.stack}.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for dimension reduction in the lagged variable matrix.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param obj.fn expression;
#' \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.ARMA.optim}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns the following matrices of forecasted variables:
#' \itemize{
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
#'
#'  ### Using lags 1:4 for each variable
#'  NNS.VAR(A, h = 12, tau = 4, status = TRUE)
#'
#'  ### Using lag 1 for variable 1, lag 3 for variable 2 and lag 3 for variable 3
#'  NNS.VAR(A, h = 12, tau = c(1,3,3), status = TRUE)
#'
#'  ### Using lags c(1,2,3) for variables 1 and 3, while using lags c(4,5,6) for variable 2
#'  NNS.VAR(A, h = 12, tau = list(c(1,2,3), c(4,5,6), c(1,2,3)), status = TRUE)
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
                    ncores = NULL){

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
    na_s <- sum(is.na(variable))
    variable <- na.omit(variable)
    periods <- NNS.seas(variable, modulo = min(tau[[min(i, length(tau))]]),
                        mod.only = FALSE, plot = FALSE)$periods

    b <- NNS.ARMA.optim(variable, seasonal.factor = periods,
                        training.set = length(variable) - 2*(h + na_s),
                        obj.fn = obj.fn,
                        objective = objective,
                        print.trace = status,
                        ncores = 1)

    nns_IVs$results <- NNS.ARMA(variable, h = (h + na_s), seasonal.factor = b$periods, weights = b$weights,
             method = b$method, ncores = 1, plot = FALSE) + b$bias.shift

    nns_IVs$obj_fn <- b$obj.fn

    return(nns_IVs)

  }

  stopCluster(cl)
  registerDoSEQ()

  nns_IVs_results <- data.frame(tail(do.call(cbind, lapply(nns_IVs, `[[`, 1)), h))
  colnames(nns_IVs_results) <- colnames(variables)
  row.names(nns_IVs_results) <- seq_len(h)

  new_values <- list()

  # Combine forecasted IVs onto training data.frame
  for(i in 1:ncol(variables)){
      new_values[[i]] <- c(na.omit(variables[,i]), nns_IVs[[i]]$results  )

  }


  new_values <- data.frame(do.call(cbind, new_values))
  colnames(new_values) <- as.character(colnames(variables))

  # Now lag new forecasted data.frame
  lagged_new_values <- lag.mtx(new_values, tau = tau)

  # Keep original variables as training set
  lagged_new_values_train <- head(lagged_new_values, dim(lagged_new_values)[1] - h)

  # Select tau = 0 as test set DVs
  DVs <- which(grepl("tau.0", colnames(lagged_new_values)))

  nns_DVs <- list()
  DV_obj_fn <- list()
  relevant_vars <- list()

  if(status){
    message("Currently generating multi-variate estimates...", "\r", appendLF = TRUE)
  }

  for(i in DVs){
    index <- which(DVs%in%i)

    if(status){
      message("Variable ", index, " of ", length(DVs), appendLF = TRUE)
    }

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
        rel.1 <- cor(cbind(lagged_new_values_train[, i],lagged_new_values_train[, -i]), method = "spearman")
    }

    if(any(dim.red.method == "nns.dep" | dim.red.method == "all")){
        rel.2 <- NNS.dep(cbind(lagged_new_values_train[, i],lagged_new_values_train[, -i]))$Dependence
    }

    if(any(dim.red.method == "nns.caus" | dim.red.method == "all")){
        rel.3 <- NNS.caus(cbind(lagged_new_values_train[, i],lagged_new_values_train[, -i]))
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

    rel_vars <- rel_vars[rel_vars>cor_threshold$NNS.dim.red.threshold]

    rel_vars <- names(lagged_new_values)[names(lagged_new_values)%in%names(rel_vars)]



    relevant_vars[[index]] <- rel_vars

    if(length(rel_vars)==0){
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
                               order = NULL, stack = FALSE)

        nns_DVs[[index]] <- DV_values$stack

        DV_obj_fn[[index]] <- sum( (c(DV_values$OBJfn.reg, DV_values$OBJfn.dim.red) / (DV_values$OBJfn.reg + DV_values$OBJfn.dim.red)) * c(DV_values$OBJfn.reg, DV_values$OBJfn.dim.red))

        } else {
            nns_DVs[[index]] <- nns_IVs_results[,index]

            DV_obj_fn[[index]] <- cor_threshold$OBJfn.dim.red
        }



  }

  nns_DVs <- data.frame(do.call(cbind, nns_DVs))
  colnames(nns_DVs) <- colnames(variables)
  row.names(nns_DVs) <- seq_len(h)

  RV <- lapply(relevant_vars, function(x) if(is.null(x)){NA} else {x})

  RV <- do.call(cbind, lapply(RV, `length<-`, max(lengths(RV))))
  colnames(RV) <- as.character(colnames(variables))

  uni <- numeric()
  multi <- numeric()

  for(i in 1:length(colnames(RV))){
      uni[i] <-  .5 + .5*((sum(unlist(strsplit(colnames(RV)[i], split = "_tau"))[1]==do.call(rbind,(strsplit(na.omit(RV[,i]), split = "_tau")))[,1]) -
                             sum(unlist(strsplit(colnames(RV)[i], split = "_tau"))[1]!=do.call(rbind,(strsplit(na.omit(RV[,i]), split = "_tau")))[,1]))
                          / max(ifelse(length(tau)>1, length(tau[[min(i, length(tau))]]), tau), length(na.omit(RV[,i]))))
      multi[i] <- 1 - uni[i]
  }



  forecasts <- data.frame(Reduce(`+`,list(t(t(nns_IVs_results)*uni) , t(t(nns_DVs)*multi))))
  colnames(forecasts) <- colnames(variables)
  row.names(forecasts) <- seq_len(h)

  return( list("relevant_variables" = RV,
               univariate = nns_IVs_results,
               multivariate = nns_DVs,
               ensemble = forecasts) )
}
