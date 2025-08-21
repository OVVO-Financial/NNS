#' NNS VAR
#'
#' Nonparametric vector autoregressive model incorporating \link{NNS.ARMA} estimates of variables into \link{NNS.reg} for a multi-variate time-series forecast.
#'
#' @param variables a numeric matrix or data.frame of contemporaneous time-series to forecast.
#' @param h integer; 1 (default) Number of periods to forecast. \code{(h = 0)} will return just the interpolated and extrapolated values.
#' @param tau positive integer [ > 0]; 1 (default) Number of lagged observations to consider for the time-series data.  Vector for single lag for each respective variable or list for multiple lags per each variable.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all") method for reducing regressors via \link{NNS.stack}.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for dimension reduction in the lagged variable matrix.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param naive.weights logical; \code{TRUE} (default) Equal weights applied to univariate and multivariate outputs in ensemble.  \code{FALSE} will apply weights based on the number of relevant variables detected. 
#' @param obj.fn expression;
#' \code{expression(mean((predicted - actual)^2)) / (Sum of NNS Co-partial moments)} (default) MSE / co-movements is the default objective function.  Any \code{expression(...)} using the specific terms \code{predicted} and \code{actual} can be used.
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
#' \item \code{"Error in { : task xx failed -}"} should be re-run with \code{NNS.VAR(..., ncores = 1)}.
#' \item Not recommended for factor variables, even after transformed to numeric.  \link{NNS.reg} is better suited for factor or binary regressor extrapolation.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric Vector Autoregression Using NNS"  \doi{10.2139/ssrn.3489550}
#'
#' Viole, F. (2020) "NOWCASTING with NNS"  \doi{10.2139/ssrn.3589816}
#'
#' Viole, F. (2019) "Forecasting Using NNS"  \doi{10.2139/ssrn.3382300}
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"  \doi{10.1007/s10614-017-9713-5}
#'
#' Vinod, H. and Viole, F. (2018) "Clustering and Curve Fitting by Line Segments"  \doi{10.20944/preprints201801.0090.v1}
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
#'  ### PREDICTION INTERVALS
#'  # Store NNS.VAR output
#'  nns_estimate <- NNS.VAR(A, h = 12, tau = 4, status = TRUE)
#'
#'  # Create bootstrap replicates using NNS.meboot
#'  replicates <- NNS.meboot(nns_estimate$ensemble[,1], rho = seq(-1,1,.25))["replicates",]
#'  replicates <- do.call(cbind, replicates)
#'
#'  # Apply UPM.VaR and LPM.VaR for desired prediction interval...95 percent illustrated
#'  # Tail percentage used in first argument per {LPM.VaR} and {UPM.VaR} functions
#'  lower_CIs <- apply(replicates, 1, function(z) LPM.VaR(0.025, 0, z))
#'  upper_CIs <- apply(replicates, 1, function(z) UPM.VaR(0.025, 0, z))
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
#'                           order = "asc", collapse = "monthly", start_date = "2000-01-01")
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
                    naive.weights = TRUE,
                    obj.fn = expression( mean((predicted - actual)^2) / (NNS::Co.LPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual)) + NNS::Co.UPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual)) )  ),
                    objective = "min",
                    status = TRUE,
                    ncores = NULL,
                    nowcast = FALSE){
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  dates <- NULL
  
  # --- helpers for data-driven blending ---
  blend_weight_mse_opt <- function(y_val, u_val, m_val, shrink = NULL) {
    d  <- u_val - m_val
    num <- stats::cov(y_val - m_val, d, use = "complete.obs")
    den <- stats::var(d, na.rm = TRUE)
    w <- if (is.finite(num) && is.finite(den) && den > 0) num / den else 0.5
    w <- max(0, min(1, w))               # clip to [0,1]
    if (!is.null(shrink) && is.finite(shrink) && shrink > 0)
      w <- (1 - shrink) * w + shrink * 0.5
    w
  }
  
  pool_logistic <- function(w_tau, w_val, gamma = 0.5, shrink = 0) {
    w_tau <- max(1e-8, min(1 - 1e-8, w_tau))
    w_val <- max(1e-8, min(1 - 1e-8, w_val))
    o_tau <- w_tau / (1 - w_tau)
    o_val <- w_val / (1 - w_val)
    o_star <- (o_tau^(1 - gamma)) * (o_val^gamma)
    w_star <- o_star / (1 + o_star)
    if (shrink > 0) w_star <- (1 - shrink) * w_star + shrink * 0.5
    max(0, min(1, w_star))
  }
  
  # Fully dynamic gamma from validation only (no constants)
  gamma_from_validation <- function(y_val, u_val, m_val, method = c("t", "wilcox")) {
    method <- match.arg(method)
    delta <- (y_val - u_val)^2 - (y_val - m_val)^2  # + => uni better
    if (all(!is.finite(delta)) || length(delta) < 2L) return(0)  # no evidence
    p <- tryCatch({
      if (method == "t") {
        m  <- mean(delta, na.rm = TRUE)
        s  <- stats::sd(delta, na.rm = TRUE)
        n  <- sum(is.finite(delta))
        if (!is.finite(s) || s == 0 || n < 2) 1 else {
          tstat <- m / (s / sqrt(n))
          2 * (1 - stats::pt(abs(tstat), df = n - 1))
        }
      } else {
        stats::wilcox.test(delta, mu = 0, alternative = "two.sided", exact = FALSE)$p.value
      }
    }, error = function(e) 1)
    gamma <- 1 - p
    if (!is.finite(gamma)) gamma <- 0
    max(0, min(1, gamma))
  }
  
  if(nowcast){
    year_mon <- zoo::as.yearmon(format(zoo::index(variables), '%Y-%m'))
    dates <- c(year_mon, tail(year_mon, h) + h/12)
  }
  
  if(any(class(variables)%in%c("tbl","data.table"))) variables <- as.data.frame(variables)
  
  dim.red.method <- tolower(dim.red.method)
  if(sum(dim.red.method%in%c("cor","nns.dep","nns.caus","all"))==0){ stop('Please ensure the dimension reduction method is set to one of "cor", "nns.dep", "nns.caus" or "all".')}
  
  if(is.null(colnames(variables))){
    colnames.list <- lapply(1 : ncol(variables), function(i) paste0("x", i))
    colnames(variables) <- as.character(colnames.list)
  }
  
  if(any(colnames(variables)=="")){
    var_names <- character()
    for(i in 1:length(which(colnames(variables)==""))){
      var_names[i] <- paste0("x",i)
    }
    colnames(variables)[which(colnames(variables)=="")] <- var_names
  }
  
  colnames(variables) <- gsub(" - ", "...", colnames(variables))
  
  # Parallel process...
  if (is.null(ncores)) {
    num_cores <- as.integer(max(2L, parallel::detectCores(), na.rm = TRUE)) - 1
  } else {
    num_cores <- ncores
  }
  
  if(num_cores > 1){
    doParallel::registerDoParallel(num_cores)
    invisible(data.table::setDTthreads(1))
  } else {
    foreach::registerDoSEQ()
    invisible(data.table::setDTthreads(0, throttle = NULL))
  }
  
  if(status) message("Currently generating univariate estimates...","\r", appendLF=TRUE)
  
  nns_IVs <- variable_interpolation <- variable_interpolation_and_extrapolation <- list(ncol(variables))
  
  nns_IVs <- foreach(i = 1:ncol(variables), .packages = c("NNS", "data.table"))%dopar%{
    n <- nrow(variables)
    index <- seq_len(n)
    a <- cbind.data.frame("index" = index, variables)
    
    # For Interpolation / Extrapolation of all missing values
    selected_variable <- a[, c(1,(i+1))]
    
    missing_index <- which(is.na(selected_variable[,2]))
    
    if(length(missing_index)==0){
      variable_interpolation <- selected_variable[,2]
    } else {
      interpolation_point <- tail(which(!is.na(selected_variable[,2])), 1)
      selected_variable <- selected_variable[complete.cases(selected_variable),]
      
      h_int <- tail(index, 1) - interpolation_point
      variable_interpolation <- variables[,i]
      
      if(h_int > 0){
        multi <- NNS.stack(cbind(selected_variable[,1], selected_variable[,1]), selected_variable[,2],
                           order = NULL, ncores = 1, status = FALSE, folds = 5,
                           IVs.test = cbind(missing_index, missing_index), method = 1)$stack
        variable_interpolation[missing_index] <- multi
      } else {
        variable_interpolation <- NNS.reg(selected_variable[,1], selected_variable[,2], order = "max", ncores = 1,
                                          point.est = index, plot = FALSE, point.only = TRUE)$Point.est
      }
    }
    
    if(h > 0){
      periods <- NNS.seas(variable_interpolation, modulo = min(tau[[min(i, length(tau))]]),
                          mod.only = FALSE, plot = FALSE)$periods
      # guard: if no seasonality detected, allow NULL
      if (!is.numeric(periods) || length(periods) == 0L) periods <- NULL
      
      b <- NNS.ARMA.optim(variable_interpolation, seasonal.factor = periods,
                          obj.fn = obj.fn,
                          objective = objective,
                          print.trace = FALSE,
                          ncores = 1,
                          negative.values = min(variable_interpolation, na.rm = TRUE)<0, h = h)
      
      variable_extrapolation <- b$results
    } else variable_extrapolation <- NULL
    
    return(list(variable_interpolation, variable_extrapolation))
  }
  
  interpolation_results <- lapply(nns_IVs, `[[`, 1)
  
  nns_IVs_interpolated_extrapolated <- data.frame(do.call(cbind, interpolation_results))
  colnames(nns_IVs_interpolated_extrapolated) <- colnames(variables)
  
  positive_values <- apply(variables, 2, function(x) min(x, na.rm = TRUE)>0)
  
  for(i in 1:length(positive_values)){
    if(positive_values[i]) nns_IVs_interpolated_extrapolated[,i] <- pmax(0, nns_IVs_interpolated_extrapolated[,i])
  }
  
  rownames(nns_IVs_interpolated_extrapolated) <- head(dates, nrow(variables))
  colnames(nns_IVs_interpolated_extrapolated) <- colnames(variables)
  
  if(h == 0) {
    options(warn = oldw)
    return(nns_IVs_interpolated_extrapolated)
  }
  
  extrapolation_results <- lapply(nns_IVs, `[[`, 2)
  nns_IVs_results <- data.frame(do.call(cbind, extrapolation_results))
  colnames(nns_IVs_results) <- colnames(variables)
  
  # Combine interpolated / extrapolated / forecasted IVs onto training data.frame
  new_values <- lapply(1:ncol(variables), function(i) c(nns_IVs_interpolated_extrapolated[,i], nns_IVs_results[,i]))
  
  new_values <- data.frame(do.call(cbind, new_values))
  colnames(new_values) <- as.character(colnames(variables))
  
  nns_IVs_interpolated_extrapolated <- head(new_values, nrow(variables))
  
  # Now lag new forecasted data.frame
  lagged_new_values <- lag.mtx(new_values, tau = tau)
  
  # Keep original variables as training set
  lagged_new_values_train <- head(lagged_new_values, nrow(lagged_new_values) - h)
  
  if(status) message("Currently generating multi-variate estimates...", "\r", appendLF = TRUE)
  
  if(num_cores > 1){
    if(status) message("Parallel process running, status unavailable... \n","\r",appendLF=FALSE)
    status <- FALSE
  }
  
  # ==== FIXED MULTIVARIATE BLOCK (shape alignment + clean IVs.test + validation weights) ====
  lists <- foreach(i = 1:ncol(variables),
                   .packages = c("NNS", "data.table"),
                   .export   = c("blend_weight_mse_opt","pool_logistic","gamma_from_validation",
                                 "nns_IVs_interpolated_extrapolated","tau","obj.fn","objective",
                                 "dim.red.method","h","lagged_new_values","nns_IVs_results")) %dopar%{
                                   if(status) message("Variable ", i, " of ", ncol(variables), appendLF = TRUE)
                                   
                                   IV <- lagged_new_values_train[, -i, drop = FALSE]
                                   DV <- lagged_new_values_train[,  i]
                                   
                                   # align by complete cases across DV + IV
                                   train_block <- data.frame(DV = DV, IV)
                                   cc <- stats::complete.cases(train_block)
                                   train_block <- train_block[cc, , drop = FALSE]
                                   
                                   DVc <- as.numeric(train_block[, 1])                    # vector
                                   IVc <- as.matrix(train_block[, -1, drop = FALSE])      # matrix
                                   
                                   # Guard: if no usable rows, fall back to univariate forecast for this target
                                   if (nrow(IVc) < 2) {
                                     return(list(nns_IVs_results[, i, drop = TRUE], colnames(lagged_new_values_train), 0.5, 0))
                                   }
                                   
                                   # ts.test must be an integer; base it on the CLEAN DV length
                                   ts <- as.integer( max(2*h, ceiling(0.2 * length(DVc))) )
                                   
                                   # ---- TEST BLOCK: last h rows of the FULL lagged panel (future slice) ----
                                   IV_test <- as.matrix(tail(lagged_new_values[, -i, drop = FALSE], h))
                                   if (anyNA(IV_test)) {
                                     for (j in seq_len(ncol(IV_test))) {
                                       if (anyNA(IV_test[, j])) {
                                         fill_val <- utils::tail(IVc[, j], 1)
                                         if (length(fill_val) == 0 || is.na(fill_val)) fill_val <- 0
                                         IV_test[is.na(IV_test[, j]), j] <- fill_val
                                       }
                                     }
                                   }
                                   
                                   # ---- Dimension reduction + stacked predictions (future horizon) ----
                                   cor_threshold <- NNS.stack(IVs.train = IVc,
                                                              DV.train  = DVc,
                                                              IVs.test  = IV_test,
                                                              ts.test   = ts,
                                                              folds = 1,
                                                              obj.fn = obj.fn,
                                                              objective = objective,
                                                              method = c(1,2),
                                                              dim.red.method = dim.red.method,
                                                              order = NULL, ncores = 1, stack = TRUE, status = FALSE)
                                   
                                   # Relevance on aligned training block
                                   if(any(dim.red.method == "cor" | dim.red.method == "all")){
                                     rel.1 <- abs(stats::cor(train_block, method = "spearman", use = "pairwise.complete.obs"))
                                   }
                                   if(any(dim.red.method == "nns.dep" | dim.red.method == "all")){
                                     rel.2 <- NNS.dep(train_block)$Dependence
                                   }
                                   if(any(dim.red.method == "nns.caus" | dim.red.method == "all")){
                                     rel.3 <- NNS.caus(train_block)
                                   }
                                   
                                   if(dim.red.method == "cor")      rel_vars_vec <- rel.1[-1,1]
                                   if(dim.red.method == "nns.dep")  rel_vars_vec <- rel.2[-1,1]
                                   if(dim.red.method == "nns.caus") rel_vars_vec <- rel.3[1,-1]
                                   if(dim.red.method == "all")      rel_vars_vec <- ((rel.1+rel.2+rel.3)/3)[1, -1]
                                   
                                   rel_vars <- names(rel_vars_vec[rel_vars_vec > cor_threshold$NNS.dim.red.threshold])
                                   rel_vars <- stats::na.omit(rel_vars)
                                   if(any(length(rel_vars)==0 | is.null(rel_vars))){
                                     rel_vars <- colnames(lagged_new_values_train)
                                   }
                                   
                                   nns_DVs <- cor_threshold$stack
                                   nns_DVs[is.na(nns_DVs)] <- nns_IVs_results[is.na(nns_DVs),i]
                                   
                                   # --- Validation predictions for weight learning ---
                                   # Multivariate validation: stack on last ts in-sample rows
                                   m_val <- NNS.stack(
                                     IVs.train = IVc, DV.train = DVc,
                                     IVs.test  = tail(IVc, ts),
                                     ts.test   = ts, folds = 1,
                                     obj.fn = obj.fn, objective = objective,
                                     method = c(1, 2), dim.red.method = dim.red.method,
                                     order = NULL, ncores = 1, stack = TRUE, status = FALSE
                                   )$stack
                                   
                                   # Univariate validation from imputed target series:
                                   x_i <- nns_IVs_interpolated_extrapolated[, i]
                                   periods_val <- NNS.seas(x_i, modulo = min(tau[[min(i, length(tau))]]),
                                                           mod.only = FALSE, plot = FALSE)$periods
                                   if (!is.numeric(periods_val) || length(periods_val) == 0L) periods_val <- NULL
                                   
                                   u_val <- NNS.ARMA.optim(
                                     variable = head(x_i, length(x_i) - ts),
                                     seasonal.factor = periods_val,
                                     obj.fn = obj.fn, objective = objective,
                                     print.trace = FALSE, ncores = 1,
                                     negative.values = (min(x_i, na.rm = TRUE) < 0),
                                     h = ts
                                   )$results
                                   
                                   y_val <- tail(DVc, ts)
                                   
                                   # dynamic, no-constant gamma from validation
                                   gamma_i_val <- gamma_from_validation(y_val, u_val, m_val, method = "wilcox")
                                   
                                   # validation-based weight (still clipped/shrunk internally if requested; set shrink=NULL for pure)
                                   w_val_i <- blend_weight_mse_opt(y_val, u_val, m_val, shrink = NULL)
                                   
                                   list(nns_DVs, rel_vars, w_val_i, gamma_i_val)
                                 }
  
  if(num_cores > 1) {
    doParallel::stopImplicitCluster()
    foreach::registerDoSEQ()
    invisible(data.table::setDTthreads(0, throttle = NULL))
    invisible(gc(verbose = FALSE))
  }
  
  nns_DVs <- lapply(lists, `[[`, 1)
  relevant_vars <- lapply(lists, `[[`, 2)
  w_val_list   <- sapply(lists, `[[`, 3)
  gamma_v_list <- sapply(lists, `[[`, 4)
  
  nns_DVs <- data.frame(do.call(cbind, nns_DVs))
  nns_DVs <- head(nns_DVs, h)
  
  RV <- lapply(relevant_vars, function(x) if(length(x)==0){NA} else {x})
  
  colnames(nns_DVs) <- colnames(variables)
  
  RV <- do.call(cbind, lapply(RV, `length<-`, max(lengths(RV))))
  colnames(RV) <- as.character(colnames(variables))
  
  multi <- uni <- numeric(length(colnames(RV)))
  
  for(i in 1:length(colnames(RV))){
    if(length(stats::na.omit(RV[,i])) > 0){
      given_var <- unlist(strsplit(colnames(RV)[i], split = "_tau"))[1]
      observed_var <- do.call(rbind,(strsplit(stats::na.omit(RV[,i]), split = "_tau")))[,1]
      
      equal_tau <- sum(given_var==observed_var)
      unequal_tau <- sum(given_var!=observed_var)
      
      # tau-based heuristic
      w_tau <- if(naive.weights) 0.5 else equal_tau/(equal_tau + unequal_tau)
      
      # validation-based weight and fully dynamic gamma
      w_val_i <- w_val_list[i]
      gamma_i <- gamma_v_list[i]
      
      # pool in log-odds; no extra shrink to avoid new constants
      w_star  <- pool_logistic(w_tau, w_val_i, gamma = gamma_i, shrink = 0)
      
      uni[i]  <- w_star
      multi[i] <- 1 - w_star
    } else {
      uni[i] <- 0.5
      multi[i] <- 0.5
    }
  }
  
  forecasts <- data.frame(Reduce(`+`,list(t(t(nns_IVs_results)*uni) , t(t(nns_DVs)*multi))))
  colnames(forecasts) <- colnames(variables)
  
  colnames(nns_IVs_results) <- colnames(variables)
  rownames(nns_IVs_results) <- tail(dates, h)
  colnames(nns_DVs) <- colnames(variables)
  rownames(nns_DVs) <- tail(dates, h)
  colnames(forecasts) <- colnames(variables)
  rownames(forecasts) <- tail(dates, h)
  rownames(nns_IVs_interpolated_extrapolated) <- head(dates, nrow(nns_IVs_interpolated_extrapolated))
  
  options(warn = oldw)
  
  return( list("interpolated_and_extrapolated" = nns_IVs_interpolated_extrapolated,
               "relevant_variables" = data.frame(RV),
               univariate = nns_IVs_results,
               multivariate = nns_DVs,
               ensemble = forecasts) )
  
}
