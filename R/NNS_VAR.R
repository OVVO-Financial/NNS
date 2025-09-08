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
  
  # --- helpers for data-driven blending via partial moments ---
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
  
  # --- lag.mtx ---
  lag.mtx <- function(x, tau) {
    max_tau <- max(unlist(tau))
    
    # --- Normalize x to a matrix and get clean base names ---
    if (is.null(dim(x))) {
      # x is a vector; try to recover the original colname from the call
      mc <- match.call()
      base_name <- NULL
      if (is.call(mc$x) && identical(mc$x[[1L]], as.name("["))) {
        # mc$x is like: new_values[, i, drop = TRUE]
        pf <- parent.frame()
        base_obj <- try(eval(mc$x[[2L]], envir = pf), silent = TRUE)  # new_values
        col_idx  <- try(eval(mc$x[[3L]], envir = pf), silent = TRUE)  # i  (may be symbol)
        if (!inherits(base_obj, "try-error") && !is.null(colnames(base_obj))) {
          col_idx <- try(as.integer(col_idx), silent = TRUE)
          if (!inherits(col_idx, "try-error") && length(col_idx) == 1L &&
              col_idx >= 1L && col_idx <= ncol(base_obj)) {
            base_name <- colnames(base_obj)[col_idx]
          }
        }
      }
      x <- matrix(x, ncol = 1L)
      colnames(x) <- if (!is.null(base_name)) base_name else "V1"
    } else {
      x <- as.matrix(x)
      if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
    }
    
    p <- ncol(x)
    j.vectors <- vector("list", p)
    
    for (j in seq_len(p)) {
      colhead <- colnames(x)[j]
      heads <- paste0(colhead, "_tau_")              # keep your internal style: *_tau_*
      heads <- gsub('"', '', heads, fixed = TRUE)
      
      x.vectors <- vector("list", max_tau + 1L)
      names(x.vectors) <- paste0(heads, 0:max_tau)   # name slots first
      
      for (i in 0:max_tau) {
        start <- max_tau - i + 1L
        end   <- nrow(x) - i
        x.vectors[[i + 1L]] <- x[start:end, j]
      }
      j.vectors[[j]] <- do.call(cbind, x.vectors)
    }
    
    mtx <- as.data.frame(do.call(cbind, j.vectors), check.names = FALSE)
    
    # If tau is a list of per-series lags, select those; else keep all 0..max_tau
    if (length(unlist(tau)) > 1L) {
      # Blocks are (max_tau+1) wide per original series
      block <- max_tau + 1L
      relevant <- unlist(lapply(seq_along(tau), function(i) {
        off <- (i - 1L) * block
        c(off + 1L, off + unlist(tau[[i]]) + 1L)   # include tau_0 plus requested lags
      }))
      mtx <- mtx[, sort(unique(relevant)), drop = FALSE]
    }
    
    # Move all tau_0 columns to the front (your original behavior)
    vars0 <- grep("tau_0$", colnames(mtx))
    rest  <- setdiff(seq_len(ncol(mtx)), vars0)
    mtx   <- mtx[, c(vars0, rest), drop = FALSE]
    
    mtx
  }
  
  # LPM/UPM-based validation weight and credibility (no constants)
  # delta = (y - u)^2 - (y - m)^2 ; negative => univariate better, positive => multivariate better
  compute_wval_gamma_lpm <- function(y_val, u_val, m_val) {
    delta <- (y_val - u_val)^2 - (y_val - m_val)^2
    L <- NNS::LPM.ratio(degree = 1, target = 0, variable = delta)  # standardized mass below 0
    U <- NNS::UPM.ratio(degree = 1, target = 0, variable = delta)  # standardized mass above 0
    if (!is.finite(L)) L <- 0
    if (!is.finite(U)) U <- 0
    S <- L + U
    if (!is.finite(S) || S == 0) return(c(w_val = 0.5, gamma = 0))
    w_val <- L / S
    # gamma reflects credibility of validation, mapped to [0,1] by the same mass S
    gamma <- S / (S + 1)
    c(w_val = w_val, gamma = gamma)
  }
  
  # ---------- Input handling ----------
  if(any(class(variables) == "ts")){
    dates <- zoo::as.yearmon(zoo::index(variables))
    variables <- data.frame(zoo::coredata(variables))
  }
  
  if(is.null(colnames(variables))) colnames(variables) <- paste0("V", 1:ncol(variables))
  dim.red.method <- tolower(dim.red.method)
  
  # cores / parallel setup
  num_cores <- if(is.null(ncores)) max(1, parallel::detectCores() - 1) else ncores
  if(num_cores > 1){
    data.table::setDTthreads(1L)
    doParallel::registerDoParallel(cores = num_cores)
  } else {
    foreach::registerDoSEQ()
  }
  
  # ---------- Interpolate NA and extrapolate to align mixed frequencies ----------
  if(status) message("Interpolating / extrapolating variables...", appendLF = TRUE)
  
  nns_IVs_interpolated_extrapolated <- variables
  
  for(i in 1:ncol(variables)){
    xi <- variables[, i]
    
    # fast linear interpolation for NAs
    if(anyNA(xi)){
      not_na <- !is.na(xi)
      xi[!not_na] <- approx(seq_along(xi)[not_na], xi[not_na], xout = which(!not_na), method = "linear", rule = 2)$y
    }
    
    # extrapolate tail using NNS.ARMA.optim
    seasonal.periods <- tryCatch({
      NNS::NNS.seas(xi, modulo = min(if(is.list(tau)) tau[[min(i, length(tau))]] else tau),
                    mod.only = FALSE, plot = FALSE)$periods
    }, error = function(e) NULL)
    if(!is.numeric(seasonal.periods) || length(seasonal.periods) == 0L) seasonal.periods <- NULL
    
    extrap <- NNS::NNS.ARMA.optim(variable = xi,
                                  seasonal.factor = seasonal.periods,
                                  obj.fn = obj.fn,
                                  objective = objective,
                                  print.trace = FALSE,
                                  ncores = 1,
                                  negative.values = (min(xi, na.rm = TRUE) < 0),
                                  h = 0)$results
    nns_IVs_interpolated_extrapolated[, i] <- extrap
  }
  
  # Early exit if only imputation/extrapolation requested
  if(h == 0){
    if(!is.null(dates)) rownames(nns_IVs_interpolated_extrapolated) <- dates
    options(warn = oldw)
    return(list("interpolated_and_extrapolated" = nns_IVs_interpolated_extrapolated))
  }
  
  # ---------- Univariate forecasts (per series) ----------
  if(status) message("Computing univariate forecasts...", appendLF = TRUE)
  
  univariate_list <- foreach::foreach(i = 1:ncol(nns_IVs_interpolated_extrapolated),
                                      .packages = c("NNS"),
                                      .export   = c()) %dopar% {
                                        
                                        xi <- nns_IVs_interpolated_extrapolated[, i]
                                        seasonal.periods <- tryCatch({
                                          NNS::NNS.seas(xi, modulo = min(if(is.list(tau)) tau[[min(i, length(tau))]] else tau),
                                                        mod.only = FALSE, plot = FALSE)$periods
                                        }, error = function(e) NULL)
                                        if(!is.numeric(seasonal.periods) || length(seasonal.periods) == 0L) seasonal.periods <- NULL
                                        
                                        fit <- NNS::NNS.ARMA.optim(variable = xi,
                                                                   seasonal.factor = seasonal.periods,
                                                                   obj.fn = obj.fn,
                                                                   objective = objective,
                                                                   print.trace = FALSE,
                                                                   ncores = 1,
                                                                   negative.values = (min(xi, na.rm = TRUE) < 0),
                                                                   h = h)$results
                                        fit
                                      }
  
  nns_IVs_results <- data.frame(do.call(cbind, univariate_list))
  colnames(nns_IVs_results) <- colnames(variables)
  
  # ---------- Build lagged panel ----------
  new_values <- nns_IVs_interpolated_extrapolated
  lagged_new_values <- do.call(cbind, lapply(seq_len(ncol(new_values)), function(i){
    k <- if (is.list(tau)) tau[[min(i, length(tau))]] else tau
    lag.mtx(new_values[, i, drop = FALSE], tau = k)  
  }))
  
 
  # training rows (drop NAs arising from lagging)
  row_keep <- stats::complete.cases(lagged_new_values)
  lagged_new_values_train <- lagged_new_values[row_keep, , drop = FALSE]
  colnames(lagged_new_values_train) <- colnames(lagged_new_values)
  
  
  # ---------- Multivariate block (alignment + NA guard + single NNS.stack call) ----------
  if(status) message("Computing multivariate stacks & relevance...", appendLF = TRUE)
  
  lists <- foreach::foreach(i = 1:ncol(variables),
                            .packages = c("NNS", "data.table"),
                            .export   = c("pool_logistic","compute_wval_gamma_lpm",
                                          "tau","obj.fn","objective",
                                          "dim.red.method","h","lagged_new_values","nns_IVs_interpolated_extrapolated","nns_IVs_results")) %dopar%{
                                            
                                            if(status) message("Variable ", i, " of ", ncol(variables), appendLF = TRUE)
                                            
                                            IV <- lagged_new_values_train[, -i, drop = FALSE]
                                            DV <- lagged_new_values_train[,  i]
                                            
                                            iv_names <- colnames(lagged_new_values_train)[-i]
                                            colnames(IV) <- iv_names
                                            
                                            # align by complete cases across DV + IV
                                            train_block <- data.frame(DV = DV, IV)
                                            cc <- stats::complete.cases(train_block)
                                            train_block <- train_block[cc, , drop = FALSE]
                                            
                                            colnames(train_block) <- c("DV", iv_names)
                                            
                                            DVc <- as.numeric(train_block[, 1])               # vector
                                            IVc <- as.matrix(train_block[, -1, drop = FALSE]) # matrix
                                            
                                            # Guard: if no usable rows, fall back to univariate forecast for this target
                                            if (nrow(IVc) < 2) {
                                              return(list(nns_IVs_results[, i, drop = TRUE],
                                                          colnames(lagged_new_values_train), 0.5, 0))
                                            }
                                            
                                            # ts.test must be an integer; base it on the CLEAN DV length
                                            ts <- as.integer( max(2*h, ceiling(0.2 * length(DVc))) )
                                            
                                            # ---- TEST BLOCKS ----
                                            # FUTURE: last h rows from full panel (may contain NA; fill from last observed in train)
                                            IV_test <- as.matrix(utils::tail(lagged_new_values[, -i, drop = FALSE], h))
                                            if (anyNA(IV_test)) {
                                              for (j in seq_len(ncol(IV_test))) {
                                                if (anyNA(IV_test[, j])) {
                                                  fill_val <- utils::tail(IVc[, j], 1)
                                                  if (length(fill_val) == 0 || is.na(fill_val)) fill_val <- 0
                                                  IV_test[is.na(IV_test[, j]), j] <- fill_val
                                                }
                                              }
                                            }
                                            
                                            # VALIDATION: last ts rows of in-sample IVs
                                            IV_val <- as.matrix(utils::tail(IVc, ts))
                                            
                                            # ---- ONE call to NNS.stack: rbind FUTURE and VALIDATION, then split ----
                                            ct <- NNS::NNS.stack(
                                              IVs.train = IVc, DV.train = DVc,
                                              IVs.test  = rbind(IV_test, IV_val),
                                              ts.test   = ts, folds = 1,
                                              obj.fn = obj.fn, objective = objective,
                                              method = c(1, 2), dim.red.method = dim.red.method,
                                              order = NULL, ncores = 1, stack = TRUE, status = FALSE
                                            )
                                            
                                            stack_all <- drop(ct$stack)                 # vector
                                            
                                            n_future <- nrow(IV_test)
                                            n_val    <- ts
                                            n_total  <- n_future + n_val
                                            n_preds  <- length(stack_all)
                                            
                                            # (Optional) defensive guard
                                            if (n_preds != n_total) {
                                              n_future <- min(n_preds, n_future)
                                              n_val    <- min(n_preds - n_future, n_val)
                                            }
                                            
                                            # Split and **drop names** so they don't carry rownames through
                                            nns_DVs <- unname(utils::head(stack_all, n_future))  # FUTURE h
                                            m_val   <- unname(utils::tail(stack_all, n_val))     # VALIDATION ts
                                           
                                            # Relevance on aligned training block (use the same threshold object)
                                            if(any(dim.red.method == "cor" | dim.red.method == "all")){
                                              rel.1 <- abs(stats::cor(train_block, method = "spearman", use = "pairwise.complete.obs"))
                                            }
                                            if(any(dim.red.method == "nns.dep" | dim.red.method == "all")){
                                              rel.2 <- NNS::NNS.dep(train_block)$Dependence
                                            }
                                            if(any(dim.red.method == "nns.caus" | dim.red.method == "all")){
                                              rel.3 <- NNS::NNS.caus(train_block)
                                            }
                                            
                                            if(dim.red.method == "cor")      rel_vars <- rel.1[-1,1]
                                            if(dim.red.method == "nns.dep")  rel_vars <- rel.2[-1,1]
                                            if(dim.red.method == "nns.caus") rel_vars <- rel.3[1,-1]
                                            if(dim.red.method == "all")      rel_vars <- ((rel.1+rel.2+rel.3)/3)[1, -1]
                                            
                                            
                                            rel_vars <- names(rel_vars[rel_vars > ct$NNS.dim.red.threshold])
                                            rel_vars <- rel_vars[rel_vars!=i]
                                            rel_vars <- na.omit(rel_vars)
                                            
                                            if(any(length(rel_vars)==0 | is.null(rel_vars))){
                                              rel_vars <- colnames(lagged_new_values_train)
                                            }
                                            
                                            # Fill any NA in FUTURE stack from univariate fallback
                                            nns_DVs[is.na(nns_DVs)] <- nns_IVs_results[is.na(nns_DVs), i]
                                            
                                            # --- Validation predictions for weight learning ---
                                            # Univariate validation from imputed target series:
                                            x_i <- nns_IVs_interpolated_extrapolated[, i]
                                            periods_val <- NNS::NNS.seas(x_i, modulo = min(if(is.list(tau)) tau[[min(i, length(tau))]] else tau),
                                                                         mod.only = FALSE, plot = FALSE)$periods
                                            if (!is.numeric(periods_val) || length(periods_val) == 0L) periods_val <- NULL
                                            
                                            u_val <- NNS::NNS.ARMA.optim(
                                              variable = utils::head(x_i, length(x_i) - ts),
                                              seasonal.factor = periods_val,
                                              obj.fn = obj.fn, objective = objective,
                                              print.trace = FALSE, ncores = 1,
                                              negative.values = (min(x_i, na.rm = TRUE) < 0),
                                              h = ts
                                            )$results
                                            
                                            y_val <- utils::tail(DVc, ts)
                                            
                                            # LPM/UPM-based validation weight and credibility (no constants)
                                            wg <- compute_wval_gamma_lpm(y_val, u_val, m_val)
                                            w_val_i   <- wg[["w_val"]]
                                            gamma_val <- wg[["gamma"]]
                                            
                                            list(nns_DVs, rel_vars, w_val_i, gamma_val)
                                          }
  
  if(num_cores > 1) {
    doParallel::stopImplicitCluster()
    foreach::registerDoSEQ()
    invisible(data.table::setDTthreads(0, throttle = NULL))
    invisible(gc(verbose = FALSE))
  }
  
  nns_DVs       <- lapply(lists, `[[`, 1)
  relevant_vars <- lapply(lists, `[[`, 2)
  w_val_list    <- as.numeric(sapply(lists, `[[`, 3))
  gamma_v_list  <- as.numeric(sapply(lists, `[[`, 4))

  nns_DVs <- data.frame(do.call(cbind, nns_DVs))
  nns_DVs <- utils::head(nns_DVs, h)
  
  colnames(nns_DVs) <- colnames(variables)
  
  
  RV <- lapply(relevant_vars, function(x) if(length(x)==0){NA} else {x})
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
      
      # validation-based weight and fully dynamic gamma (from LPM/UPM)
      w_val_i <- w_val_list[i]
      gamma_i <- gamma_v_list[i]
      
      # pool in log-odds; no extra shrink to avoid constants
      w_star  <- pool_logistic(w_tau, w_val_i, gamma = gamma_i, shrink = 0)
      
      uni[i]   <- w_star
      multi[i] <- 1 - w_star
    } else {
      uni[i] <- 0.5
      multi[i] <- 0.5
    }
  }
  
  forecasts <- data.frame(Reduce(`+`, list(t(t(nns_IVs_results)*uni),
                                           t(t(nns_DVs)*multi))))
  colnames(forecasts) <- colnames(variables)
  
  if (!is.null(dates)) {
    rownames(nns_IVs_results) <- tail(dates, h)
    rownames(nns_DVs)         <- tail(dates, h)
    rownames(forecasts)       <- tail(dates, h)
  } else {
    rn <- seq_len(h)
    rownames(nns_IVs_results) <- rn
    rownames(nns_DVs)         <- rn
    rownames(forecasts)       <- rn
  }
  
  options(warn = oldw)
  
  return( list("interpolated_and_extrapolated" = nns_IVs_interpolated_extrapolated,
               "relevant_variables" = data.frame(RV),
               univariate  = nns_IVs_results,
               multivariate = nns_DVs,
               ensemble = forecasts) )
}