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
#' @param nowcast logical; \code{FALSE} (default) internal call for frequency alignment in downstream nowcasting applications.
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
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
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
                    obj.fn = expression( {
                      .denom <- NNS::Co.LPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual)) +
                                NNS::Co.UPM(1, predicted, actual, target_x = mean(predicted), target_y = mean(actual))
                      if (.denom == 0) Inf else mean((predicted - actual)^2) / .denom
                    } ),
                    objective = "min",
                    status = TRUE,
                    ncores = NULL,
                    nowcast = FALSE){
  
  # Capture whether the built-in ratio objective is in force before anything
  # touches obj.fn. Only then may NNS.VAR fall back to plain MSE when the
  # co-movement denominator is degenerate for every candidate.
  use_default_objective <- missing(obj.fn) && identical(objective, "min")

  oldw <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = oldw), add = TRUE)

  dates <- NULL
  
  # ===================== Lag builder (robust names) =====================
  lag.mtx <- function(x, tau) {
    max_tau <- max(unlist(tau))
    if (is.null(dim(x))) {
      mc <- match.call(); base_name <- NULL
      if (is.call(mc$x) && identical(mc$x[[1L]], as.name("["))) {
        pf <- parent.frame()
        base_obj <- try(eval(mc$x[[2L]], envir = pf), silent = TRUE)
        col_idx  <- try(eval(mc$x[[3L]], envir = pf), silent = TRUE)
        if (!inherits(base_obj, "try-error") && !is.null(colnames(base_obj))) {
          col_idx <- try(as.integer(col_idx), silent = TRUE)
          if (!inherits(col_idx, "try-error") && length(col_idx) == 1L &&
              col_idx >= 1L && col_idx <= ncol(base_obj)) {
            base_name <- colnames(base_obj)[col_idx]
          }
        }
      }
      x <- matrix(x, ncol = 1L); colnames(x) <- if (!is.null(base_name)) base_name else "V1"
    } else {
      x <- as.matrix(x); if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
    }
    p <- ncol(x); j.vectors <- vector("list", p)
    for (j in seq_len(p)) {
      colhead <- colnames(x)[j]
      heads <- gsub('"','', paste0(colhead, "_tau_"), fixed = TRUE)
      x.vectors <- vector("list", max_tau + 1L); names(x.vectors) <- paste0(heads, 0:max_tau)
      for (i in 0:max_tau) {
        start <- max_tau - i + 1L; end <- nrow(x) - i
        x.vectors[[i + 1L]] <- x[start:end, j]
      }
      j.vectors[[j]] <- do.call(cbind, x.vectors)
    }
    mtx <- as.data.frame(do.call(cbind, j.vectors), check.names = FALSE)
    if (length(unlist(tau)) > 1L) {
      block <- max_tau + 1L
      relevant <- unlist(lapply(seq_along(tau), function(i) {
        off <- (i - 1L) * block
        c(off + 1L, off + unlist(tau[[i]]) + 1L)
      }))
      mtx <- mtx[, sort(unique(relevant)), drop = FALSE]
    }
    vars0 <- grep("tau_0$", colnames(mtx)); rest <- setdiff(seq_len(ncol(mtx)), vars0)
    mtx[, c(vars0, rest), drop = FALSE]
  }
  
  # --------- Time-index dates (base R, no external date packages) ---------
  # `dates` are used only as row labels on the returned matrices, so they are
  # produced as "Mon YYYY" character labels.  Numerical forecasts never depend
  # on them.  A yearmon-equivalent numeric (year + (month-1)/12) is used for the
  # month arithmetic, then formatted -- matching zoo::as.yearmon's appearance.
  .index_to_yearmon <- function(idx) {
    if (inherits(idx, c("Date", "POSIXct", "POSIXt"))) {
      lt <- as.POSIXlt(idx)
      (lt$year + 1900) + lt$mon / 12
    } else {
      suppressWarnings(as.numeric(idx))
    }
  }
  .format_yearmon <- function(ym) {
    yr <- floor(ym + 1e-8)
    mo <- as.integer(round((ym - yr) * 12)) + 1L
    roll <- !is.na(mo) & mo > 12L
    yr[roll] <- yr[roll] + 1L; mo[roll] <- 1L
    out <- rep(NA_character_, length(ym))
    ok  <- !is.na(mo) & mo >= 1L & mo <= 12L
    out[ok] <- paste(month.abb[mo[ok]], yr[ok])
    out
  }

  # Extract a time index (best effort) without zoo/xts.
  idx <- NULL
  if (inherits(variables, "ts")) {
    idx <- as.numeric(stats::time(variables))
  } else if (inherits(variables, c("xts", "zoo"))) {
    raw <- attr(variables, "index")
    if (is.null(raw)) raw <- suppressWarnings(try(stats::time(variables), silent = TRUE))
    if (!is.null(raw) && !inherits(raw, "try-error")) {
      # xts/zoo store the index as numeric seconds since the epoch.
      if (is.numeric(raw)) idx <- as.POSIXct(raw, origin = "1970-01-01", tz = "UTC") else idx <- raw
    }
  }

  if (nowcast && !is.null(idx)) {
    ym <- .index_to_yearmon(idx)
    if (!any(is.na(ym))) dates <- .format_yearmon(c(ym, tail(ym, h) + h/12))
  }

  if(any(class(variables)%in%c("tbl","data.table"))) variables <- as.data.frame(variables)
  if (inherits(variables, "xts")) {
    if (is.null(dates) && !is.null(idx)) dates <- .format_yearmon(.index_to_yearmon(idx))
    variables <- data.frame(matrix(as.numeric(variables), nrow = NROW(variables),
                                   dimnames = list(NULL, colnames(variables))),
                            check.names = FALSE)
  }
  if (inherits(variables, "ts")) {
    if (is.null(dates) && !is.null(idx)) dates <- .format_yearmon(.index_to_yearmon(idx))
    variables <- data.frame(as.matrix(variables), check.names = FALSE)
  }
  
  dim.red.method <- tolower(dim.red.method)
  # dim.red.method is scalar: the downstream selectors use scalar if() (which
  # errors on a length > 1 condition in R >= 4.2). Use "all" to combine methods.
  if(length(dim.red.method) != 1L || !(dim.red.method %in% c("cor","nns.dep","nns.caus","all"))){ stop('Please ensure the dimension reduction method is set to a single value: one of "cor", "nns.dep", "nns.caus" or "all".')}
  
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
  auto_cores <- is.null(ncores)
  if (auto_cores) {
    num_cores <- max(1L, as.integer(parallel::detectCores()) - 1L, na.rm = TRUE)
  } else {
    num_cores <- max(1L, as.integer(ncores), na.rm = TRUE)
  }

  # Workload gate: standing up a process cluster costs ~1s (spawn + exporting
  # data to workers). NNS.VAR parallelizes across VARIABLES -- one heavy task per
  # column (NNS.stack + dimension reduction over the ~ncol*tau lagged regressors)
  # -- so the parallel work grows with BOTH the number of variables (task count)
  # and the series length. A wide VAR (many variables) therefore has ample
  # parallel work even on a short series: e.g. a 35-variable monthly nowcast with
  # ~300 rows is 35 heavy tasks that parallelize well, yet nrow < 500. Gate the
  # AUTOMATIC default (ncores = NULL) on either dimension -- enough variables OR a
  # long series -- and always honor an explicit ncores. (Thresholds are heuristic;
  # tune ncol_parallel_min / the nrow cutoff against representative workloads.)
  ncol_parallel_min <- 6L
  use_parallel <- num_cores > 1 &&
    ncol(variables) >= 2L &&
    (!auto_cores ||
       ncol(variables) >= ncol_parallel_min ||
       nrow(variables) >= 500L)

  # Create a single cluster shared by both parallel sections (fork-first, with a
  # PSOCK fallback). Stopped at the end of the multi-variate estimate section.
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
    on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  }
  
  if(status) message("Currently interpolating/extrapolating variables...","\r", appendLF=TRUE)
  
  # ===================== Interpolation / Extrapolation  =====================
  .interp_worker <- function(i) {
    n <- nrow(variables)
    index <- seq_len(n)
    a <- cbind.data.frame("index" = index, variables)

    # For Interpolation / Extrapolation of all missing values
    selected_variable <- a[, c(1,(i+1))]

    interpolation_point <- tail(which(!is.na(selected_variable[,2])), 1)
    
    missing_index <- which(is.na(selected_variable[,2]))
    selected_variable <- selected_variable[complete.cases(selected_variable), , drop = FALSE]
    
    h_int <- tail(index, 1) - interpolation_point
    # ensure plain numeric to avoid classed assignment issues
    variable_interpolation <- as.numeric(variables[,i])
    
    if (length(missing_index) == 0L) {
      # ---- FIX: dataset is complete -> DO NOT SMOOTH ----
      # keep the original series exactly
      variable_interpolation <- as.numeric(variables[, i])
      
    } else if (h_int > 0) {
      # trailing NA(s): estimate them using NNS.stack on the index (as in original).
      # The index is duplicated into two identical predictors so the multivariate
      # machinery has >1 column; strip cbind's auto-generated column names so the
      # duplicate labels do not trip NNS.stack's unique-predictor-name guard.
      multi <- NNS.stack(unname(cbind(selected_variable[,1], selected_variable[,1])), selected_variable[,2],
                         order = NULL, ncores = 1, status = FALSE, folds = 5,
                         IVs.test = unname(cbind(missing_index, missing_index)), method = 1)$stack
      variable_interpolation[missing_index] <- as.numeric(multi)
      
    } else {
      # interior NA(s) only: fit on index, but fill ONLY the missing indices (no global smoothing)
      fitted_missing <- NNS.reg(selected_variable[,1], selected_variable[,2],
                                order = "max", ncores = 1,
                                point.est = missing_index, plot = FALSE, point.only = TRUE)$Point.est
      if (length(missing_index)) variable_interpolation[missing_index] <- as.numeric(fitted_missing)
    }
    
    if(h > 0){
      # robust tau selection without changing flow
      tau_i <- if (is.list(tau)) tau[[min(i, length(tau))]] else tau
      periods <- tryCatch(NNS.seas(variable_interpolation, modulo = min(tau_i),
                                   mod.only = FALSE, plot = FALSE)$periods,
                          error = function(e) NULL)
      if (!is.numeric(periods) || length(periods) == 0L) periods <- NULL
      
      b <- NNS.ARMA.optim(variable_interpolation, seasonal.factor = periods,
                          obj.fn = obj.fn,
                          objective = objective,
                          print.trace = FALSE,
                          ncores = 1,
                          negative.values = min(variable_interpolation, na.rm = TRUE) < 0, h = h)
      
      variable_extrapolation <- b$results
      
    } else variable_extrapolation <- NULL

    list(variable_interpolation, variable_extrapolation)
  }

  nns_IVs <- if (use_parallel) {
    parallel::clusterExport(
      cl,
      varlist = c("variables", "h", "tau", "obj.fn", "objective"),
      envir = environment()
    )
    parallel::parLapply(cl, 1:ncol(variables), .interp_worker)
  } else {
    lapply(1:ncol(variables), .interp_worker)
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
  
  if(h == 0) return(.NNS.df(nns_IVs_interpolated_extrapolated))
  
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
  
  
  if(use_parallel){
    if(status) message("Parallel process running, status unavailable... \n","\r",appendLF=FALSE)
    status <- FALSE
  }
  
  
  .model_worker <- function(i) {
    if(status) message("Variable ", i, " of ", ncol(variables), appendLF = TRUE)
    
    IV <- lagged_new_values_train[, -i]
    DV <- lagged_new_values_train[, i]
    
    ts <- 2*h
    ts <- max(ts, .2*length(DV))
    
    # Dimension reduction NNS.reg to reduce variables
    .run_var_stack <- function(stack_obj_fn, stack_objective) {
      NNS.stack(IVs.train = IV,
                DV.train = DV,
                IVs.test = tail(IV, h),
                ts.test = ts,
                folds = 1,
                obj.fn = stack_obj_fn,
                objective = stack_objective,
                method = c(1, 2),
                dim.red.method = dim.red.method,
                order = NULL, ncores = 1, stack = TRUE, status = FALSE)
    }

    cor_threshold <- tryCatch(
      .run_var_stack(obj.fn, objective),
      error = function(e) {
        no_finite_method1 <- identical(
          conditionMessage(e),
          paste0("No Method 1 candidate produced a finite ",
                 "complete-coverage OOF objective.")
        )
        if (!use_default_objective || !no_finite_method1) stop(e)
        if (status) {
          message(paste0("NNS.VAR: the co-movement objective was undefined for ",
                         "every candidate; retrying this variable using MSE."))
        }
        .run_var_stack(expression(mean((predicted - actual)^2)), "min")
      }
    )
    
    
    
    if(dim.red.method %in% c("cor", "all")){
      rel.1 <- abs(cor(cbind(DV, IV), method = "spearman"))
    }
    
    if(dim.red.method %in% c("nns.dep", "all")){
      rel.2 <- NNS.dep(cbind(DV, IV))$Dependence
    }
    
    if(dim.red.method %in% c("nns.caus", "all")){
      rel.3 <- NNS.caus(cbind(DV, IV))
    }
    
    if(dim.red.method == "cor") rel_vars <- rel.1[-1,1]
    
    if(dim.red.method == "nns.dep") rel_vars <- rel.2[-1,1]
    
    if(dim.red.method == "nns.caus") rel_vars <- rel.3[1,-1]
    
    if(dim.red.method == "all") rel_vars <- ((rel.1+rel.2+rel.3)/3)[1, -1]
    
    rel_vars <- names(rel_vars[rel_vars > cor_threshold$NNS.dim.red.threshold])
    # (The self-variable's tau_0 column is already excluded from IV via [, -i]
    # above; the former rel_vars[rel_vars != i] compared character names to the
    # integer index i and matched nothing -- a no-op -- so it is removed.)
    rel_vars <- na.omit(rel_vars)
    
    if(length(rel_vars) == 0){
      rel_vars <- colnames(lagged_new_values_train)
    }
    
    nns_DVs <- cor_threshold$stack
    nns_DVs[is.na(nns_DVs)] <- nns_IVs_results[is.na(nns_DVs),i]

    list(nns_DVs, rel_vars)
  }

  lists <- if (use_parallel) {
    parallel::clusterExport(
      cl,
      varlist = c("status", "variables", "lagged_new_values_train", "h",
                  "obj.fn", "objective", "use_default_objective",
                  "dim.red.method", "nns_IVs_results"),
      envir = environment()
    )
    parallel::parLapply(cl, 1:ncol(variables), .model_worker)
  } else {
    lapply(1:ncol(variables), .model_worker)
  }

  if(use_parallel) {
    parallel::stopCluster(cl)
    cl <- NULL
    invisible(gc(verbose = FALSE))
  }

  nns_DVs <- lapply(lists, `[[`, 1)
  relevant_vars <- lapply(lists, `[[`, 2)
  
  
  nns_DVs <- data.frame(do.call(cbind, nns_DVs))
  nns_DVs <- head(nns_DVs, h)
  
  RV <- lapply(relevant_vars, function(x) if(length(x)==0){NA} else {x})
  
  colnames(nns_DVs) <- colnames(variables)
  
  RV <- do.call(cbind, lapply(RV, `length<-`, max(lengths(RV))))
  colnames(RV) <- as.character(colnames(variables))
  
  multi <- uni <- numeric(length(colnames(RV)))
  
  for(i in 1:length(colnames(RV))){
    if(length(na.omit(RV[,i])) > 0){
      given_var <- unlist(strsplit(colnames(RV)[i], split = "_tau"))[1]
      observed_var <- do.call(rbind,(strsplit(na.omit(RV[,i]), split = "_tau")))[,1]
      
      equal_tau <- sum(given_var==observed_var)
      unequal_tau <- sum(given_var!=observed_var)
      
      if(naive.weights) uni[i] <- 0.5 else uni[i] <- equal_tau/(equal_tau + unequal_tau)
      multi[i] <- 1 - uni[i]
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
  
  
  return( list("interpolated_and_extrapolated" = .NNS.df(nns_IVs_interpolated_extrapolated),
               "relevant_variables" = .NNS.df(data.frame(RV)),
               univariate = .NNS.df(nns_IVs_results),
               multivariate = .NNS.df(nns_DVs),
               ensemble = .NNS.df(forecasts)) )
  
}