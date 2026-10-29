#' NNS VAR
#'
#' Nonparametric vector autoregressive model incorporating \link{NNS.ARMA} estimates of variables into \link{NNS.reg} for a multi-variate time-series forecast.
#'
#' @param variables a numeric matrix, data.frame, \code{ts}, \code{xts} or \code{zoo} object of contemporaneous time-series to forecast.
#' @param h integer; 1 (default) Number of periods to forecast. \code{(h = 0)} will return just the interpolated and extrapolated values.
#' @param tau positive integer [ > 0]; 1 (default) Number of lagged observations to consider for the time-series data.  A single integer uses lags \code{1:tau} for every variable.  A vector with one entry per variable uses that single lag for each respective variable, and a list with one lag vector per variable (or one lag vector shared by every variable) uses multiple lags per variable.  The smallest lag of each variable is also the \code{modulo} for its \link{NNS.seas} seasonal-period search.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all") method for reducing regressors via \link{NNS.stack}.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for dimension reduction in the lagged variable matrix.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param naive.weights logical; \code{TRUE} (default) Equal weights applied to univariate and multivariate outputs in ensemble.  \code{FALSE} will apply weights based on the number of relevant variables detected.
#' @param obj.fn expression;
#' \code{expression(mean((predicted - actual)^2)) / (Sum of NNS Co-partial moments)} (default) MSE / co-movements is the default objective function.  Any \code{expression(...)} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.ARMA.optim}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param nowcast logical; \code{FALSE} (default) internal call for frequency alignment in downstream nowcasting applications.  When \code{TRUE}, rows of a time-indexed \code{variables} are labelled by month (\code{"Mon YYYY"}) whatever the input frequency, with the forecast rows labelled by the following \code{h} months.
#'
#' @return Returns the following matrices of forecasted variables:
#' \itemize{
#'  \item{\code{"interpolated_and_extrapolated"}} Returns a \code{data.frame} of the linear interpolated and \link{NNS.ARMA} extrapolated values to replace \code{NA} values in the original \code{variables} argument.  This is required for working with variables containing different frequencies, e.g. where \code{NA} would be reported for intra-quarterly data when indexed with monthly periods.
#'  \item{\code{"relevant_variables"}} Returns, for each variable, the lagged regressors selected by the \link{NNS.stack} dimension reduction, ranked by absolute coefficient.  These are exactly the regressors the \code{"multivariate"} estimate is built from.
#'
#'  \item{\code{"univariate"}} Returns the univariate \link{NNS.ARMA} forecasts.
#'
#'  \item{\code{"multivariate"}} Returns the multi-variate \link{NNS.reg} forecasts.
#'
#'  \item{\code{"ensemble"}} Returns the ensemble of both \code{"univariate"} and \code{"multivariate"} forecasts.
#'  }
#'
#' When \code{variables} is a \code{ts}, \code{xts} or \code{zoo} object, rows are labelled with its time index and the forecast rows with the following \code{h} periods.
#'
#' @note
#' \itemize{
#' \item \code{"Error in { : task xx failed -}"} should be re-run with \code{NNS.VAR(..., ncores = 1)}.
#' \item Not recommended for factor variables, even after transformed to numeric.  \link{NNS.reg} is better suited for factor or binary regressor extrapolation.
#' \item The lagged training sample, \code{nrow(variables) - max(tau)}, must contain at least \code{2 * h + 3} observations.
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
                    h = 1,
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
  default_obj_fn <- missing(obj.fn)

  # ===================== Argument validation =====================
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h < 0 || h != floor(h)) {
    stop("[h] must be a single non-negative integer.", call. = FALSE)
  }
  h <- as.integer(h)

  if (!is.character(objective) || length(objective) != 1L || is.na(objective) ||
      !(tolower(objective) %in% c("min", "max"))) {
    stop("[objective] must be 'min' or 'max'.", call. = FALSE)
  }
  objective <- tolower(objective)
  use_default_objective <- default_obj_fn && objective == "min"

  # dim.red.method is scalar: use "all" to combine methods.
  if (!is.character(dim.red.method) || length(dim.red.method) != 1L || is.na(dim.red.method) ||
      !(tolower(dim.red.method) %in% c("cor", "nns.dep", "nns.caus", "all"))) {
    stop('Please ensure the dimension reduction method is set to a single value: one of "cor", "nns.dep", "nns.caus" or "all".', call. = FALSE)
  }
  dim.red.method <- tolower(dim.red.method)

  .check_flag <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("[%s] must be TRUE or FALSE.", name), call. = FALSE)
    }
  }
  .check_flag(naive.weights, "naive.weights")
  .check_flag(status, "status")
  .check_flag(nowcast, "nowcast")

  # ===================== Inputs, names and time index =====================
  time_index <- .nns_var_time_index(variables)
  variables <- .nns_var_as_frame(variables)
  p <- ncol(variables)
  n <- nrow(variables)

  lag_spec <- .nns_var_lags(tau, p)

  if (h > 0L) {
    # NNS.stack holds out max(2h, 20%) of the lagged training rows and needs
    # at least three fitting rows ahead of that window.
    n_lagged <- n - lag_spec$max
    n_needed <- max(2L * h + 3L, 4L)
    if (n_lagged < n_needed) {
      stop(sprintf(paste0("Not enough observations: with h = %d and max(tau) = %d the lagged ",
                          "training sample has %d rows but needs at least %d."),
                   h, lag_spec$max, max(n_lagged, 0L), n_needed),
           call. = FALSE)
    }
  }

  labels <- .nns_var_row_labels(time_index, n, h, nowcast)

  # Variables never observed below zero keep every output non-negative.
  nonneg <- vapply(variables, function(v) min(v, na.rm = TRUE) >= 0, logical(1L))

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
    p >= 2L &&
    (!auto_cores ||
       p >= ncol_parallel_min ||
       n >= 500L)

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
    x <- variables[[i]]
    len <- length(x)
    observed <- which(!is.na(x))
    last_observed <- max(observed)
    filled <- x

    # Leading and interior NA(s): fit on the index and fill ONLY the missing
    # indices; observed values are never smoothed.
    inner_missing <- which(is.na(x) & seq_len(len) < last_observed)
    if (length(inner_missing)) {
      filled[inner_missing] <- .nns_var_index_fit(observed, x[observed], inner_missing)
    }

    # Trailing NA(s): NNS.ARMA extrapolation of the gap-filled history.
    if (last_observed < len) {
      trailing <- (last_observed + 1L):len
      filled[trailing] <- .nns_var_extrapolate(filled[seq_len(last_observed)], length(trailing),
                                               lag_spec$modulo[i], nonneg[i], obj.fn, objective)
    }

    if (nonneg[i]) filled <- pmax(0, filled)

    forecast <- if (h > 0L) {
      .nns_var_arma(filled, h, lag_spec$modulo[i], nonneg[i], obj.fn, objective)
    } else NULL

    list(filled, forecast)
  }

  nns_IVs <- if (use_parallel) {
    parallel::clusterExport(
      cl,
      varlist = c("variables", "h", "lag_spec", "nonneg", "obj.fn", "objective"),
      envir = environment()
    )
    parallel::parLapply(cl, seq_len(p), .interp_worker)
  } else {
    lapply(seq_len(p), .interp_worker)
  }

  .as_columns <- function(columns) {
    out <- as.data.frame(setNames(columns, colnames(variables)), check.names = FALSE)
    rownames(out) <- NULL
    out
  }

  interpolation_results <- lapply(nns_IVs, `[[`, 1L)
  nns_IVs_interpolated_extrapolated <- .as_columns(interpolation_results)
  rownames(nns_IVs_interpolated_extrapolated) <- labels$observed

  if(h == 0L) return(.NNS.df(nns_IVs_interpolated_extrapolated))

  extrapolation_results <- lapply(nns_IVs, `[[`, 2L)
  nns_IVs_results <- .as_columns(extrapolation_results)

  # Combine interpolated / extrapolated / forecasted IVs onto training data.frame
  new_values <- .as_columns(Map(c, interpolation_results, extrapolation_results))

  # Now lag new forecasted data.frame. The first nrow - h lagged rows are the
  # observed training sample; the last h rows are the forecast horizon, whose
  # regressors carry the univariate forecasts of every variable.
  lagged <- .nns_var_lag_matrix(new_values, lag_spec$lags)
  lagged_new_values <- lagged$values
  lagged_new_values_train <- head(lagged_new_values, nrow(lagged_new_values) - h)
  lagged_new_values_test <- tail(lagged_new_values, h)
  lag_source <- lagged$source


  if(status) message("Currently generating multi-variate estimates...", "\r", appendLF = TRUE)


  if(use_parallel){
    if(status) message("Parallel process running, status unavailable... \n","\r",appendLF=FALSE)
    status <- FALSE
  }


  .model_worker <- function(i) {
    if(status) message("Variable ", i, " of ", p, appendLF = TRUE)

    # Column i is variable i's contemporaneous (tau_0) column.
    IV <- lagged_new_values_train[, -i, drop = FALSE]
    DV <- lagged_new_values_train[[i]]
    IV.test <- lagged_new_values_test[, -i, drop = FALSE]

    # ts.test must be a whole number: NNS.stack validates it as an integer >= 1.
    ts <- max(2L * h, as.integer(ceiling(.2 * length(DV))))

    # The dimension reduction inside NNS.stack selects the regressors; with
    # stack = TRUE both of its methods are fitted on exactly that selection.
    .run_var_stack <- function(stack_obj_fn, stack_objective) {
      suppressWarnings(
        NNS.stack(IVs.train = IV,
                  DV.train = DV,
                  IVs.test = IV.test,
                  ts.test = ts,
                  folds = 1,
                  obj.fn = stack_obj_fn,
                  objective = stack_objective,
                  method = c(1, 2),
                  dim.red.method = dim.red.method,
                  order = NULL, ncores = 1, stack = TRUE, status = FALSE)
      )
    }

    fit <- tryCatch(
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

    # Report exactly the regressors the estimate was built from. Without a
    # dimension reduction step (a single regressor) every regressor is used.
    rel_vars <- fit$NNS.dim.red.variables
    if (is.null(rel_vars) || !length(rel_vars)) rel_vars <- names(IV)

    nns_DVs <- as.numeric(fit$stack)
    missing_estimate <- !is.finite(nns_DVs)
    nns_DVs[missing_estimate] <- nns_IVs_results[missing_estimate, i]
    if (nonneg[i]) nns_DVs <- pmax(0, nns_DVs)

    list(nns_DVs, rel_vars)
  }

  lists <- if (use_parallel) {
    parallel::clusterExport(
      cl,
      varlist = c("status", "p", "lagged_new_values_train", "lagged_new_values_test",
                  "h", "obj.fn", "objective", "use_default_objective",
                  "dim.red.method", "nns_IVs_results", "nonneg"),
      envir = environment()
    )
    parallel::parLapply(cl, seq_len(p), .model_worker)
  } else {
    lapply(seq_len(p), .model_worker)
  }

  if(use_parallel) {
    parallel::stopCluster(cl)
    cl <- NULL
    invisible(gc(verbose = FALSE))
  }

  nns_DVs <- .as_columns(lapply(lists, `[[`, 1L))
  relevant_vars <- lapply(lists, `[[`, 2L)

  RV <- do.call(cbind, lapply(relevant_vars, `length<-`, max(lengths(relevant_vars))))
  colnames(RV) <- colnames(variables)

  # Ensemble weights: the univariate share is the fraction of the relevant
  # regressors that are the variable's own lags. lag_source maps every lagged
  # column to its variable exactly, so variable names are never parsed.
  uni <- rep(0.5, p)
  if (!naive.weights) {
    for (i in seq_len(p)) {
      observed_var <- lag_source[match(relevant_vars[[i]], names(lagged_new_values))]
      observed_var <- observed_var[!is.na(observed_var)]
      if (length(observed_var)) uni[i] <- mean(observed_var == colnames(variables)[i])
    }
  }
  multi <- 1 - uni

  forecasts <- .as_columns(Map(function(u, m, a, b) u * a + m * b,
                               uni, multi, nns_IVs_results, nns_DVs))

  rownames(nns_IVs_results) <- labels$forecast
  rownames(nns_DVs) <- labels$forecast
  rownames(forecasts) <- labels$forecast

  return( list("interpolated_and_extrapolated" = .NNS.df(nns_IVs_interpolated_extrapolated),
               "relevant_variables" = .NNS.df(data.frame(RV, check.names = FALSE)),
               univariate = .NNS.df(nns_IVs_results),
               multivariate = .NNS.df(nns_DVs),
               ensemble = .NNS.df(forecasts)) )

}


# ===================== NNS.VAR internals =====================

# Coerce the supported inputs to a numeric data.frame with unique, non-empty
# column names.
.nns_var_as_frame <- function(variables) {
  if (any(class(variables) %in% c("tbl", "data.table"))) variables <- as.data.frame(variables)
  if (inherits(variables, c("xts", "zoo", "ts"))) {
    variables <- matrix(as.numeric(variables), nrow = NROW(variables),
                        dimnames = list(NULL, colnames(variables)))
  }
  if (is.null(dim(variables))) {
    if (!is.numeric(variables) && !is.logical(variables)) {
      stop("[variables] must be a numeric matrix or data.frame of time-series.", call. = FALSE)
    }
    variables <- matrix(as.numeric(variables), ncol = 1L)
  }

  var_names <- colnames(variables)
  variables <- as.data.frame(variables, check.names = FALSE, stringsAsFactors = FALSE)

  if (ncol(variables) < 1L || nrow(variables) < 1L) {
    stop("[variables] must contain at least one variable and one observation.", call. = FALSE)
  }

  numeric_cols <- vapply(variables, function(v) is.numeric(v) || is.logical(v), logical(1L))
  if (!all(numeric_cols)) {
    stop(sprintf("All [variables] must be numeric; non-numeric column(s): %s.",
                 paste(which(!numeric_cols), collapse = ", ")), call. = FALSE)
  }
  variables[] <- lapply(variables, as.numeric)
  if (any(vapply(variables, function(v) any(is.infinite(v)), logical(1L)))) {
    stop("[variables] contains infinite values.", call. = FALSE)
  }
  empty <- vapply(variables, function(v) all(is.na(v)), logical(1L))
  if (any(empty)) {
    stop(sprintf("Variable(s) %s contain only missing values.",
                 paste(which(empty), collapse = ", ")), call. = FALSE)
  }

  if (is.null(var_names)) var_names <- rep("", ncol(variables))
  blank <- is.na(var_names) | var_names == ""
  var_names[blank] <- paste0("x", which(blank))
  var_names <- gsub(" - ", "...", var_names)
  colnames(variables) <- make.unique(var_names)
  rownames(variables) <- NULL
  variables
}

# Validate tau and expand it to one lag vector per variable. `modulo` is the
# NNS.seas modulo for each variable: its smallest lag (a scalar tau keeps tau).
.nns_var_lags <- function(tau, p) {
  invalid <- function() {
    stop(paste0("[tau] must be a positive integer, a vector of positive integers with one lag ",
                "per variable, or a list of positive-integer lag vectors (one per variable, or ",
                "one shared by every variable)."), call. = FALSE)
  }
  as_lags <- function(v) {
    v <- unlist(v, use.names = FALSE)
    if (!is.numeric(v) || !length(v) || any(!is.finite(v)) || any(v < 1) || any(v != floor(v))) {
      invalid()
    }
    sort(unique(as.integer(v)))
  }

  if (is.list(tau)) {
    if (length(tau) == 1L) tau <- rep(tau, p)
    if (length(tau) != p) invalid()
    lags <- lapply(tau, as_lags)
    modulo <- vapply(lags, min, integer(1L))
  } else {
    v <- as_lags(tau)
    if (length(tau) == 1L) {
      lags <- rep(list(seq_len(v)), p)
      modulo <- rep(v, p)
    } else if (length(tau) == p) {
      lags <- as.list(as.integer(tau))
      modulo <- as.integer(tau)
    } else {
      invalid()
    }
  }

  list(lags = lags, modulo = as.integer(modulo), max = max(unlist(lags)))
}

# Lag matrix: every variable's contemporaneous (tau_0) column first, in
# variable order, then each variable's requested lags. `source` names the
# variable behind every column.
.nns_var_lag_matrix <- function(x, lags) {
  x <- as.matrix(x)
  var_names <- gsub('"', '', colnames(x), fixed = TRUE)
  max_tau <- max(unlist(lags))
  rows <- (max_tau + 1L):nrow(x)

  columns <- list(); col_names <- character(); col_source <- character()
  for (j in seq_len(ncol(x))) {
    columns[[length(columns) + 1L]] <- x[rows, j]
    col_names <- c(col_names, paste0(var_names[j], "_tau_0"))
    col_source <- c(col_source, colnames(x)[j])
  }
  for (j in seq_len(ncol(x))) {
    for (l in lags[[j]]) {
      columns[[length(columns) + 1L]] <- x[rows - l, j]
      col_names <- c(col_names, paste0(var_names[j], "_tau_", l))
      col_source <- c(col_source, colnames(x)[j])
    }
  }

  values <- as.data.frame(do.call(cbind, columns), check.names = FALSE)
  names(values) <- make.unique(col_names)
  list(values = values, source = col_source)
}

# Fit a series on its index and evaluate it at `at` (interpolation inside the
# observed range, linear extrapolation outside it).
.nns_var_index_fit <- function(index, values, at) {
  if (length(index) < 2L) return(rep(values[1L], length(at)))
  as.numeric(suppressWarnings(
    NNS.reg(index, values, order = "max", ncores = 1, point.est = at,
            plot = FALSE, point.only = TRUE)$Point.est
  ))
}

# NNS.ARMA.optim forecast with a seasonal-period fallback chain: the modulo
# periods from NNS.seas, then the unrestricted NNS.seas periods, then a plain
# lag-1 model. NNS.ARMA.optim rejects periods too long for its training set,
# so one unusable candidate set no longer aborts the whole call.
.nns_var_arma <- function(series, horizon, modulo, nonneg, obj.fn, objective) {
  seas_periods <- function(...) {
    periods <- tryCatch(suppressWarnings(NNS.seas(series, mod.only = FALSE, plot = FALSE, ...)$periods),
                        error = function(e) NULL)
    periods <- suppressWarnings(as.integer(periods))
    periods[!is.na(periods) & periods >= 1L]
  }
  candidates <- list(function() seas_periods(modulo = modulo),
                     function() seas_periods(),
                     function() 1L)

  failure <- NULL
  for (candidate in candidates) {
    periods <- candidate()
    if (!length(periods)) next
    fit <- tryCatch(
      suppressWarnings(NNS.ARMA.optim(series, h = horizon, seasonal.factor = periods,
                                      obj.fn = obj.fn, objective = objective,
                                      print.trace = FALSE, ncores = 1,
                                      negative.values = !nonneg)),
      error = function(e) {
        failure <<- e
        NULL
      }
    )
    if (!is.null(fit)) {
      forecast <- as.numeric(fit$results)
      if (nonneg) forecast <- pmax(0, forecast)
      return(forecast)
    }
  }
  stop(paste0("NNS.ARMA.optim failed for every seasonal-period candidate: ",
              if (is.null(failure)) "no candidate periods" else conditionMessage(failure)),
       call. = FALSE)
}

# Trailing NA(s): NNS.ARMA extrapolation, falling back to a linear index fit
# when the observed history is too short for NNS.ARMA.optim.
.nns_var_extrapolate <- function(history, horizon, modulo, nonneg, obj.fn, objective) {
  tryCatch(
    .nns_var_arma(history, horizon, modulo, nonneg, obj.fn, objective),
    error = function(e) {
      out <- .nns_var_index_fit(seq_along(history), history, length(history) + seq_len(horizon))
      if (nonneg) pmax(0, out) else out
    }
  )
}

# Best-effort time index of a ts / xts / zoo input, without zoo or xts.
.nns_var_time_index <- function(variables) {
  if (inherits(variables, "ts")) {
    return(list(type = "time", time = as.numeric(stats::time(variables)),
                frequency = stats::frequency(variables)))
  }
  if (!inherits(variables, c("xts", "zoo"))) return(NULL)

  raw <- attr(variables, "index")
  if (is.null(raw)) return(NULL)

  if (inherits(variables, "xts")) {
    # xts stores its index as numeric seconds since the epoch.
    tz <- attr(raw, "tzone")
    if (is.null(tz) || !nzchar(tz[1L])) tz <- "UTC"
    raw <- as.POSIXct(as.numeric(raw), origin = "1970-01-01", tz = tz[1L])
  }

  if (inherits(raw, "yearmon")) {
    return(list(type = "time", time = as.numeric(raw), frequency = 12))
  }
  if (inherits(raw, "yearqtr")) {
    return(list(type = "time", time = as.numeric(raw), frequency = 4))
  }
  if (inherits(raw, "POSIXt")) {
    lt <- as.POSIXlt(raw)
    if (all(lt$hour == 0 & lt$min == 0 & lt$sec == 0)) {
      return(list(type = "date", date = as.Date(format(raw, "%Y-%m-%d"))))
    }
    return(list(type = "datetime", datetime = as.POSIXct(raw)))
  }
  if (inherits(raw, "Date")) return(list(type = "date", date = raw))
  if (is.numeric(raw) && length(raw) >= 2L) {
    spacing <- unique(round(diff(as.numeric(raw)), 8))
    if (length(spacing) == 1L && spacing > 0) {
      return(list(type = "time", time = as.numeric(raw), frequency = 1 / spacing))
    }
  }
  NULL
}

# "Mon YYYY" for monthly, "YYYY Qn" for quarterly, "YYYY" for annual and
# "YYYY(cycle)" for any other frequency.
.nns_var_cycle_labels <- function(time, frequency) {
  year <- floor(time + 1e-8)
  cycle <- as.integer(round((time - year) * frequency)) + 1L
  roll <- !is.na(cycle) & cycle > round(frequency)
  year[roll] <- year[roll] + 1; cycle[roll] <- 1L
  if (isTRUE(all.equal(frequency, 12))) return(paste(month.abb[cycle], year))
  if (isTRUE(all.equal(frequency, 4))) return(paste0(year, " Q", cycle))
  if (isTRUE(all.equal(frequency, 1))) return(as.character(year))
  paste0(year, "(", cycle, ")")
}

# Row labels for the observed rows and the h forecast rows. Labels that are not
# unique are dropped (with a warning) rather than failing the call.
.nns_var_row_labels <- function(time_index, n, h, nowcast) {
  none <- list(observed = NULL, forecast = NULL)
  if (is.null(time_index)) return(none)
  step <- seq_len(h)

  as_months <- function(date) {
    lt <- as.POSIXlt(date)
    (lt$year + 1900) + lt$mon / 12
  }
  monthly_labels <- function(months) {
    .nns_var_cycle_labels(c(months, tail(months, 1L) + step / 12), 12)
  }

  labels <- switch(
    time_index$type,
    time = {
      tt <- time_index$time
      f <- if (nowcast) 12 else time_index$frequency
      .nns_var_cycle_labels(c(tt, tail(tt, 1L) + step / f), f)
    },
    date = {
      d <- time_index$date
      spacing <- as.numeric(diff(d))
      months <- as_months(d)
      if (nowcast || (length(d) >= 2L && all(spacing >= 28 & spacing <= 31))) {
        monthly_labels(months)
      } else if (length(d) >= 2L && all(spacing >= 89 & spacing <= 92)) {
        quarters <- floor(months + 1e-8) + floor(round((months %% 1) * 12) / 3) / 4
        .nns_var_cycle_labels(c(quarters, tail(quarters, 1L) + step / 4), 4)
      } else if (length(d) >= 2L && all(spacing >= 365 & spacing <= 366)) {
        years <- floor(months + 1e-8)
        as.character(c(years, tail(years, 1L) + step))
      } else {
        last <- tail(d, 1L)
        weekdays_only <- length(d) >= 2L && !any(as.POSIXlt(d)$wday %in% c(0L, 6L))
        future <- if (length(d) < 2L) {
          last + step
        } else if (weekdays_only && stats::median(spacing) <= 1) {
          days <- last + seq_len(2L * h + 7L)
          head(days[!(as.POSIXlt(days)$wday %in% c(0L, 6L))], h)
        } else {
          last + stats::median(spacing) * step
        }
        format(c(d, future), "%Y-%m-%d")
      }
    },
    datetime = {
      dt <- time_index$datetime
      if (nowcast) {
        monthly_labels(as_months(dt))
      } else {
        spacing <- if (length(dt) >= 2L) stats::median(as.numeric(diff(dt), units = "secs")) else 1
        format(c(dt, tail(dt, 1L) + spacing * step), "%Y-%m-%d %H:%M:%S")
      }
    },
    NULL
  )

  if (is.null(labels) || length(labels) != n + h || anyNA(labels) || anyDuplicated(labels)) {
    if (!is.null(labels)) {
      warning("NNS.VAR: the time index does not give unique row labels; rows are left unlabelled.",
              call. = FALSE)
    }
    return(none)
  }
  list(observed = labels[seq_len(n)], forecast = labels[n + step])
}
