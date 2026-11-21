.nns_has_factor <- function(x) {
  is.factor(x) || (is.list(x) && !inherits(x, "POSIXlt") &&
                     any(vapply(x, .nns_has_factor, logical(1))))
}

.nns_no_factor <- function(x, arg) {
  if (.nns_has_factor(x)) stop(sprintf("`%s` must be numeric, not a factor.", arg))
}

.nns_finite <- function(x, arg) {
  if (!all(is.finite(x))) stop(sprintf("`%s` must be finite.", arg))
  x
}

.nns_sample <- function(x, arg, na = "keep", inf = "allow") {
  if (is.double(x) && is.null(attributes(x)) && na == "keep" && inf == "allow") return(x)
  if (is.factor(x)) .nns_no_factor(x, arg)
  value <- as.numeric(x)
  if (is.list(x)) .nns_no_factor(x, arg)
  if (na == "omit") value <- value[!is.na(value)]
  if ((na == "reject" && anyNA(value)) || (inf == "reject" && any(is.infinite(value))))
    stop(sprintf("`%s` must be finite.", arg))
  value
}

.nns_var_sample <- function(x, arg) {
  .nns_no_factor(x, arg)
  dated <- is.data.frame(x) && ncol(x) > 1L &&
    any(vapply(x, inherits, logical(1), c("Date", "POSIXt")))
  if (inherits(x, c("tbl", "data.table"))) x <- .nns_unlist_measures(x)
  value <- .nns_sample(x, arg, na = "omit")
  if (!length(value)) stop("x must contain at least one non-NA numeric value.")
  if (dated) stop(sprintf("`%s` must have exactly one column when a date is present.", arg))
  value
}

# Existing callers rely on the native conversion diagnostic.
.nns_native_type_error <- function(x) {
  stop(sprintf("Not compatible with requested type: [type=%s; target=double].", typeof(x)))
}

.nns_native_sample <- function(x, arg) {
  if (is.data.frame(x)) {
    flat <- as.vector(.nns_unlist_measures(x))
    if (!typeof(flat) %in% c("double", "integer", "logical", "complex", "raw"))
      .nns_native_type_error(flat)
    .nns_no_factor(x, arg)
    x <- flat
  } else if (!typeof(x) %in% c("double", "integer")) {
    stop(sprintf("%s should be numeric vector, or data table", arg))
  }
  .nns_sample(x, arg)
}

.nns_sd_column <- function(x, arg) {
  if ((is.data.frame(x) || is.matrix(x) || inherits(x, c("xts", "zoo"))) &&
      !is.null(dim(x)) && ncol(x) > 1L)
    stop(sprintf("%s must be a single-column object or numeric vector.", arg))
  if (is.factor(x)) stop(sprintf("%s must be numeric, not a factor.", arg))
  if (is.data.frame(x) || is.matrix(x) || inherits(x, c("xts", "zoo")))
    x <- as.vector(unlist(x, use.names = FALSE))
  .nns_sample(x, arg)
}

.nns_numeric_labels <- function(x, arg) {
  if (is.data.frame(x)) {
    if (ncol(x) != 1L)
      stop(sprintf("%s must be a vector or single-column object in the bivariate path.", arg))
    x <- x[[1L]]
  }
  if (is.factor(x)) x <- as.character(x)
  .nns_sample(x, arg)
}

.nns_is_measure_column <- function(x) {
  is.object(x) && is.numeric(x) && is.double(x) && is.null(dim(x))
}

.nns_measure_column <- function(x) {
  if (.nns_is_measure_column(x) && !inherits(x, c("ts", "zoo"))) {
    value <- as.numeric(x)
    names(value) <- names(x)
    return(value)
  }
  x
}

.nns_measure_frame <- function(x) {
  if (!is.data.frame(x)) return(x)
  columns <- vapply(x, .nns_is_measure_column, logical(1))
  if (any(columns)) {
    x <- as.data.frame(x)
    x[columns] <- lapply(x[columns], as.numeric)
  }
  x
}

.nns_sample_matrix <- function(x, arg, finite = FALSE) {
  if (!is.matrix(x) || finite) x <- as.matrix(.nns_measure_frame(x))
  if (finite) .nns_finite(x, arg) else x
}

.nns_copula_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.double(x)) storage.mode(x) <- "double"
  x
}

.nns_target <- function(x, arg = "target", native = FALSE, means = FALSE) {
  if (is.null(x) || (means && (identical(x, "mean") || identical(x, "means"))))
    return(NULL)
  if (native && (is.factor(x) || !typeof(x) %in% c("double", "integer", "logical")))
    stop("`target` must be numeric or NULL.")
  .nns_sample(x, arg)
}

.nns_thresholds <- function(x, arg = "target") {
  .nns_sample(x, arg, na = "reject", inf = "reject")
}

.nns_degree <- function(x, arg = "degree", form = "native") {
  .nns_no_factor(x, arg)
  if (form == "quantile") {
    if (is.object(x) && is.double(x) && is.numeric(x)) x <- as.numeric(x)
    if (x == 0) return(0)
    # Only these four integer degrees have a compiled quantile solver.
    if (length(x) == 1L && is.finite(x) && x == as.integer(x) && x >= 1 && x <= 4)
      return(as.integer(x))
  }
  if (form == "coerce") return(as.numeric(x))
  if (length(x) != 1L) stop(sprintf("Expecting a single value: [extent=%s].", length(x)))
  if (!typeof(x) %in% c("integer", "double", "logical", "raw", "complex"))
    .nns_native_type_error(x)
  as.numeric(x)
}

.nns_probability <- function(x) {
  pmin(pmax(.nns_sample(x, "percentile"), 0), 1)
}

.nns_column <- function(x, arg) {
  if (is.data.frame(x) || is.matrix(x)) {
    if (NCOL(x) != 1L)
      stop(sprintf("[%s] must be a vector or one-column object.", arg), call. = FALSE)
    if (is.data.frame(x)) x <- x[[1L]] else x <- x[, 1L]
  }
  x
}

.nns_moment_pair <- function(x, y) {
  list(x = .nns_native_sample(x, "x"), y = .nns_native_sample(y, "y"))
}

.nns_dependence_column <- function(x, arg) {
  if (is.data.frame(x) && ncol(x) == 1L) x <- x[[1L]]
  if (is.factor(x)) {
    labels <- suppressWarnings(as.numeric(levels(x)))[as.integer(x)]
    if (any(is.na(labels) & !is.na(x))) return(as.numeric(x))
    return(labels)
  }
  .nns_numeric_labels(x, arg)
}

.nns_dependence_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.data.frame(x)
    factors <- vapply(x, is.factor, logical(1))
    x[factors] <- lapply(x[factors], .nns_dependence_column, arg = "x")
    x <- data.matrix(.nns_measure_frame(x))
  }
  if (anyNA(x)) stop("x has missing values, please address.")
  if (any(is.infinite(x))) stop("x must be finite.")
  x
}

.nns_dependence_pair <- function(x, y) {
  x <- .nns_dependence_column(x, "x")
  y <- .nns_dependence_column(y, "y")
  if (anyNA(x)) stop("x has missing values, please address.")
  if (anyNA(y)) stop("y has missing values, please address.")
  if (any(is.infinite(x))) stop("x must be finite.")
  if (any(is.infinite(y))) stop("y must be finite.")
  if (length(x) != length(y)) stop("x and y must have the same length.")
  list(x = x, y = y)
}

.nns_partition_pair <- function(x, y) {
  x <- .nns_column(x, "x")
  y <- .nns_column(y, "y")
  if (length(x) != length(y)) stop("[x] and [y] must have the same length.", call. = FALSE)
  if (!length(x)) stop("[x] and [y] must not be empty.", call. = FALSE)
  # Partition predictors use their level ordering.
  x <- if (is.factor(x)) as.numeric(x) else .nns_sample(x, "x")
  y <- .nns_sample(y, "y")
  keep <- !(is.na(x) | is.na(y))
  if (!any(keep))
    stop("[x] and [y] must contain at least one complete (non-missing) pair.", call. = FALSE)
  x <- x[keep]
  y <- y[keep]
  list(x = x, y = y)
}

.nns_series_values <- function(x, arg) {
  if (!is.null(dim(x)) && NCOL(x) != 1L)
    stop(sprintf("`%s` must have exactly one column.", arg))
  if (inherits(x, c("tbl", "data.table"))) {
    .nns_no_factor(x, arg)
    x <- as.vector(.nns_unlist_measures(x))
  }
  .nns_no_factor(x, arg)
  if (anyNA(x)) stop("You have some missing values, please address.")
  x
}

.nns_bootstrap_series <- function(x, arg) {
  if (is.data.frame(x) && ncol(x) == 1L) x <- x[[1L]]
  .nns_series_values(.nns_measure_column(x), arg)
}

.nns_forecast_series <- function(x, arg) {
  as.numeric(.nns_series_values(x, arg))
}

.nns_seasonal_series <- function(x, arg) {
  if (!inherits(x, c("tbl", "data.table")) && !is.numeric(x))
    stop("Variable must be numeric")
  x <- .nns_series_values(x, arg)
  if (!is.numeric(x)) stop("Variable must be numeric")
  if (any(is.infinite(x))) stop("Infinite values not allowed")
  as.numeric(x)
}

.nns_observed <- function(x, rows = FALSE) {
  if (inherits(x, c("tbl", "data.table"))) {
    if (!rows || ncol(x) == 1L) x <- .nns_table_values(x) else x <- as.data.frame(x)
  }
  if (rows && is.factor(x)) stop("`variable` must be numeric, not a factor.")
  if (!rows || is.null(dim(x))) x[!is.na(x)] else x[stats::complete.cases(x), , drop = FALSE]
}

.nns_point <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("'point' must be a single finite numeric value.")
  x
}

.nns_coordinates <- function(x, arg) {
  .nns_no_factor(x, arg)
  x <- if (is.list(x) && !inherits(x, "POSIXlt")) .nns_unlist_measures(x) else unlist(x)
  values <- as.numeric(x)
  names(values) <- names(x)
  values
}

.nns_table_values <- function(x) {
  if (inherits(x, c("tbl", "data.table"))) as.vector(.nns_unlist_measures(x)) else .nns_measure_column(x)
}

.nns_evaluation_values <- function(x) {
  if (is.list(x)) x <- .nns_unlist_measures(x)
  .nns_sample(x, "eval.points")
}

.nns_derivative_x <- function(x) {
  x <- .nns_table_values(x)
  if (!is.null(ncol(x)) && is.null(colnames(x))) x <- unlist(data.frame(x))
  x
}

.nns_derivative_point <- function(x) {
  if (!length(x)) stop("eval.point must be provided: numeric point(s) or 'overall'")
  if (identical(unname(x), "overall")) return("overall")
  if (!is.numeric(x)) stop("eval.point must be numeric point(s) or 'overall'")
  value <- as.numeric(x)
  if (any(!is.finite(value))) stop("eval.point must contain only finite numeric values")
  value
}

.nns_unlist_measures <- function(x) {
  unlist(lapply(x, .nns_measure_column))
}

.nns_causation_inputs <- function(x, y, factor.2.dummy, tau) {
  if(!is.null(y))  if(anyNA(x) || anyNA(y)) stop("You have some missing values, please address.")
  if(is.null(y))  if(anyNA(x)) stop("You have some missing values, please address.")

  if(is.data.frame(x) && ncol(x) == 1) x <- x[[1L]]
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
  if(is.data.frame(y) && ncol(y) == 1) y <- y[[1L]]

  if(factor.2.dummy){
    if(!is.null(dim(x))){
      if(!is.numeric(x)){
        x <- do.call(cbind, lapply(x, factor_2_dummy_FR))
      } else {
        x <- apply(x, 2, as.double)
      }
      if(is.list(x)){
        x <- do.call(cbind, x)
        x <- apply(x, 2, as.double)
      }
    } else {
      x <- factor_2_dummy(x)
      if(is.null(dim(x))){
        x <- as.double(x)
      } else {
        x <- apply(x, 2, as.double)
      }
    }
  }

  if (is.null(y)) return(list(x = x, y = NULL))
  if(is.factor(y)) y <- as.numeric(y)

  # ---- input validation ---------------------------------------------------
  if(NCOL(x) > 1 || NCOL(y) > 1)
    stop("x and y must each be a single variable (factor.2.dummy = TRUE turns a factor with k levels into k - 1 columns).")
  if(length(x) != length(y)) stop("x and y must have the same length.")
  if(any(is.infinite(x)) || any(is.infinite(y))) stop("x and y must be finite.")
  if(length(unique(x)) < 2 || length(unique(y)) < 2)
    stop("x and y must not be constant: causation is undefined for a variable with no variation.")
  if(!(identical(tau, "cs") || identical(tau, "ts"))){
    if(!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau < 0 ||
       tau != round(tau) || tau >= length(x) - 2)
      stop("tau must be \"cs\", \"ts\" or a single non-negative integer smaller than length(x) - 2.")
  }

  x <- if (is.factor(x)) as.numeric(x) else .nns_sample(x, "x")
  y <- .nns_sample(y, "y")

  list(x = x, y = y)
}

.nns_distance_response <- function(x, class) {
  if (!is.null(class) && is.factor(x)) return(as.numeric(x))
  .nns_sample(x, "y.hat")
}

.nns_normalization_input <- function(X) {
  if(anyNA(X, recursive = TRUE)) stop("You have some missing values, please address.")

  if(any(class(X)%in%c("tbl","data.table"))) X <- as.data.frame(X)
  if(any(class(X)%in%"list")) return(X)
  columns <- apply(.nns_measure_frame(X), 2, unlist)
  if(is.null(dim(columns)) && nrow(X) == 1L){
    columns <- matrix(columns, nrow = 1L, dimnames = list(rownames(X), names(columns)))
  }
  columns
}

.nns_superiority_sample <- function(x, arg) {
  if (is.factor(x) || (is.list(x) && any(vapply(x, is.factor, logical(1)))))
    stop(sprintf("`%s` must be numeric, not a factor.", arg))
  x <- stats::na.omit(x)
  value <- if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x)
  if (anyNA(value)) stop(sprintf("`%s` must contain only numeric values.", arg))
  value
}

.nns_forecast_method <- function(method) {
  method <- tolower(method)
  if(length(method) != 1L || !method %in% c("lin", "nonlin", "both", "means"))
    stop("[method] must be one of 'lin', 'nonlin', 'both' or 'means'.")
  method
}

.arma_check_horizon <- function(h){
  if(!is.numeric(h) || length(h) != 1L || !is.finite(h) || h < 1 || h != round(h))
    stop("[h] must be a positive integer.")
}

.arma_check_periods <- function(periods){
  if(0 %in% periods)
    stop("NNS.seas returns 0 when the history is too short to scan for a period. ",
         "Pass a positive lag supported by the training data, e.g. seasonal.factor = 1.")
}

.nns_forecast_objective <- function(objective) {
  objective <- tolower(objective)
  if(length(objective) != 1L || !objective %in% c("min", "max"))
    stop("[objective] must be 'min' or 'max'.")
  objective
}
