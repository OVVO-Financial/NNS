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

.nns_predictor_levels <- function(x) {
  if (is.ordered(x)) return(levels(x))
  if (is.factor(x)) return(sort(levels(x)))
  sort(unique(as.character(x)))
}

.nns_predictor_frames <- function(train, test) {
  schema <- .nns_predictor_schema(train)
  for (j in seq_along(train)) {
    column <- schema[[j]]
    nm <- column$name
    tr <- train[[j]]
    te <- test[[j]]

    if (column$kind == "Date") {
      if (!inherits(te, "Date")) {
        stop(sprintf("Test predictor [%s] must also be a Date.", nm),
             call. = FALSE)
      }
      train[[j]] <- as.numeric(tr)
      test[[j]] <- as.numeric(te)
    } else if (column$kind == "datetime") {
      if (!inherits(te, c("POSIXct", "POSIXlt"))) {
        stop(sprintf("Test predictor [%s] must also be a date-time value.", nm),
             call. = FALSE)
      }
      train[[j]] <- as.numeric(tr)
      test[[j]] <- as.numeric(te)
    } else if (column$categorical && column$kind != "logical") {
      tr_chr <- as.character(tr)
      te_chr <- as.character(te)
      lev <- column$levels
      unseen <- setdiff(unique(te_chr), lev)
      if (length(unseen)) {
        stop(sprintf(
          "Test predictor [%s] contains unseen level(s): %s.",
          nm, paste(unseen, collapse = ", ")
        ), call. = FALSE)
      }
      if (is.factor(tr)) train[[j]] <- factor(tr_chr, levels = lev, ordered = is.ordered(tr))
      test[[j]] <- factor(te_chr, levels = lev, ordered = is.ordered(tr))
    } else if (column$kind == "logical") {
      if (!is.logical(te)) {
        stop(sprintf("Test predictor [%s] must also be logical.", nm),
             call. = FALSE)
      }
    } else if (is.numeric(tr) || is.integer(tr)) {
      if (!(is.numeric(te) || is.integer(te))) {
        stop(sprintf("Test predictor [%s] must be numeric.", nm),
             call. = FALSE)
      }
      train[[j]] <- as.numeric(tr)
      test[[j]] <- as.numeric(te)
    } else {
      stop(sprintf("Unsupported predictor type for [%s].", nm),
           call. = FALSE)
    }
  }
  list(train = train, test = test, metadata = schema)
}

.nns_check_predictors <- function(x, name) {
  if (anyNA(x)) stop(sprintf("[%s] contains missing values.", name), call. = FALSE)
  for (j in seq_along(x)) {
    if (is.numeric(x[[j]]) && any(!is.finite(x[[j]]))) {
      stop(sprintf("[%s] predictor [%s] contains non-finite values.",
                   name, names(x)[j]), call. = FALSE)
    }
  }
}

.nns_declare_characters <- function(x) {
  chars <- vapply(x, is.character, logical(1))
  x[chars] <- lapply(x[chars], factor)
  x
}

.nns_predictor_schema <- function(train) {
  unname(Map(function(x, name) {
    kind <- if (is.ordered(x)) "ordered" else if (is.factor(x)) "factor" else
      if (is.character(x)) "character" else if (is.logical(x)) "logical" else
      if (inherits(x, "Date")) "Date" else if (inherits(x, "POSIXt")) "datetime" else "numeric"
    categorical <- kind %in% c("ordered", "factor", "character", "logical")
    list(name = name, kind = kind, categorical = categorical,
         levels = if (categorical) .nns_predictor_levels(x))
  }, train, names(train)))
}

.nns_stack_layout <- function(schema) {
  blocks <- lapply(schema, function(col) {
    if (col$categorical && col$kind != "logical")
      paste0(col$name, "__", make.names(col$levels, unique = TRUE))
    else col$name
  })
  names <- make.unique(unlist(blocks, use.names = FALSE), sep = "_")
  list(names = names, source = rep(vapply(schema, `[[`, "", "name"), lengths(blocks)))
}

.nns_design_matrix <- function(frame, schema, dummy, label) {
  if (is.null(frame)) return(NULL)
  blocks <- Map(function(x, col, expand) {
    nm <- col$name
    if (!col$categorical) {
      value <- if (label == "point.est") suppressWarnings(.nns_numeric_labels(x, label)) else
        suppressWarnings(.nns_sample(x, nm))
      if (length(value) != length(x) || any(!is.finite(value)))
        stop(sprintf("%s '%s' must contain only finite numeric values.",
                     if (label == "point.est") "[point.est] predictor" else "Predictor", nm), call. = FALSE)
    } else {
      labels <- as.character(x)
      prefix <- if (label == "point.est") "[point.est] predictor" else "Predictor"
      if (anyNA(labels)) stop(sprintf("%s '%s' contains missing values.", prefix, nm), call. = FALSE)
      code <- match(labels, col$levels)
      if (anyNA(code)) stop(sprintf("%s '%s' contains unseen level(s): %s", prefix, nm,
                                    paste(unique(labels[is.na(code)]), collapse = ", ")), call. = FALSE)
      if (expand) {
        block <- outer(code, seq_along(col$levels), "==") * 1
        colnames(block) <- paste0(nm, "_", make.names(col$levels, unique = TRUE))
        return(block)
      }
      value <- if (is.null(col$values)) code else col$values[code]
    }
    matrix(as.numeric(value), ncol = 1L, dimnames = list(NULL, nm))
  }, frame, schema, dummy)
  do.call(cbind, unname(blocks))
}

.nns_encode_predictors <- function(train, test, dummy, schema) {
  list(x = .nns_design_matrix(train, schema, dummy, "x"),
       point.est = .nns_design_matrix(test, schema, dummy, "point.est"),
       metadata = schema)
}

.nns_stack_predictors <- function(train, test, schema = .nns_predictor_schema(train)) {
  layout <- .nns_stack_layout(schema)
  logicals <- vapply(schema, function(col) col$kind == "logical", logical(1))
  for (i in which(logicals)) {
    schema[[i]]$levels <- c("FALSE", "TRUE")
    schema[[i]]$values <- c(0, 1)
  }
  design <- .nns_encode_predictors(train, test, !logicals, schema)
  colnames(design$x) <- layout$names
  if (!is.null(design$point.est)) colnames(design$point.est) <- layout$names
  design$source <- layout$source
  design
}

.nns_reg_encode_predictors <- function(x, point.est = NULL,
                                       factor.2.dummy = TRUE, coefficient.names = NULL) {
  factor.2.dummy <- .nns_reg_scalar_logical(factor.2.dummy, "factor.2.dummy")
  # A plain numeric vector is the single numeric predictor "x", as
  # data.frame(x = x) below names it; building that data frame was most of the
  # encoding time for short series. Plain numeric points are its values.
  if (is.numeric(x) && is.null(dim(x)) && is.null(oldClass(x)) &&
      length(x) && all(is.finite(x)) &&
      (is.null(point.est) ||
         (is.numeric(point.est) && is.null(dim(point.est)) &&
            is.null(oldClass(point.est))))) {
    point.matrix <- NULL
    if (!is.null(point.est)) {
      if (any(!is.finite(point.est))) {
        stop("[point.est] predictor 'x' must contain only finite numeric values.",
             call. = FALSE)
      }
      point.matrix <- matrix(.nns_sample(point.est, "point.est"), ncol = 1L,
                             dimnames = list(NULL, "x"))
    }
    return(list(
      x = matrix(.nns_sample(x, "x"), ncol = 1L, dimnames = list(NULL, "x")),
      point.est = point.matrix,
      stack.index = FALSE,
      raw.names = "x",
      encoded.names = "x",
      metadata = list(list(name = "x", kind = "numeric", categorical = FALSE, levels = NULL))
    ))
  }
  train <- .nns_reg_as_frame(x)
  points <- .nns_reg_prepare_points(point.est, names(train))

  schema <- .nns_predictor_schema(train)
  stack.index <- FALSE
  if (!is.null(coefficient.names) && factor.2.dummy && ncol(train) > 1L &&
      any(vapply(schema, `[[`, logical(1), "categorical"))) {
    stack.names <- .nns_stack_layout(schema)$names
    stack.index <- length(coefficient.names) == length(stack.names) && setequal(coefficient.names, stack.names)
  }
  design <- if (stack.index) .nns_stack_predictors(train, points, schema) else
    .nns_encode_predictors(train, points, factor.2.dummy, schema)
  list(x = design$x, point.est = design$point.est, stack.index = stack.index,
       raw.names = names(train), encoded.names = colnames(design$x), metadata = design$metadata)
}

.nns_encode_response <- function(y, type) {
  is.class <- identical(type, "class")
  nominal <- is.factor(y) || is.character(y)
  if (nominal) {
    if (is.factor(y)) {
      labels <- levels(y)
      # An explicit NA level has a code distinct from a missing observation.
      value <- as.numeric(match(levels(y), labels)[as.integer(y)])
    } else {
      labels <- sort(unique(y))
      value <- as.numeric(match(y, labels))
    }
  } else {
    value <- .nns_sample(y, "y")
    labels <- NULL
  }
  values <- codes <- NULL
  if (is.class || nominal || is.logical(y)) {
    values <- sort(unique(value))
    codes <- as.numeric(match(value, values))
    if (is.class && is.logical(y)) values <- as.logical(values)
    if (!nominal) labels <- values
  }
  list(type = type, is.class = is.class, is.xonly = identical(type, "xonly"),
       y = value, codes = codes, class.levels = labels, class.values = values,
       labels = if (nominal) labels[values] else as.character(values),
       decode = function(code) values[code])
}

.nns_response_vector <- function(y, name = "y") {
  if (!is.data.frame(y) && !is.matrix(y)) return(y)
  if (ncol(y) != 1L) {
    message <- if (name == "y") "[y] must be a vector or one-column object." else
      "[DV.train] must contain exactly one response column."
    stop(message, call. = FALSE)
  }
  if (is.data.frame(y)) y[[1L]] else y[, 1L]
}

.nns_fold_columns <- function(fold, train, source) {
  # Two level labels can share a syntactic name, so match their original levels.
  by.predictor <- split(seq_along(source), source)
  columns <- lapply(names(train), function(nm) {
    idx <- by.predictor[[nm]]
    if (is.factor(fold[[nm]])) return(idx[match(levels(fold[[nm]]), levels(train[[nm]]))])
    if (is.factor(train[[nm]])) return(NA_integer_)
    idx
  })
  unlist(columns, use.names = FALSE)
}

.nns_reg_response <- function(y, type) {
  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || is.na(type)) {
      stop("[type] must be NULL, 'REG', 'CLASS', or 'XONLY'.", call. = FALSE)
    }
    type <- tolower(type)
    if (!type %in% c("reg", "class", "xonly")) {
      stop("[type] must be NULL, 'REG', 'CLASS', or 'XONLY'.", call. = FALSE)
    }
  }
  if (identical(type, "reg") && !is.numeric(y)) {
    stop("[type = 'REG'] requires a numeric response.", call. = FALSE)
  }

  # The automatic test hashes y, so it only runs when no type was given.
  is.class <- identical(type, "class") ||
    (is.null(type) &&
       (is.factor(y) || is.character(y) || is.logical(y) ||
          (is.numeric(y) && length(unique(y)) == 2L &&
             setequal(sort(unique(as.numeric(y))), c(0, 1)))))
  task <- .nns_encode_response(y, if (is.class) "class" else type)
  if (!is.null(task$class.levels)) task$class.levels <- as.character(task$class.levels)
  task
}

.nns_boost_response <- function(y, type) {
  if (identical(type, "class") && !is.numeric(y) && !is.factor(y) &&
      !is.logical(y) && !is.character(y)) y <- factor(y)
  .nns_encode_response(y, type)
}
