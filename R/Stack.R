#' NNS Stack
#'
#' Cross-validated ensemble of the full multivariate and synthetic-dimension
#' \link{NNS.reg} models.
#'
#' @param IVs.train a vector, matrix, or data frame of numeric, logical,
#'   character, factor, Date, or date-time predictors.
#' @param DV.train a numeric, logical, character, or factor response with one
#'   value per row of \code{IVs.train}.
#' @param IVs.test a vector, matrix, or data frame with the same predictors as
#'   \code{IVs.train}. If \code{NULL}, \code{IVs.train} is used.
#' @param type \code{NULL} (default) for regression or \code{"CLASS"} for
#'   classification. Factor, character, logical, and two-level numeric responses
#'   automatically select classification.
#' @param obj.fn an expression using \code{predicted} and \code{actual}.
#'   Sum of squared errors is the regression default. For classification, the
#'   untouched default is replaced by mean classification accuracy.
#' @param objective one of \code{"min"} or \code{"max"}.
#' @param optimize.threshold logical; optimize the class-rounding threshold from
#'   out-of-fold predictions. If \code{FALSE}, use 0.5.
#' @param dist distance option. The corrected implementation currently accepts
#'   only \code{"L2"}, because the production multivariate \code{NNS.reg}
#'   path does not presently implement distinct L1, DTW, or FACTOR estimators.
#' @param CV.size optional validation fraction in \code{(0, 1)}. If supplied,
#'   \code{folds} repeated stratified/random holdouts are used. If \code{NULL},
#'   disjoint k-fold cross-validation is used.
#' @param balance logical; balance only each fitting partition and the final
#'   fitting data. Validation observations are never resampled.
#' @param ts.test positive integer; validation-block length for chronological
#'   rolling-origin cross-validation. The final block always contains the most
#'   recent observations.
#' @param folds positive integer; number of ordinary or rolling-origin folds.
#' @param order integer, \code{"max"}, or \code{NULL}; passed unchanged to
#'   \link{NNS.reg}.
#' @param method any unique combination of \code{1} and \code{2}. Method 1 is
#'   the full multivariate \code{NNS.reg} model with cross-validated
#'   \code{n.best}; Method 2 is a synthetic X* dimension-reduction model.
#' @param stack logical; when both methods are requested, use fold-local and
#'   training-only X* as Method 1's input. If \code{FALSE}, Method 1 uses the
#'   full independently encoded predictor matrix.
#' @param dim.red.method one of \code{"cor"}, \code{"NNS.dep"},
#'   \code{"NNS.caus"}, \code{"equal"}, \code{"all"}, or a numeric
#'   coefficient vector aligned to the encoded design columns.
#' @param pred.int numeric in \code{(0, 1)} or \code{NULL}; prediction interval
#'   level for the final component fits.
#' @param status logical; print progress messages.
#' @param ncores positive integer or \code{NULL}; native thread count.
#' @param seed non-negative integer or \code{NULL}; local random seed. The
#'   caller's random-number state is restored on exit.
#'
#' @return A list retaining the historical fields:
#' \itemize{
#'   \item \code{OBJfn.reg}: selected Method 1 out-of-fold objective.
#'   \item \code{NNS.reg.n.best}: selected Method 1 \code{n.best}.
#'   \item \code{probability.threshold}: threshold optimized directly on the
#'     out-of-fold weighted ensemble, or 0.5 for regression.
#'   \item \code{OBJfn.dim.red}: selected Method 2 out-of-fold objective.
#'   \item \code{NNS.dim.red.threshold}: full-data coefficient-magnitude
#'     cutoff corresponding to the selected active-dimension count.
#'   \item \code{reg}, \code{dim.red}, and \code{stack}: final predictions.
#'   \item component and stacked prediction intervals.
#' }
#' Classification predictions are returned as numeric class codes, matching the
#' historical NNS.stack interface. Additional fields \code{weights} and
#' \code{class.levels} report the out-of-fold blend and the code-to-label map.
#'
#' @note
#' Categorical encoding and min-max normalization are fitted on each training
#' partition only and then applied unchanged to its validation partition. The
#' final transformations are fitted on the complete training data only; external
#' test observations never affect training scales.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#'   \doi{10.2139/ssrn.2864711}
#'
#' @examples
#' \dontrun{
#' fit <- NNS.stack(
#'   iris[1:140, 1:4], iris[1:140, 5],
#'   IVs.test = iris[141:150, 1:4],
#'   type = "CLASS", balance = TRUE
#' )
#' fit$stack
#' }
#'
#' @export


NNS.stack <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      type = NULL,
                      obj.fn = expression(sum((predicted - actual)^2)),
                      objective = "min",
                      optimize.threshold = TRUE,
                      dist = "L2",
                      CV.size = NULL,
                      balance = FALSE,
                      ts.test = NULL,
                      folds = 5,
                      order = NULL,
                      method = c(1, 2),
                      stack = TRUE,
                      dim.red.method = "cor",
                      pred.int = NULL,
                      status = TRUE,
                      ncores = NULL,
                      seed = 123L) {
  
  # -------------------------------------------------------------------------
  # Validation, coercion, and state helpers
  # -------------------------------------------------------------------------
  
  .scalar_logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("[%s] must be TRUE or FALSE.", name), call. = FALSE)
    }
    x
  }
  
  .scalar_integer <- function(x, name, minimum = 0L, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x)) {
      stop(sprintf("[%s] must be an integer >= %d.", name, minimum),
           call. = FALSE)
    }
    as.integer(x)
  }
  
  .as_train_frame <- function(x) {
    if (any(class(x) %in% c("tbl", "data.table"))) x <- as.data.frame(x)
    if (is.null(dim(x))) x <- data.frame(X1 = x, check.names = FALSE)
    x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) < 1L) {
      stop("[IVs.train] must contain at least one predictor.", call. = FALSE)
    }
    if (is.null(names(x)) || any(names(x) == "")) {
      names(x) <- paste0("X", seq_len(ncol(x)))
    }
    if (anyDuplicated(names(x))) {
      stop("[IVs.train] predictor names must be unique.", call. = FALSE)
    }
    x
  }
  
  .as_test_frame <- function(x, train_names) {
    p <- length(train_names)
    had_names <- !is.null(dim(x)) && !is.null(colnames(x)) &&
      length(colnames(x)) == NCOL(x) && all(colnames(x) != "")
    
    if (any(class(x) %in% c("tbl", "data.table"))) {
      had_names <- !is.null(names(x)) && all(names(x) != "")
      x <- as.data.frame(x)
    }
    
    if (is.null(dim(x))) {
      if (p == 1L) {
        x <- data.frame(x, check.names = FALSE)
        names(x) <- train_names
      } else if (length(x) == p) {
        supplied <- names(x)
        if (!is.null(supplied) && all(nzchar(supplied)) &&
            !identical(make.unique(supplied, sep = "."), train_names)) {
          if (anyDuplicated(supplied) || !setequal(supplied, train_names)) {
            stop("Named [IVs.test] values must exactly match the training predictors.",
                 call. = FALSE)
          }
          x <- x[train_names]
        }
        x <- as.data.frame(as.list(x), check.names = FALSE,
                           stringsAsFactors = FALSE)
        names(x) <- train_names
      } else {
        stop(paste0(
          "A vector [IVs.test] must contain one complete test row, unless ",
          "[IVs.train] has one predictor."
        ), call. = FALSE)
      }
    } else {
      x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
    if (ncol(x) != p) {
      stop("[IVs.test] must have the same number of predictors as [IVs.train].",
           call. = FALSE)
    }
    
    if (!had_names) {
      names(x) <- train_names
    } else {
      if (anyDuplicated(names(x))) {
        stop("[IVs.test] predictor names must be unique.", call. = FALSE)
      }
      missing_names <- setdiff(train_names, names(x))
      extra_names <- setdiff(names(x), train_names)
      if (length(missing_names) || length(extra_names)) {
        stop(sprintf(
          paste0("[IVs.test] columns must exactly match [IVs.train]. ",
                 "Missing: %s; extra: %s."),
          if (length(missing_names)) paste(missing_names, collapse = ", ") else "none",
          if (length(extra_names)) paste(extra_names, collapse = ", ") else "none"
        ), call. = FALSE)
      }
      x <- x[, train_names, drop = FALSE]
    }
    
    x
  }
  
  .align_predictors <- function(train, test) {
    for (j in seq_along(train)) {
      nm <- names(train)[j]
      tr <- train[[j]]
      te <- test[[j]]
      
      if (inherits(tr, "Date")) {
        if (!inherits(te, "Date")) {
          stop(sprintf("Test predictor [%s] must also be a Date.", nm),
               call. = FALSE)
        }
        train[[j]] <- as.numeric(tr)
        test[[j]] <- as.numeric(te)
      } else if (inherits(tr, c("POSIXct", "POSIXlt"))) {
        if (!inherits(te, c("POSIXct", "POSIXlt"))) {
          stop(sprintf("Test predictor [%s] must also be a date-time value.", nm),
               call. = FALSE)
        }
        train[[j]] <- as.numeric(tr)
        test[[j]] <- as.numeric(te)
      } else if (is.factor(tr) || is.character(tr)) {
        tr_chr <- as.character(tr)
        te_chr <- as.character(te)
        lev <- if (is.factor(tr)) levels(tr) else unique(tr_chr)
        unseen <- setdiff(unique(te_chr), lev)
        if (length(unseen)) {
          stop(sprintf(
            "Test predictor [%s] contains unseen level(s): %s.",
            nm, paste(unseen, collapse = ", ")
          ), call. = FALSE)
        }
        train[[j]] <- factor(tr_chr, levels = lev, ordered = is.ordered(tr))
        test[[j]] <- factor(te_chr, levels = lev, ordered = is.ordered(tr))
      } else if (is.logical(tr)) {
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
    list(train = train, test = test)
  }
  
  .check_predictors <- function(x, name) {
    if (anyNA(x)) stop(sprintf("[%s] contains missing values.", name), call. = FALSE)
    for (j in seq_along(x)) {
      if (is.numeric(x[[j]]) && any(!is.finite(x[[j]]))) {
        stop(sprintf("[%s] predictor [%s] contains non-finite values.",
                     name, names(x)[j]), call. = FALSE)
      }
    }
  }
  
  .restore_rng <- local({
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old <- if (existed) get(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE) else NULL
    function() {
      if (existed) {
        assign(".Random.seed", old, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  })
  on.exit(.restore_rng(), add = TRUE)
  
  optimize.threshold <- .scalar_logical(optimize.threshold, "optimize.threshold")
  balance <- .scalar_logical(balance, "balance")
  stack <- .scalar_logical(stack, "stack")
  status <- .scalar_logical(status, "status")
  folds <- .scalar_integer(folds, "folds", minimum = 1L)
  ts.test <- .scalar_integer(ts.test, "ts.test", minimum = 1L,
                             allow_null = TRUE)
  
  if (!is.null(seed)) {
    seed <- .scalar_integer(seed, "seed", minimum = 0L)
    set.seed(seed)
  }
  
  if (is.null(obj.fn) || !(is.expression(obj.fn) || is.call(obj.fn))) {
    stop("[obj.fn] must be a non-NULL expression or call.", call. = FALSE)
  }
  if (is.expression(obj.fn) && length(obj.fn) != 1L) {
    stop("[obj.fn] must contain exactly one expression.", call. = FALSE)
  }
  
  if (!is.character(objective) || length(objective) != 1L || is.na(objective)) {
    stop("[objective] must be exactly 'min' or 'max'.", call. = FALSE)
  }
  objective <- match.arg(tolower(objective), c("min", "max"))
  
  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || is.na(type)) {
      stop("[type] must be NULL or 'CLASS'.", call. = FALSE)
    }
    type <- match.arg(tolower(type), "class")
  }
  
  if (!is.null(order)) {
    if (is.character(order)) {
      if (length(order) != 1L || tolower(order) != "max") {
        stop("[order] must be NULL, a positive integer, or 'max'.",
             call. = FALSE)
      }
      order <- "max"
    } else {
      order <- .scalar_integer(order, "order", minimum = 1L)
    }
  }
  
  if (!is.numeric(method) || !length(method) || any(!is.finite(method)) ||
      any(method != floor(method))) {
    stop("[method] must contain only 1 and/or 2.", call. = FALSE)
  }
  method <- sort(unique(as.integer(method)))
  if (!all(method %in% c(1L, 2L))) {
    stop("[method] must contain only 1 and/or 2.", call. = FALSE)
  }
  
  if (!is.null(CV.size)) {
    if (!is.numeric(CV.size) || length(CV.size) != 1L ||
        !is.finite(CV.size) || CV.size <= 0 || CV.size >= 1) {
      stop("[CV.size] must be a finite scalar strictly between 0 and 1.",
           call. = FALSE)
    }
    CV.size <- as.numeric(CV.size)
  }
  
  if (!is.null(pred.int)) {
    if (!is.numeric(pred.int) || length(pred.int) != 1L ||
        !is.finite(pred.int) || pred.int <= 0 || pred.int >= 1) {
      stop("[pred.int] must be a finite scalar strictly between 0 and 1.",
           call. = FALSE)
    }
    pred.int <- as.numeric(pred.int)
  }
  
  if (is.null(ncores)) {
    detected <- suppressWarnings(parallel::detectCores())
    if (!is.finite(detected)) detected <- 2L
    ncores <- max(1L, as.integer(detected) - 1L)
  } else {
    ncores <- .scalar_integer(ncores, "ncores", minimum = 1L)
  }
  
  if (!is.character(dist) || length(dist) != 1L || is.na(dist)) {
    stop("[dist] must be one character value.", call. = FALSE)
  }
  dist <- match.arg(tolower(dist), c("l2", "l1", "dtw", "factor"))
  if (dist != "l2") {
    stop(paste0(
      "The corrected NNS.stack currently supports dist = 'L2' only. ",
      "The production multivariate NNS.reg path does not yet implement ",
      "distinct L1, DTW, or FACTOR estimators, so those values are rejected ",
      "rather than silently treated as L2."
    ), call. = FALSE)
  }
  dist <- "L2"
  
  # -------------------------------------------------------------------------
  # Input data and response coding
  # -------------------------------------------------------------------------
  
  x <- .as_train_frame(IVs.train)
  
  if (any(class(DV.train) %in% c("tbl", "data.table"))) {
    DV.train <- as.vector(unlist(DV.train))
  }
  if (is.data.frame(DV.train) || is.matrix(DV.train)) {
    if (NCOL(DV.train) != 1L) {
      stop("[DV.train] must contain exactly one response column.",
           call. = FALSE)
    }
    DV.train <- as.vector(unlist(DV.train))
  }
  
  if (length(DV.train) != nrow(x)) {
    stop("nrow(IVs.train) must equal length(DV.train).", call. = FALSE)
  }
  if (length(DV.train) < 4L) {
    stop("NNS.stack requires at least four training observations.",
         call. = FALSE)
  }
  if (anyNA(DV.train)) {
    stop("[DV.train] contains missing values.", call. = FALSE)
  }
  if ((is.numeric(DV.train) || is.integer(DV.train)) &&
      any(!is.finite(DV.train))) {
    stop("[DV.train] contains non-finite values.", call. = FALSE)
  }
  
  response_was_factor <- is.factor(DV.train)
  response_was_ordered <- is.ordered(DV.train)
  response_was_character <- is.character(DV.train)
  response_was_logical <- is.logical(DV.train)
  response_was_numeric <- is.numeric(DV.train) || is.integer(DV.train)
  
  auto_class <- response_was_factor || response_was_character ||
    response_was_logical ||
    (response_was_numeric && length(unique(DV.train)) == 2L)
  
  if (auto_class && is.null(type)) type <- "class"
  if (balance && is.null(type)) {
    warning("type = 'CLASS' selected because balance = TRUE.", call. = FALSE)
    type <- "class"
  }
  is_class <- identical(type, "class")
  
  original_response <- DV.train
  class_values <- NULL
  
  if (is_class) {
    if (response_was_factor) {
      response_factor <- droplevels(DV.train)
      class_values <- levels(response_factor)
      y <- as.integer(response_factor)
    } else if (response_was_character) {
      class_values <- unique(DV.train)
      y <- match(DV.train, class_values)
    } else if (response_was_logical) {
      class_values <- sort(unique(DV.train))
      y <- match(DV.train, class_values)
    } else if (response_was_numeric) {
      class_values <- sort(unique(as.numeric(DV.train)))
      y <- match(as.numeric(DV.train), class_values)
    } else {
      stop("Unsupported classification response type.", call. = FALSE)
    }
    
    y <- as.numeric(y)
    if (length(unique(y)) < 2L) {
      stop("Classification requires at least two response classes.",
           call. = FALSE)
    }
    if (min(table(y)) < 2L) {
      stop("Each response class requires at least two observations for cross-validation.",
           call. = FALSE)
    }
    
    if (identical(obj.fn, expression(sum((predicted - actual)^2)))) {
      obj.fn <- expression(mean(predicted == actual))
      objective <- "max"
    }
  } else {
    if (!response_was_numeric) {
      stop("A nonnumeric response requires type = 'CLASS'.", call. = FALSE)
    }
    y <- as.numeric(DV.train)
  }
  
  if (balance && !is_class) {
    stop("[balance = TRUE] requires classification.", call. = FALSE)
  }
  
  if (is.null(IVs.test)) {
    z <- x
  } else {
    z <- .as_test_frame(IVs.test, names(x))
  }
  
  aligned <- .align_predictors(x, z)
  x <- aligned$train
  z <- aligned$test
  .check_predictors(x, "IVs.train")
  .check_predictors(z, "IVs.test")
  
  n_obs <- nrow(x)
  original_p <- ncol(x)
  
  if (!is.null(ts.test) && ts.test >= n_obs) {
    stop("[ts.test] must be smaller than the number of training observations.",
         call. = FALSE)
  }
  
  if (original_p == 1L && 2L %in% method) {
    warning("Method 2 was removed because dimension reduction requires more than one original predictor.",
            call. = FALSE)
    method <- setdiff(method, 2L)
    if (!length(method)) method <- 1L
  }
  
  if (2L %in% method) {
    if (is.numeric(dim.red.method)) {
      if (!length(dim.red.method) || any(!is.finite(dim.red.method))) {
        stop("Numeric [dim.red.method] coefficients must be finite.",
             call. = FALSE)
      }
    } else {
      if (!is.character(dim.red.method) || length(dim.red.method) != 1L ||
          is.na(dim.red.method)) {
        stop("[dim.red.method] must be one supported character value or a numeric vector.",
             call. = FALSE)
      }
      dim.red.method <- match.arg(
        tolower(dim.red.method),
        c("cor", "nns.dep", "nns.caus", "equal", "all")
      )
    }
  }
  
  # Add a constant non-integer translation to every fitting response. This keeps
  # NNS.reg from silently auto-switching integer-valued regression or class-code
  # targets into its internal classification path. Predictions and intervals are
  # translated back before scoring or returning.
  response_offset <- 0.123456789
  y_fit_all <- y + response_offset
  
  # -------------------------------------------------------------------------
  # Shared helpers
  # -------------------------------------------------------------------------
  
  .score <- function(predicted, actual) {
    if (length(predicted) != length(actual)) {
      stop("The objective received predicted and actual vectors of different lengths.",
           call. = FALSE)
    }
    value <- eval(obj.fn,
                  envir = list(predicted = predicted, actual = actual),
                  enclos = parent.frame())
    if (!is.numeric(value) || length(value) != 1L) {
      stop("[obj.fn] must return one numeric scalar.", call. = FALSE)
    }
    value <- as.numeric(value)
    if (!is.finite(value)) return(NA_real_)
    value
  }
  
  .sanitize_raw <- function(predicted, fallback_y) {
    predicted <- as.numeric(predicted)
    bad <- !is.finite(predicted)
    if (any(bad)) {
      good <- predicted[!bad]
      replacement <- if (length(good)) gravity(good) else gravity(fallback_y)
      if (!is.finite(replacement)) replacement <- mean(fallback_y)
      if (!is.finite(replacement)) {
        stop("NNS.reg returned no finite predictions.", call. = FALSE)
      }
      predicted[bad] <- replacement
    }
    predicted
  }
  
  .round_codes <- function(raw, threshold) {
    raw <- pmin(pmax(as.numeric(raw), 1), length(class_values))
    lo <- floor(raw)
    hi <- ceiling(raw)
    frac <- raw - lo
    out <- ifelse(frac < threshold, lo, hi)
    as.integer(pmin(pmax(out, 1L), length(class_values)))
  }
  
  .best_threshold <- function(raw, actual) {
    if (!is_class || !optimize.threshold) return(0.5)
    grid <- seq(0.01, 0.99, by = 0.01)
    scores <- vapply(grid, function(th) {
      .score(.round_codes(raw, th), actual)
    }, numeric(1L))
    valid <- which(is.finite(scores))
    if (!length(valid)) return(0.5)
    best_value <- if (objective == "min") {
      min(scores[valid])
    } else {
      max(scores[valid])
    }
    tied <- valid[scores[valid] == best_value]
    grid[tied[ceiling(length(tied) / 2)]]
  }
  
  .evaluate_raw <- function(raw, actual) {
    threshold <- if (is_class) .best_threshold(raw, actual) else 0.5
    predicted <- if (is_class) .round_codes(raw, threshold) else raw
    list(score = .score(predicted, actual), threshold = threshold)
  }
  
  .decode_codes <- function(code) {
    as.numeric(as.integer(pmin(pmax(code, 1L), length(class_values))))
  }
  
  .decode_interval <- function(interval, threshold) {
    if (is.null(interval)) return(NULL)
    interval <- as.data.frame(interval, check.names = FALSE)
    interval[] <- lapply(interval, function(v) {
      .decode_codes(.round_codes(v, threshold))
    })
    interval
  }
  
  .translate_interval <- function(interval) {
    if (is.null(interval)) return(NULL)
    interval <- as.data.frame(interval, check.names = FALSE)
    interval[] <- lapply(interval, function(v) as.numeric(v) - response_offset)
    interval
  }
  
  .has_all_classes <- function(train_y) {
    !is_class || identical(sort(unique(train_y)), sort(unique(y)))
  }
  
  .balance_indices <- function(train_y) {
    if (!balance) return(seq_along(train_y))
    groups <- split(seq_along(train_y), train_y)
    if (length(groups) < 2L || any(lengths(groups) == 0L)) {
      stop("Balancing requires at least two non-empty classes in the fitting split.",
           call. = FALSE)
    }
    smallest <- min(lengths(groups))
    largest <- max(lengths(groups))
    down_idx <- unlist(lapply(groups, function(g) {
      sample(g, smallest, replace = FALSE)
    }), use.names = FALSE)
    up_idx <- unlist(lapply(groups, function(g) {
      sample(g, largest, replace = TRUE)
    }), use.names = FALSE)
    sample(c(down_idx, up_idx), replace = FALSE)
  }
  
  .numeric_design <- function(train, test) {
    train_blocks <- vector("list", ncol(train))
    test_blocks <- vector("list", ncol(train))
    
    for (j in seq_along(train)) {
      nm <- names(train)[j]
      tr <- train[[j]]
      te <- test[[j]]
      
      if (is.factor(tr)) {
        lev <- levels(tr)
        tr_chr <- as.character(tr)
        te_chr <- as.character(te)
        tr_block <- vapply(lev, function(level) {
          as.numeric(tr_chr == level)
        }, numeric(length(tr_chr)))
        te_block <- vapply(lev, function(level) {
          as.numeric(te_chr == level)
        }, numeric(length(te_chr)))
        if (is.null(dim(tr_block))) tr_block <- matrix(tr_block, ncol = 1L)
        if (is.null(dim(te_block))) te_block <- matrix(te_block, ncol = 1L)
        block_names <- paste0(nm, "__", make.names(lev, unique = TRUE))
        colnames(tr_block) <- block_names
        colnames(te_block) <- block_names
      } else {
        tr_block <- matrix(as.numeric(tr), ncol = 1L,
                           dimnames = list(NULL, nm))
        te_block <- matrix(as.numeric(te), ncol = 1L,
                           dimnames = list(NULL, nm))
      }
      
      train_blocks[[j]] <- tr_block
      test_blocks[[j]] <- te_block
    }
    
    train_matrix <- do.call(cbind, train_blocks)
    test_matrix <- do.call(cbind, test_blocks)
    storage.mode(train_matrix) <- "double"
    storage.mode(test_matrix) <- "double"
    colnames(train_matrix) <- make.unique(colnames(train_matrix), sep = "_")
    colnames(test_matrix) <- colnames(train_matrix)
    
    train_min <- apply(train_matrix, 2L, min)
    train_max <- apply(train_matrix, 2L, max)
    train_range <- train_max - train_min
    train_range[!is.finite(train_range) | train_range == 0] <- 1
    
    train_scaled <- sweep(train_matrix, 2L, train_min, "-")
    train_scaled <- sweep(train_scaled, 2L, train_range, "/")
    test_scaled <- sweep(test_matrix, 2L, train_min, "-")
    test_scaled <- sweep(test_scaled, 2L, train_range, "/")
    
    train_scaled <- as.data.frame(train_scaled, check.names = FALSE)
    test_scaled <- as.data.frame(test_scaled, check.names = FALSE)
    
    list(train = train_scaled, test = test_scaled,
         minimum = train_min, range = train_range)
  }
  
  .make_splits <- function() {
    all_index <- seq_len(n_obs)
    
    if (!is.null(ts.test)) {
      possible <- floor((n_obs - 1L) / ts.test)
      if (possible < 1L) {
        stop("Not enough observations for the requested [ts.test].",
             call. = FALSE)
      }
      use_folds <- min(folds, possible)
      if (use_folds < folds) {
        warning(sprintf(
          "Only %d non-overlapping chronological fold(s) are available.",
          use_folds
        ), call. = FALSE)
      }
      starts <- n_obs - (use_folds:1L) * ts.test + 1L
      out <- lapply(starts, function(start) {
        validation <- seq.int(start, start + ts.test - 1L)
        training <- seq_len(start - 1L)
        list(train = training, validation = validation)
      })
      keep <- vapply(out, function(s) {
        length(s$train) >= 3L && .has_all_classes(y[s$train])
      }, logical(1L))
      out <- out[keep]
      if (!length(out)) {
        stop(paste0(
          "No chronological fold retained enough training observations and ",
          "all response classes."
        ), call. = FALSE)
      }
      return(out)
    }
    
    # Historical compatibility: folds = 1 means one repeated holdout even when
    # CV.size is omitted. The original implementation drew a validation
    # fraction between 0.20 and 1/3, so retain that behavior under the local
    # seed while still keeping validation observations completely unbalanced.
    holdout_size <- CV.size
    if (is.null(holdout_size) && folds == 1L) {
      holdout_size <- round(stats::runif(1L, 0.20, 1 / 3), 3L)
    }
    
    if (!is.null(holdout_size)) {
      out <- vector("list", folds)
      for (b in seq_len(folds)) {
        if (is_class) {
          groups <- split(all_index, y)
          validation <- unlist(lapply(groups, function(g) {
            size <- min(length(g) - 1L,
                        max(1L, as.integer(round(holdout_size * length(g)))))
            if (size <= 0L) integer() else sample(g, size, replace = FALSE)
          }), use.names = FALSE)
          validation <- sort(unique(validation))
        } else {
          size <- max(1L, min(n_obs - 1L,
                              as.integer(round(holdout_size * n_obs))))
          validation <- sort(sample.int(n_obs, size, replace = FALSE))
        }
        training <- setdiff(all_index, validation)
        if (length(training) < 3L || !.has_all_classes(y[training])) {
          stop(paste0(
            "Unable to create a repeated holdout retaining enough fitting ",
            "observations and every response class."
          ), call. = FALSE)
        }
        out[[b]] <- list(train = training, validation = validation)
      }
      return(out)
    }
    
    use_folds <- min(folds, n_obs)
    if (is_class) {
      class_counts <- table(y)
      if (min(class_counts) < 2L) {
        stop("Each response class requires at least two observations for cross-validation.",
             call. = FALSE)
      }
      use_folds <- min(use_folds, max(class_counts))
    }
    if (use_folds < 2L) {
      stop("At least two folds are required when [CV.size] and [ts.test] are NULL.",
           call. = FALSE)
    }
    if (use_folds < folds) {
      warning(sprintf("Cross-validation folds reduced to %d for the available data.",
                      use_folds), call. = FALSE)
    }
    
    fold_id <- integer(n_obs)
    if (is_class) {
      groups <- split(all_index, y)
      for (g in groups) {
        shuffled <- sample(g, length(g), replace = FALSE)
        fold_id[shuffled] <- rep(seq_len(use_folds), length.out = length(g))
      }
    } else {
      shuffled <- sample(all_index, n_obs, replace = FALSE)
      fold_id[shuffled] <- rep(seq_len(use_folds), length.out = n_obs)
    }
    
    out <- lapply(seq_len(use_folds), function(b) {
      validation <- which(fold_id == b)
      training <- setdiff(all_index, validation)
      list(train = training, validation = validation)
    })
    keep <- vapply(out, function(s) {
      length(s$validation) > 0L && length(s$train) >= 3L &&
        .has_all_classes(y[s$train])
    }, logical(1L))
    out <- out[keep]
    if (length(out) < 2L) {
      stop("Unable to create at least two valid cross-validation folds.",
           call. = FALSE)
    }
    out
  }
  
  splits <- .make_splits()
  
  .coefficient_vector <- function(X, response) {
    X <- as.matrix(X)
    p <- ncol(X)
    
    if (is.numeric(dim.red.method)) {
      coef <- as.numeric(dim.red.method)
      if (!is.null(names(dim.red.method))) {
        if (!setequal(names(dim.red.method), colnames(X))) {
          stop(paste0(
            "Named numeric [dim.red.method] coefficients must exactly match ",
            "the encoded design columns."
          ), call. = FALSE)
        }
        coef <- as.numeric(dim.red.method[colnames(X)])
      } else if (length(coef) != p) {
        stop(sprintf(
          "Numeric [dim.red.method] must contain %d encoded coefficients.", p
        ), call. = FALSE)
      }
      coef[!is.finite(coef)] <- 0
      return(coef)
    }
    
    cor_coef <- function() {
      out <- vapply(seq_len(p), function(j) {
        suppressWarnings(stats::cor(X[, j], response, method = "spearman"))
      }, numeric(1L))
      out[!is.finite(out)] <- 0
      out
    }
    
    dep_coef <- function() {
      out <- vapply(seq_len(p), function(j) {
        tryCatch(
          as.numeric(NNS.dep(X[, j], response,
                             print.map = FALSE, asym = TRUE)$Dependence)[1L],
          error = function(e) 0
        )
      }, numeric(1L))
      out[!is.finite(out)] <- 0
      out
    }
    
    caus_coef <- function() {
      tau_value <- if (is.null(ts.test)) "cs" else "ts"
      out <- vapply(seq_len(p), function(j) {
        tryCatch(
          as.numeric(Uni.caus(response, X[, j],
                              tau = tau_value, plot = FALSE))[1L],
          error = function(e) 0
        )
      }, numeric(1L))
      out[!is.finite(out)] <- 0
      out
    }
    
    coef <- switch(
      dim.red.method,
      "cor" = cor_coef(),
      "nns.dep" = dep_coef(),
      "nns.caus" = caus_coef(),
      "equal" = rep(1, p),
      "all" = rowMeans(cbind(caus_coef(), dep_coef(), cor_coef(), rep(1, p))),
      stop("Unsupported [dim.red.method].", call. = FALSE)
    )
    
    coef[!is.finite(coef)] <- 0
    if (!any(abs(coef) > 0)) coef <- rep(1, p)
    coef
  }
  
  .active_coefficients <- function(coef, count) {
    coef <- as.numeric(coef)
    count <- max(1L, min(as.integer(count), length(coef)))
    ord <- order(abs(coef), decreasing = TRUE, na.last = NA, method = "radix")
    active <- ord[seq_len(min(count, length(ord)))]
    out <- rep(0, length(coef))
    out[active] <- coef[active]
    if (!any(abs(out) > 0)) {
      out[active] <- 1
    }
    out
  }
  
  .project_xstar <- function(train_design, test_design, coef) {
    denom <- sum(abs(coef) > 0)
    if (denom < 1L) stop("No active dimension-reduction coefficients.",
                         call. = FALSE)
    list(
      train = as.numeric(as.matrix(train_design) %*% coef / denom),
      test = as.numeric(as.matrix(test_design) %*% coef / denom)
    )
  }
  
  .xstar_path <- function(train_design, test_design, coef) {
    ord <- order(abs(coef), decreasing = TRUE, na.last = NA, method = "radix")
    if (isTRUE(getOption("NNS.native.stack", TRUE))) {
      tryCatch(NNS_xstar_path_cpp(as.matrix(train_design), as.matrix(test_design),
                                  as.numeric(coef), as.integer(ord), as.integer(ncores)),
               error = function(e) NULL)
    } else NULL
  }
  
  .fit_univariate_raw <- function(train_x, train_y_fit, test_x,
                                  confidence.interval = NULL,
                                  point.only = TRUE,
                                  allow_failure = FALSE) {
    fit <- tryCatch(
      suppressWarnings({
        if (point.only && isTRUE(getOption("NNS.native.univariate", TRUE))) {
          fast <- .nns_reg_univariate_fast(as.numeric(train_x), as.numeric(train_y_fit),
                                           as.numeric(test_x), order = order,
                                           noise.reduction = "off", is.class = FALSE)
          list(Point.est = fast$prediction)
        } else {
          NNS.reg(
            as.numeric(train_x),
            as.numeric(train_y_fit),
            point.est = as.numeric(test_x),
            plot = FALSE,
            residual.plot = FALSE,
            order = order,
            type = NULL,
            factor.2.dummy = FALSE,
            dist = dist,
            ncores = ncores,
            point.only = point.only,
            confidence.interval = confidence.interval
          )
        }
      }),
      error = function(e) {
        if (allow_failure) return(NULL)
        stop(e)
      }
    )
    if (is.null(fit)) {
      return(list(raw = rep(NA_real_, length(test_x)), fit = NULL))
    }
    raw <- .sanitize_raw(as.numeric(fit$Point.est) - response_offset,
                         train_y_fit - response_offset)
    list(raw = raw, fit = fit)
  }
  
  
  # Score every k with exactly the estimator the final NNS.reg fit uses:
  # NNS_mreg_predict_path_cpp implements the repaired multivariate prediction
  # rule (range-normalized metric, stable ties, ensemble weights) for
  # k = 1..kmax. The repaired NNS.M.reg no longer extrapolates outside the
  # training support, so no gradient extension is applied here either.
  .production_multivariate_path <- function(rpm, Xtest, train_design, kmax = NULL) {
    Xtest <- as.data.frame(Xtest, check.names = FALSE)
    feature_names <- setdiff(names(rpm), "y.hat")
    Xtest <- Xtest[, feature_names, drop = FALSE]
    train_design <- as.data.frame(train_design, check.names = FALSE)
    train_design <- train_design[, feature_names, drop = FALSE]
    if (is.null(kmax)) kmax <- nrow(rpm)
    kmax <- min(kmax, nrow(rpm))
    
    rpm_x <- as.matrix(rpm[, feature_names, drop = FALSE])
    storage.mode(rpm_x) <- "double"
    test_matrix <- as.matrix(Xtest)
    storage.mode(test_matrix) <- "double"
    minimums <- vapply(train_design, min, numeric(1L))
    maximums <- vapply(train_design, max, numeric(1L))
    dist_code <- match(dist, c("L2", "L1", "FACTOR")) - 1L
    
    path <- if (isTRUE(getOption("NNS.native.stack", TRUE))) {
      NNS_mreg_predict_path_v2_cpp(
        rpm_x, as.numeric(rpm$y.hat), test_matrix, as.integer(kmax),
        dist_code, as.numeric(minimums), as.numeric(maximums), FALSE,
        as.integer(ncores)
      )
    } else {
      NNS_mreg_predict_path_cpp(
        rpm_x, as.numeric(rpm$y.hat), test_matrix, as.integer(kmax),
        dist_code, as.numeric(minimums), as.numeric(maximums), FALSE
      )
    }
    if (nrow(path) != nrow(test_matrix)) {
      stop("The production distance path returned an invalid row count.",
           call. = FALSE)
    }
    path
  }
  
  .candidate_from_oof <- function(sum_matrix, count_matrix, candidate) {
    count <- count_matrix[, candidate]
    valid <- count > 0L
    raw <- rep(NA_real_, n_obs)
    raw[valid] <- sum_matrix[valid, candidate] / count[valid]
    valid <- valid & is.finite(raw)
    if (!any(valid)) {
      return(list(score = NA_real_, threshold = 0.5, raw = raw))
    }
    evaluation <- .evaluate_raw(raw[valid], y[valid])
    list(score = evaluation$score,
         threshold = evaluation$threshold,
         raw = raw)
  }
  
  .select_best <- function(scores) {
    valid <- which(is.finite(scores))
    if (!length(valid)) {
      stop("No candidate produced a finite out-of-fold objective.",
           call. = FALSE)
    }
    best_value <- if (objective == "min") {
      min(scores[valid])
    } else {
      max(scores[valid])
    }
    tied <- valid[scores[valid] == best_value]
    as.integer(tied[ceiling(length(tied) / 2)])
  }
  
  # Full-data encoded design establishes the fixed encoded column domain.
  full_design <- .numeric_design(x, z)
  encoded_p <- ncol(full_design$train)
  
  if (2L %in% method && is.numeric(dim.red.method) &&
      is.null(names(dim.red.method)) && length(dim.red.method) != encoded_p) {
    stop(sprintf(
      "Numeric [dim.red.method] must contain %d encoded coefficients.",
      encoded_p
    ), call. = FALSE)
  }
  
  # -------------------------------------------------------------------------
  # Method 2: select active dimension count from OOF predictions
  # -------------------------------------------------------------------------
  
  dim_best_count <- NA_integer_
  dim_best_score <- NA_real_
  dim_component_threshold <- 0.5
  dim_oof_raw <- rep(NA_real_, n_obs)
  dim_threshold_report <- NA_real_
  dim_full_coef <- NULL
  dim_full_xstar_train <- NULL
  dim_full_xstar_test <- NULL
  
  if (2L %in% method) {
    dim_sum <- matrix(0, nrow = n_obs, ncol = encoded_p)
    dim_count <- matrix(0L, nrow = n_obs, ncol = encoded_p)
    
    for (b in seq_along(splits)) {
      split <- splits[[b]]
      train_idx <- split$train
      valid_idx <- split$validation
      fold_design <- .numeric_design(
        x[train_idx, , drop = FALSE],
        x[valid_idx, , drop = FALSE]
      )
      coef_fold <- .coefficient_vector(fold_design$train, y[train_idx])
      balanced_idx <- .balance_indices(y[train_idx])
      
      if (status) message(sprintf("Method 2 fold %d/%d: generating %d cumulative projections", b, length(splits), encoded_p))
      native_path <- .xstar_path(fold_design$train, fold_design$test, coef_fold)
      representatives <- if (!is.null(native_path)) as.integer(native_path$representative) else seq_len(encoded_p)
      unique_candidates <- unique(representatives)
      if (status) message(sprintf("Method 2 fold %d/%d: evaluating %d unique candidates", b, length(splits), length(unique_candidates)))
      candidate_raw <- vector("list", encoded_p)
      for (m in unique_candidates) {
        if (!is.null(native_path)) {
          projected <- list(train = native_path$train[, m], test = native_path$test[, m])
        } else {
          active_coef <- .active_coefficients(coef_fold, m)
          projected <- .project_xstar(fold_design$train, fold_design$test, active_coef)
        }
        fit <- .fit_univariate_raw(projected$train[balanced_idx],
                                   y_fit_all[train_idx][balanced_idx],
                                   projected$test,
                                   point.only = TRUE,
                                   allow_failure = TRUE)
        candidate_raw[[m]] <- fit$raw
      }
      for (m in seq_len(encoded_p)) {
        raw <- candidate_raw[[representatives[m]]]
        good <- is.finite(raw)
        if (any(good)) {
          rows <- valid_idx[good]
          dim_sum[rows, m] <- dim_sum[rows, m] + raw[good]
          dim_count[rows, m] <- dim_count[rows, m] + 1L
        }
      }
      if (status) message(sprintf("Method 2 fold %d/%d complete", b, length(splits)))
    }
    
    dim_scores <- rep(NA_real_, encoded_p)
    dim_thresholds <- rep(0.5, encoded_p)
    dim_raw_candidates <- vector("list", encoded_p)
    
    for (m in seq_len(encoded_p)) {
      candidate <- .candidate_from_oof(dim_sum, dim_count, m)
      dim_scores[m] <- candidate$score
      dim_thresholds[m] <- candidate$threshold
      dim_raw_candidates[[m]] <- candidate$raw
      if (status) {
        message(sprintf(
          paste0("Current dimension count = %d | OOF eval(obj.fn) = %s | ",
                 "Iterations remaining = %d"),
          m,
          if (is.finite(dim_scores[m]))
            format(dim_scores[m], digits = 6) else "NA",
          encoded_p - m
        ))
      }
    }
    
    dim_best_count <- .select_best(dim_scores)
    dim_best_score <- dim_scores[dim_best_count]
    dim_component_threshold <- dim_thresholds[dim_best_count]
    dim_oof_raw <- dim_raw_candidates[[dim_best_count]]
    
    dim_full_coef_original <- .coefficient_vector(full_design$train, y)
    dim_full_coef <- .active_coefficients(dim_full_coef_original,
                                          dim_best_count)
    dim_full_projection <- .xstar_path(full_design$train, full_design$test,
                                       dim_full_coef)
    if (is.null(dim_full_projection)) {
      dim_full_projection <- .project_xstar(
        full_design$train,
        full_design$test,
        dim_full_coef
      )
      dim_full_xstar_train <- dim_full_projection$train
      dim_full_xstar_test <- dim_full_projection$test
    } else {
      dim_full_xstar_train <- dim_full_projection$train[, dim_best_count]
      dim_full_xstar_test <- dim_full_projection$test[, dim_best_count]
    }
    
    active_magnitudes <- abs(dim_full_coef[abs(dim_full_coef) > 0])
    dim_threshold_report <- if (dim_best_count >= length(dim_full_coef)) {
      0
    } else if (length(active_magnitudes)) {
      min(active_magnitudes)
    } else {
      0
    }
    if (!is.finite(dim_threshold_report)) dim_threshold_report <- 0
  }
  
  # -------------------------------------------------------------------------
  # Method 1: bounded small-k search + mandatory all-points candidate
  # -------------------------------------------------------------------------
  
  reg_best_k <- NA_integer_
  reg_best_score <- NA_real_
  reg_component_threshold <- 0.5
  reg_oof_raw <- rep(NA_real_, n_obs)
  
  .method1_design_for_split <- function(train_idx, valid_idx) {
    fold_design <- .numeric_design(
      x[train_idx, , drop = FALSE],
      x[valid_idx, , drop = FALSE]
    )
    
    if (stack && 2L %in% method) {
      coef_fold <- .coefficient_vector(fold_design$train, y[train_idx])
      active_coef <- .active_coefficients(coef_fold, dim_best_count)
      projected <- .project_xstar(
        fold_design$train,
        fold_design$test,
        active_coef
      )
      list(
        train = data.frame(Xstar = projected$train,
                           Xstar2 = projected$train,
                           check.names = FALSE),
        test = data.frame(Xstar = projected$test,
                          Xstar2 = projected$test,
                          check.names = FALSE)
      )
    } else {
      list(train = fold_design$train, test = fold_design$test)
    }
  }
  
  if (1L %in% method) {
    # Compute the small-k limit based on training observations count
    l <- max(1L, floor(sqrt(n_obs)))
    
    # Structures to aggregate OOF predictions per candidate ID
    # Candidate IDs: integer 1..l and the special "all"
    candidate_ids <- c(as.character(seq_len(l)), "all")
    # We'll store sums and counts in lists keyed by candidate ID
    sum_list <- setNames(vector("list", length(candidate_ids)), candidate_ids)
    count_list <- setNames(vector("list", length(candidate_ids)), candidate_ids)
    for (id in candidate_ids) {
      sum_list[[id]] <- rep(0, n_obs)
      count_list[[id]] <- rep(0L, n_obs)
    }
    
    for (b in seq_along(splits)) {
      split <- splits[[b]]
      train_idx <- split$train
      valid_idx <- split$validation
      if (status) message(sprintf("Method 1 fold %d/%d: preparing fold design", b, length(splits)))
      design <- .method1_design_for_split(train_idx, valid_idx)
      balanced_idx <- .balance_indices(y[train_idx])
      train_design <- design$train[balanced_idx, , drop = FALSE]
      train_y_fit <- y_fit_all[train_idx][balanced_idx]
      valid_design <- design$test
      
      # We'll collect predictions for small candidates and the all candidate
      # For univariate vs multivariate
      if (ncol(train_design) == 1L) {
        # Univariate case: compute distances and cumulative averages for k=1..l
        train_x <- as.numeric(train_design[[1L]])
        test_x <- as.numeric(valid_design[[1L]])
        # We'll compute predictions for k=1..small_kmax (where small_kmax = min(l, length(train_x)))
        small_kmax <- min(l, length(train_x))
        # Sort distances and get indices
        # For each test point, we need sorted distances; we'll do a loop or use R's order
        # Efficient: compute distance matrix? For large train, we can compute per test point.
        # Since this is univariate, we can simply loop over test points and sort distances.
        # We'll compute a matrix of predictions for k=1..small_kmax.
        pred_mat <- matrix(NA, nrow = length(test_x), ncol = small_kmax)
        for (i in seq_along(test_x)) {
          dists <- abs(train_x - test_x[i])
          ord <- order(dists, method = "radix")
          cumsum_y <- cumsum(train_y_fit[ord])
          pred_mat[i, ] <- cumsum_y[seq_len(small_kmax)] / seq_len(small_kmax)
        }
        # All candidate prediction: mean of train_y_fit (constant)
        all_pred <- mean(train_y_fit)  # in offset scale
        # Convert to raw (subtract offset) for scoring
        # We'll subtract offset when evaluating
        # For small candidates, we need to evaluate each k sequentially
        # We'll store raw predictions for each candidate
        small_raw <- vector("list", small_kmax)
        for (k in seq_len(small_kmax)) {
          raw <- pred_mat[, k] - response_offset
          small_raw[[k]] <- raw
        }
        # All candidate raw
        all_raw <- rep(all_pred - response_offset, length(test_x))
        # Now evaluate small candidates with early stop
        small_scores <- numeric(small_kmax)
        small_thresholds <- numeric(small_kmax)
        stopped_at <- small_kmax
        for (k in seq_len(small_kmax)) {
          raw <- small_raw[[k]]
          eval <- .evaluate_raw(raw, y[valid_idx])
          small_scores[k] <- eval$score
          small_thresholds[k] <- eval$threshold
          if (status && (k %% 10 == 0 || k == small_kmax)) {
            message(sprintf("  k = %d, score = %s", k, format(eval$score, digits = 6)))
          }
          if (k >= 4) {
            # early stop condition
            cond <- if (objective == "min") {
              small_scores[k] >= small_scores[k-1] && small_scores[k] >= small_scores[k-2]
            } else {
              small_scores[k] <= small_scores[k-1] && small_scores[k] <= small_scores[k-2]
            }
            if (cond) {
              stopped_at <- k
              if (status) message(sprintf("  early stopping at k = %d", k))
              break
            }
          }
        }
        # Evaluate all candidate
        all_eval <- .evaluate_raw(all_raw, y[valid_idx])
        all_score <- all_eval$score
        all_threshold <- all_eval$threshold
        
        # Aggregate small candidates that were evaluated (up to stopped_at)
        for (k in seq_len(stopped_at)) {
          id <- as.character(k)
          raw <- small_raw[[k]]
          good <- is.finite(raw)
          if (any(good)) {
            rows <- valid_idx[good]
            sum_list[[id]][rows] <- sum_list[[id]][rows] + raw[good]
            count_list[[id]][rows] <- count_list[[id]][rows] + 1L
          }
        }
        # Aggregate all candidate
        id <- "all"
        raw <- all_raw
        good <- is.finite(raw)
        if (any(good)) {
          rows <- valid_idx[good]
          sum_list[[id]][rows] <- sum_list[[id]][rows] + raw[good]
          count_list[[id]][rows] <- count_list[[id]][rows] + 1L
        }
        
      } else {
        # Multivariate case
        if (status) message(sprintf("Method 1 fold %d/%d: building partitions", b, length(splits)))
        setup <- .nns_mreg_prepare_model(train_design, train_y_fit, order = order,
                                         noise.reduction = "off", is.class = FALSE,
                                         use.native = isTRUE(getOption("NNS.native.mreg", TRUE)))
        rpm <- setup$RPM
        if (is.null(rpm) || !is.data.frame(rpm) || nrow(rpm) < 1L || !"y.hat" %in% names(rpm)) {
          stop("NNS.reg did not return a usable regression-point matrix.", call. = FALSE)
        }
        nRPM <- nrow(rpm)
        if (status) message(sprintf("RPM rows = %d; validation rows = %d", nRPM, length(valid_idx)))
        
        small_kmax <- min(l, nRPM)
        if (small_kmax >= 1) {
          if (status) message(sprintf("  small candidates = 1...%d", small_kmax))
          # Compute small path for k=1..small_kmax
          small_path <- .production_multivariate_path(rpm = rpm, Xtest = valid_design,
                                                      train_design = train_design,
                                                      kmax = small_kmax)
          # subtract offset
          small_path <- small_path - response_offset
          # Evaluate small candidates with early stop
          small_scores <- numeric(small_kmax)
          small_thresholds <- numeric(small_kmax)
          stopped_at <- small_kmax
          for (k in seq_len(small_kmax)) {
            raw <- small_path[, k]
            raw <- .sanitize_raw(raw, y[train_idx])
            eval <- .evaluate_raw(raw, y[valid_idx])
            small_scores[k] <- eval$score
            small_thresholds[k] <- eval$threshold
            if (status && (k %% 10 == 0 || k == small_kmax)) {
              message(sprintf("  k = %d, score = %s", k, format(eval$score, digits = 6)))
            }
            if (k >= 4) {
              cond <- if (objective == "min") {
                small_scores[k] >= small_scores[k-1] && small_scores[k] >= small_scores[k-2]
              } else {
                small_scores[k] <= small_scores[k-1] && small_scores[k] <= small_scores[k-2]
              }
              if (cond) {
                stopped_at <- k
                if (status) message(sprintf("  early stopping at k = %d", k))
                break
              }
            }
          }
          # Aggregate small candidates that were evaluated (up to stopped_at)
          for (k in seq_len(stopped_at)) {
            id <- as.character(k)
            raw <- small_path[, k]
            good <- is.finite(raw)
            if (any(good)) {
              rows <- valid_idx[good]
              sum_list[[id]][rows] <- sum_list[[id]][rows] + raw[good]
              count_list[[id]][rows] <- count_list[[id]][rows] + 1L
            }
          }
        } else {
          stopped_at <- 0
        }
        
        # All candidate: prediction = mean of RPM$y.hat (or mean of train_y_fit)
        # Use mean of train_y_fit (balanced) for consistency
        all_raw <- rep(mean(train_y_fit) - response_offset, length(valid_idx))
        # Evaluate all candidate
        all_eval <- .evaluate_raw(all_raw, y[valid_idx])
        all_score <- all_eval$score
        all_threshold <- all_eval$threshold
        if (status) message(sprintf("  limit candidate = all (k = %d), score = %s", nRPM, format(all_score, digits = 6)))
        
        # Aggregate all candidate
        id <- "all"
        raw <- all_raw
        good <- is.finite(raw)
        if (any(good)) {
          rows <- valid_idx[good]
          sum_list[[id]][rows] <- sum_list[[id]][rows] + raw[good]
          count_list[[id]][rows] <- count_list[[id]][rows] + 1L
        }
      } # end multivariate
    } # end folds
    
    # Now compute scores for all candidates (small and all) from aggregated OOF predictions
    candidate_scores <- setNames(rep(NA_real_, length(candidate_ids)), candidate_ids)
    candidate_thresholds <- setNames(rep(0.5, length(candidate_ids)), candidate_ids)
    candidate_raws <- setNames(vector("list", length(candidate_ids)), candidate_ids)
    
    for (id in candidate_ids) {
      sum_vec <- sum_list[[id]]
      count_vec <- count_list[[id]]
      # Aggregate raw predictions for this candidate
      raw <- rep(NA_real_, n_obs)
      valid <- count_vec > 0L
      raw[valid] <- sum_vec[valid] / count_vec[valid]
      candidate_raws[[id]] <- raw
      eval <- .evaluate_raw(raw[valid], y[valid])  # only evaluate where valid
      candidate_scores[id] <- eval$score
      candidate_thresholds[id] <- eval$threshold
    }
    
    # Select best among all candidates (including "all")
    # We need to find the best score; objective min or max
    valid_ids <- names(candidate_scores)[is.finite(candidate_scores)]
    if (length(valid_ids) == 0) {
      stop("No Method 1 candidate produced a finite OOF objective.", call. = FALSE)
    }
    best_val <- if (objective == "min") {
      min(candidate_scores[valid_ids])
    } else {
      max(candidate_scores[valid_ids])
    }
    best_ids <- valid_ids[candidate_scores[valid_ids] == best_val]
    # if tie, take the first (or middle) - we'll take first
    best_id <- best_ids[1]
    
    # Set reg_best_k and other outputs
    if (best_id == "all") {
      # For all candidate, we need to know the full-data RPM row count for final fitting
      # We'll compute later; for now set a placeholder, but we need to store that it's "all"
      reg_best_k <- NA_integer_  # we'll set after full-data fit
      reg_best_label <- "all"
    } else {
      reg_best_k <- as.integer(best_id)
      reg_best_label <- "integer"
    }
    reg_best_score <- as.numeric(candidate_scores[best_id])
    reg_component_threshold <- as.numeric(candidate_thresholds[best_id])
    reg_oof_raw <- candidate_raws[[best_id]]
    
    # For final fitting, we need to know the full-data RPM row count when "all" wins
    # We'll compute that later in the final fit section.
    # Store best_id for use later.
    .reg_best_id <- best_id
    .reg_best_label <- if (best_id == "all") "all" else "integer"
    
    if (status) {
      if (best_id == "all") {
        message(sprintf("Best Method 1 candidate: all (k = all RPM rows), score = %s", format(reg_best_score, digits = 6)))
      } else {
        message(sprintf("Best Method 1 candidate: k = %d, score = %s", reg_best_k, format(reg_best_score, digits = 6)))
      }
    }
  } else {
    .reg_best_id <- NA_character_
    .reg_best_label <- NA_character_
  }
  
  # -------------------------------------------------------------------------
  # Optimize the actual OOF blend and its final classification threshold
  # -------------------------------------------------------------------------
  
  component_weights <- c(reg = 0, dim.red = 0)
  probability.threshold <- if (is_class) 0.5 else 0.5
  
  if (identical(method, 1L)) {
    component_weights["reg"] <- 1
    probability.threshold <- reg_component_threshold
  } else if (identical(method, 2L)) {
    component_weights["dim.red"] <- 1
    probability.threshold <- dim_component_threshold
  } else {
    valid <- is.finite(reg_oof_raw) & is.finite(dim_oof_raw)
    if (!any(valid)) {
      stop("The two component models have no common finite OOF predictions.",
           call. = FALSE)
    }
    
    weight_grid <- seq(0, 1, by = 0.01)
    blend_scores <- rep(NA_real_, length(weight_grid))
    blend_thresholds <- rep(0.5, length(weight_grid))
    
    for (i in seq_along(weight_grid)) {
      w <- weight_grid[i]
      raw <- w * reg_oof_raw[valid] + (1 - w) * dim_oof_raw[valid]
      evaluation <- .evaluate_raw(raw, y[valid])
      blend_scores[i] <- evaluation$score
      blend_thresholds[i] <- evaluation$threshold
    }
    
    valid_grid <- which(is.finite(blend_scores))
    if (!length(valid_grid)) {
      stop("No ensemble weight produced a finite OOF objective.",
           call. = FALSE)
    }
    best_value <- if (objective == "min") {
      min(blend_scores[valid_grid])
    } else {
      max(blend_scores[valid_grid])
    }
    tied <- valid_grid[blend_scores[valid_grid] == best_value]
    selected_index <- as.integer(tied[ceiling(length(tied) / 2)])
    selected_weight <- weight_grid[selected_index]
    
    component_weights <- c(reg = selected_weight,
                           dim.red = 1 - selected_weight)
    probability.threshold <- blend_thresholds[selected_index]
  }
  
  # -------------------------------------------------------------------------
  # Final production fits on complete training data
  # -------------------------------------------------------------------------
  
  if (status) message("Generating final estimates")
  
  reg_raw_final <- NULL
  reg_pred_int_raw <- NULL
  dim_raw_final <- NULL
  dim_pred_int_raw <- NULL
  
  if (2L %in% method) {
    balanced_idx <- .balance_indices(y)
    dim_final_fit <- .fit_univariate_raw(
      dim_full_xstar_train[balanced_idx],
      y_fit_all[balanced_idx],
      dim_full_xstar_test,
      confidence.interval = pred.int,
      point.only = FALSE
    )
    dim_raw_final <- dim_final_fit$raw
    dim_pred_int_raw <- .translate_interval(dim_final_fit$fit$pred.int)
  }
  
  .method1_full_design <- function() {
    if (stack && 2L %in% method) {
      list(
        train = data.frame(Xstar = dim_full_xstar_train,
                           Xstar2 = dim_full_xstar_train,
                           check.names = FALSE),
        test = data.frame(Xstar = dim_full_xstar_test,
                          Xstar2 = dim_full_xstar_test,
                          check.names = FALSE)
      )
    } else {
      list(train = full_design$train, test = full_design$test)
    }
  }
  
  if (1L %in% method) {
    final_design <- .method1_full_design()
    balanced_idx <- .balance_indices(y)
    train_design <- final_design$train[balanced_idx, , drop = FALSE]
    train_y_fit <- y_fit_all[balanced_idx]
    test_design <- final_design$test
    
    # Determine n.best for final fit
    if (.reg_best_label == "all") {
      # For all candidate, we need n.best = "all" (or numeric nrow(RPM))
      # We'll use n.best = "all" in NNS.reg
      n.best <- "all"
    } else {
      n.best <- reg_best_k
    }
    
    if (ncol(train_design) == 1L) {
      # Univariate final fit
      # Use NNS.reg with n.best = n.best
      if (is.character(n.best) && n.best == "all") {
        # For univariate, "all" means use all training points; we can pass n.best = "all"
        # but NNS.reg univariate might not accept "all"? We'll use NNS.reg with n.best = "all"
        reg_final_fit <- .fit_univariate_raw(
          train_design[[1L]],
          train_y_fit,
          test_design[[1L]],
          confidence.interval = pred.int,
          point.only = FALSE
        )
        # But that uses default n.best? Actually .fit_univariate_raw does not pass n.best.
        # We need to handle separately.
        # For univariate, we can just use NNS.reg with n.best = "all" explicitly.
        reg_final_fit <- suppressWarnings(
          NNS.reg(
            as.numeric(train_design[[1L]]),
            as.numeric(train_y_fit),
            point.est = as.numeric(test_design[[1L]]),
            plot = FALSE,
            residual.plot = FALSE,
            order = order,
            n.best = "all",
            type = NULL,
            factor.2.dummy = FALSE,
            dist = dist,
            ncores = ncores,
            point.only = FALSE,
            confidence.interval = pred.int
          )
        )
        reg_raw_final <- .sanitize_raw(
          as.numeric(reg_final_fit$Point.est) - response_offset,
          y
        )
        reg_pred_int_raw <- .translate_interval(reg_final_fit$pred.int)
      } else {
        # numeric n.best
        reg_final_fit <- suppressWarnings(
          NNS.reg(
            as.numeric(train_design[[1L]]),
            as.numeric(train_y_fit),
            point.est = as.numeric(test_design[[1L]]),
            plot = FALSE,
            residual.plot = FALSE,
            order = order,
            n.best = n.best,
            type = NULL,
            factor.2.dummy = FALSE,
            dist = dist,
            ncores = ncores,
            point.only = FALSE,
            confidence.interval = pred.int
          )
        )
        reg_raw_final <- .sanitize_raw(
          as.numeric(reg_final_fit$Point.est) - response_offset,
          y
        )
        reg_pred_int_raw <- .translate_interval(reg_final_fit$pred.int)
      }
    } else {
      # Multivariate final fit
      if (is.character(n.best) && n.best == "all") {
        # Use NNS.reg with n.best = "all"
        reg_final_fit <- suppressWarnings(
          NNS.reg(
            train_design,
            train_y_fit,
            point.est = test_design,
            plot = FALSE,
            residual.plot = FALSE,
            n.best = "all",
            order = order,
            type = NULL,
            factor.2.dummy = FALSE,
            dist = dist,
            ncores = ncores,
            point.only = FALSE,
            confidence.interval = pred.int
          )
        )
        reg_raw_final <- .sanitize_raw(
          as.numeric(reg_final_fit$Point.est) - response_offset,
          y
        )
        reg_pred_int_raw <- .translate_interval(reg_final_fit$pred.int)
        # Also set reg_best_k to the full RPM row count for return field
        # We need to get the RPM row count from the final fit? We can compute from the fit object.
        # But we can also compute by calling .nns_mreg_prepare_model on full data.
        # We'll compute later if needed.
        # For now, we'll set a placeholder.
        # Actually we'll compute after the fit.
        # We'll store the nRPM in a variable.
        # Let's compute RPM row count for full data:
        setup_full <- .nns_mreg_prepare_model(train_design, train_y_fit, order = order,
                                              noise.reduction = "off", is.class = FALSE,
                                              use.native = isTRUE(getOption("NNS.native.mreg", TRUE)))
        if (!is.null(setup_full$RPM)) {
          reg_best_k <- nrow(setup_full$RPM)
        } else {
          reg_best_k <- nrow(train_design)  # fallback
        }
      } else {
        # numeric n.best
        reg_final_fit <- suppressWarnings(
          NNS.reg(
            train_design,
            train_y_fit,
            point.est = test_design,
            plot = FALSE,
            residual.plot = FALSE,
            n.best = n.best,
            order = order,
            type = NULL,
            factor.2.dummy = FALSE,
            dist = dist,
            ncores = ncores,
            point.only = FALSE,
            confidence.interval = pred.int
          )
        )
        reg_raw_final <- .sanitize_raw(
          as.numeric(reg_final_fit$Point.est) - response_offset,
          y
        )
        reg_pred_int_raw <- .translate_interval(reg_final_fit$pred.int)
      }
    }
  }
  
  # Component predictions retain their own OOF-optimized thresholds. The stack
  # uses the threshold optimized directly on the OOF weighted ensemble.
  if (is_class) {
    reg_output <- if (!is.null(reg_raw_final)) {
      .decode_codes(.round_codes(reg_raw_final, reg_component_threshold))
    } else {
      NA
    }
    dim_output <- if (!is.null(dim_raw_final)) {
      .decode_codes(.round_codes(dim_raw_final, dim_component_threshold))
    } else {
      NA
    }
  } else {
    reg_output <- if (!is.null(reg_raw_final)) reg_raw_final else NA_real_
    dim_output <- if (!is.null(dim_raw_final)) dim_raw_final else NA_real_
  }
  
  if (identical(method, 1L)) {
    stacked_raw <- reg_raw_final
  } else if (identical(method, 2L)) {
    stacked_raw <- dim_raw_final
  } else {
    reg_bad <- !is.finite(reg_raw_final)
    dim_bad <- !is.finite(dim_raw_final)
    if (any(reg_bad & dim_bad)) {
      stop("Both final component models failed for at least one test observation.",
           call. = FALSE)
    }
    reg_use <- reg_raw_final
    dim_use <- dim_raw_final
    reg_use[reg_bad] <- dim_use[reg_bad]
    dim_use[dim_bad] <- reg_use[dim_bad]
    stacked_raw <- component_weights["reg"] * reg_use +
      component_weights["dim.red"] * dim_use
  }
  
  stacked_output <- if (is_class) {
    .decode_codes(.round_codes(stacked_raw, probability.threshold))
  } else {
    stacked_raw
  }
  
  reg_pred_int <- if (is_class) {
    .decode_interval(reg_pred_int_raw, reg_component_threshold)
  } else {
    reg_pred_int_raw
  }
  dim_pred_int <- if (is_class) {
    .decode_interval(dim_pred_int_raw, dim_component_threshold)
  } else {
    dim_pred_int_raw
  }
  
  if (is.null(pred.int)) {
    stacked_pred_int_raw <- NULL
  } else if (identical(method, 1L)) {
    stacked_pred_int_raw <- reg_pred_int_raw
  } else if (identical(method, 2L)) {
    stacked_pred_int_raw <- dim_pred_int_raw
  } else if (is.null(reg_pred_int_raw)) {
    stacked_pred_int_raw <- dim_pred_int_raw
  } else if (is.null(dim_pred_int_raw)) {
    stacked_pred_int_raw <- reg_pred_int_raw
  } else {
    reg_pi <- as.data.frame(reg_pred_int_raw, check.names = FALSE)
    dim_pi <- as.data.frame(dim_pred_int_raw, check.names = FALSE)
    if (!identical(dim(reg_pi), dim(dim_pi))) {
      stop("Component prediction intervals have incompatible dimensions.",
           call. = FALSE)
    }
    stacked_pred_int_raw <- as.data.frame(
      component_weights["reg"] * as.matrix(reg_pi) +
        component_weights["dim.red"] * as.matrix(dim_pi),
      check.names = FALSE
    )
    names(stacked_pred_int_raw) <- names(reg_pi)
  }
  
  stacked_pred_int <- if (is_class) {
    .decode_interval(stacked_pred_int_raw, probability.threshold)
  } else {
    stacked_pred_int_raw
  }
  
  result <- list(
    OBJfn.reg = unname(as.numeric(reg_best_score)),
    NNS.reg.n.best = unname(as.integer(reg_best_k)),
    probability.threshold = unname(as.numeric(probability.threshold)),
    OBJfn.dim.red = unname(as.numeric(dim_best_score)),
    NNS.dim.red.threshold = unname(as.numeric(dim_threshold_report)),
    reg = reg_output,
    reg.pred.int = .NNS.df(reg_pred_int),
    dim.red = dim_output,
    dim.red.pred.int = .NNS.df(dim_pred_int),
    stack = stacked_output,
    pred.int = .NNS.df(stacked_pred_int),
    weights = component_weights,
    class.levels = if (is_class) class_values else NULL
  )
  
  .NNS.out(result)
}