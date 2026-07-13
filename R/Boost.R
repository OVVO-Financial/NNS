#' NNS Boost
#'
#' Ensemble feature-selection method using \link{NNS.reg} as the base learner.
#'
#' @param IVs.train a vector, matrix, or data frame of numeric, logical, character,
#'   or factor predictors.
#' @param DV.train a numeric, logical, character, or factor response with one value
#'   per row of \code{IVs.train}.
#' @param IVs.test a vector, matrix, or data frame with the same predictor columns
#'   as \code{IVs.train}. If \code{NULL}, \code{IVs.train} is used.
#' @param type \code{NULL} (default) for regression, or \code{"CLASS"} for
#'   classification. Factor, character, and logical responses automatically select
#'   classification.
#' @param depth integer, \code{NULL}, or \code{"max"}; passed to the \code{order}
#'   argument of \link{NNS.reg}.
#' @param learner.trials positive integer; maximum number of feature subsets used to
#'   estimate the learner threshold. If every possible subset can be evaluated within
#'   this limit, all subsets are evaluated.
#' @param epochs non-negative integer; number of weighted feature subsets evaluated
#'   after the learner stage. Defaults to \code{2 * length(DV.train)}. Set to zero to
#'   use the surviving learner subsets directly.
#' @param CV.size numeric in \code{(0, 1)}; validation fraction for non-time-series
#'   data. If \code{NULL}, one value between 0.2 and 1/3 is drawn under the local seed.
#' @param balance logical; if \code{TRUE}, down- and up-sampling are applied only to
#'   the fitting portion of each split. Validation observations are never resampled.
#' @param ts.test positive integer smaller than the training sample size; the final
#'   \code{ts.test} observations are used as the chronological validation block.
#' @param threshold finite numeric scalar or \code{NULL}; objective cutoff used to
#'   retain feature subsets. If \code{NULL}, the lower quartile is used for a
#'   minimization objective and the upper quartile for a maximization objective.
#' @param obj.fn expression using the names \code{predicted} and \code{actual}.
#'   Defaults to sum of squared errors. For explicit classification, the untouched
#'   default is replaced by mean classification accuracy.
#' @param objective one of \code{"min"} or \code{"max"}; defaults to \code{"min"}.
#' @param extreme logical; if \code{TRUE}, use the best learner score rather than a
#'   quartile cutoff.
#' @param features.only logical; return only feature weights and frequencies.
#' @param feature.importance logical; plot up to the ten most frequently retained
#'   features.
#' @param pred.int numeric in \code{(0, 1)} or \code{NULL}; prediction interval level
#'   passed to the final \link{NNS.reg} fit.
#' @param status logical; print progress messages.
#' @param seed integer or \code{NULL}; local random seed. The caller's RNG state is
#'   restored when the function exits.
#'
#' @return A list containing \code{results}, \code{pred.int},
#'   \code{feature.weights}, and \code{feature.frequency}. With
#'   \code{features.only = TRUE}, only the last two elements are returned.
#'
#' @note
#' \itemize{
#'   \item Numeric class labels are returned on their original scale. Factor,
#'     character, and logical responses retain the historical integer-code output.
#'   \item Categorical predictors are aligned to training levels. Unseen test levels
#'     cause an explicit error rather than silent recoding.
#'   \item Incorporate an objective from another package with, for example,
#'     \code{obj.fn = expression(Metrics::mape(actual, predicted))} and
#'     \code{objective = "min"}.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#'   \doi{10.2139/ssrn.2864711}
#' @examples
#' \dontrun{
#' a <- NNS.boost(
#'   iris[1:140, 1:4], iris[1:140, 5],
#'   IVs.test = iris[141:150, 1:4],
#'   epochs = 100, learner.trials = 100,
#'   type = "CLASS", balance = TRUE
#' )
#'
#' mean(a$results == as.numeric(iris[141:150, 5]))
#' }
#'
#' @export
NNS.boost <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      type = NULL,
                      depth = NULL,
                      learner.trials = 100,
                      epochs = NULL,
                      CV.size = NULL,
                      balance = FALSE,
                      ts.test = NULL,
                      threshold = NULL,
                      obj.fn = expression(sum((predicted - actual)^2)),
                      objective = "min",
                      extreme = FALSE,
                      features.only = FALSE,
                      feature.importance = TRUE,
                      pred.int = NULL,
                      status = TRUE,
                      seed = 123L) {

  # ---------------------------------------------------------------------------
  # Local validation and coercion helpers
  # ---------------------------------------------------------------------------

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
      stop(sprintf("[%s] must be an integer >= %d.", name, minimum), call. = FALSE)
    }
    as.integer(x)
  }

  .as_train_frame <- function(x) {
    if (any(class(x) %in% c("tbl", "data.table"))) x <- as.data.frame(x)
    if (is.null(dim(x))) x <- data.frame(X1 = x, check.names = FALSE)
    x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) < 1L) stop("[IVs.train] must contain at least one predictor.", call. = FALSE)
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
    had_column_names <- !is.null(dim(x)) && !is.null(colnames(x)) &&
      length(colnames(x)) == NCOL(x) && all(colnames(x) != "")
    if (any(class(x) %in% c("tbl", "data.table"))) {
      had_column_names <- !is.null(names(x)) && all(names(x) != "")
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
          # Align a named test row by the training predictor names rather
          # than silently renaming positionally supplied values.
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
        stop("A vector [IVs.test] must contain one complete test row, unless the training data have one predictor.",
             call. = FALSE)
      }
    } else {
      x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    }

    if (ncol(x) != p) {
      stop("[IVs.test] must have the same number of predictors as [IVs.train].",
           call. = FALSE)
    }

    if (!had_column_names) {
      names(x) <- train_names
    } else {
      if (anyDuplicated(names(x))) {
        stop("[IVs.test] predictor names must be unique.", call. = FALSE)
      }
      missing_names <- setdiff(train_names, names(x))
      extra_names <- setdiff(names(x), train_names)
      if (length(missing_names) || length(extra_names)) {
        stop(sprintf(
          "[IVs.test] columns must exactly match [IVs.train]. Missing: %s; extra: %s.",
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
          stop(sprintf("Test predictor [%s] must also be a Date.", nm), call. = FALSE)
        }
        train[[j]] <- as.numeric(tr)
        test[[j]] <- as.numeric(te)
      } else if (inherits(tr, c("POSIXct", "POSIXlt"))) {
        if (!inherits(te, c("POSIXct", "POSIXlt"))) {
          stop(sprintf("Test predictor [%s] must also be a date-time value.", nm), call. = FALSE)
        }
        train[[j]] <- as.numeric(tr)
        test[[j]] <- as.numeric(te)
      } else if (is.factor(tr) || is.character(tr)) {
        tr_chr <- as.character(tr)
        te_chr <- as.character(te)
        lev <- if (is.factor(tr)) levels(droplevels(tr)) else sort(unique(tr_chr))
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
          stop(sprintf("Test predictor [%s] must also be logical.", nm), call. = FALSE)
        }
      } else if (is.numeric(tr) || is.integer(tr)) {
        if (!(is.numeric(te) || is.integer(te))) {
          stop(sprintf("Test predictor [%s] must be numeric.", nm), call. = FALSE)
        }
        train[[j]] <- as.numeric(tr)
        test[[j]] <- as.numeric(te)
      } else {
        stop(sprintf("Unsupported predictor type for [%s].", nm), call. = FALSE)
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
    old <- if (existed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    function() {
      if (existed) {
        assign(".Random.seed", old, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  })
  on.exit(.restore_rng(), add = TRUE)

  # ---------------------------------------------------------------------------
  # Validate arguments and establish response coding
  # ---------------------------------------------------------------------------

  balance <- .scalar_logical(balance, "balance")
  extreme <- .scalar_logical(extreme, "extreme")
  features.only <- .scalar_logical(features.only, "features.only")
  feature.importance <- .scalar_logical(feature.importance, "feature.importance")
  status <- .scalar_logical(status, "status")

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

  objective <- match.arg(tolower(as.character(objective)[1L]), c("min", "max"))

  if (!is.null(type)) {
    type <- match.arg(tolower(as.character(type)[1L]), "class")
  }
  if (balance && is.null(type)) {
    warning("type = 'CLASS' selected because balance = TRUE.", call. = FALSE)
    type <- "class"
  }

  if (!is.null(depth)) {
    if (is.character(depth)) {
      if (length(depth) != 1L || tolower(depth) != "max") {
        stop("[depth] must be NULL, a positive integer, or 'max'.", call. = FALSE)
      }
      depth <- "max"
    } else {
      depth <- .scalar_integer(depth, "depth", minimum = 1L)
    }
  }

  learner.trials <- .scalar_integer(learner.trials, "learner.trials", minimum = 1L)
  epochs <- .scalar_integer(epochs, "epochs", minimum = 0L, allow_null = TRUE)

  if (!is.null(CV.size)) {
    if (!is.numeric(CV.size) || length(CV.size) != 1L || !is.finite(CV.size) ||
        CV.size <= 0 || CV.size >= 1) {
      stop("[CV.size] must be a finite scalar strictly between 0 and 1.", call. = FALSE)
    }
    CV.size <- as.numeric(CV.size)
  }

  ts.test <- .scalar_integer(ts.test, "ts.test", minimum = 1L, allow_null = TRUE)

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
      stop("[threshold] must be a finite numeric scalar or NULL.", call. = FALSE)
    }
    threshold <- as.numeric(threshold)
  }

  if (!is.null(pred.int)) {
    if (!is.numeric(pred.int) || length(pred.int) != 1L || !is.finite(pred.int) ||
        pred.int <= 0 || pred.int >= 1) {
      stop("[pred.int] must be a finite scalar strictly between 0 and 1.", call. = FALSE)
    }
    pred.int <- as.numeric(pred.int)
  }

  x <- .as_train_frame(IVs.train)

  if (any(class(DV.train) %in% c("tbl", "data.table"))) {
    DV.train <- as.vector(unlist(DV.train))
  }
  if (is.data.frame(DV.train) || is.matrix(DV.train)) {
    if (NCOL(DV.train) != 1L) {
      stop("[DV.train] must contain exactly one response column.", call. = FALSE)
    }
    DV.train <- as.vector(unlist(DV.train))
  }
  if (length(DV.train) != nrow(x)) {
    stop("nrow(IVs.train) must equal length(DV.train).", call. = FALSE)
  }
  if (length(DV.train) < 4L) {
    stop("NNS.boost requires at least four training observations.", call. = FALSE)
  }
  if (anyNA(DV.train)) stop("[DV.train] contains missing values.", call. = FALSE)
  if (is.numeric(DV.train) && any(!is.finite(DV.train))) {
    stop("[DV.train] contains non-finite values.", call. = FALSE)
  }

  response_was_numeric <- is.numeric(DV.train) || is.integer(DV.train)
  auto_class <- is.factor(DV.train) || is.character(DV.train) || is.logical(DV.train)
  if (auto_class && is.null(type)) type <- "class"
  is_class <- identical(type, "class")

  original_response <- DV.train
  class_values <- NULL

  if (is_class) {
    if (response_was_numeric) {
      class_values <- sort(unique(as.numeric(DV.train)))
      y <- match(as.numeric(DV.train), class_values)
    } else {
      class_factor <- if (is.factor(DV.train)) droplevels(DV.train) else factor(DV.train)
      class_values <- levels(class_factor)
      y <- as.integer(class_factor)
    }
    y <- as.numeric(y)
    if (length(unique(y)) < 2L) {
      stop("Classification requires at least two response classes.", call. = FALSE)
    }

    if (identical(obj.fn, expression(sum((predicted - actual)^2)))) {
      obj.fn <- expression(mean(predicted == actual))
      objective <- "max"
    }
  } else {
    if (!(is.numeric(DV.train) || is.integer(DV.train))) {
      stop("A nonnumeric response requires type = 'CLASS'.", call. = FALSE)
    }
    y <- as.numeric(DV.train)
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
  n_features <- ncol(x)

  if (!is.null(ts.test) && ts.test >= n_obs) {
    stop("[ts.test] must be smaller than the number of training observations.",
         call. = FALSE)
  }

  if (is.null(epochs)) epochs <- as.integer(2L * n_obs)
  cv_fraction <- if (is.null(CV.size)) stats::runif(1L, 0.2, 1 / 3) else CV.size

  # ---------------------------------------------------------------------------
  # Scoring, splitting, balancing, and prediction helpers
  # ---------------------------------------------------------------------------

  .score <- function(predicted, actual) {
    if (length(predicted) != length(actual)) {
      stop("The objective received predicted and actual vectors of different lengths.",
           call. = FALSE)
    }
    value <- eval(obj.fn,
                  envir = list(predicted = predicted, actual = actual),
                  enclos = parent.frame())
    if (!is.numeric(value) || length(value) != 1L) {
      stop("[obj.fn] must return one numeric value.", call. = FALSE)
    }
    value <- as.numeric(value)
    if (!is.finite(value)) return(NA_real_)
    value
  }

  .central_value <- function(v, classification = FALSE) {
    finite <- v[is.finite(v)]
    if (!length(finite)) return(NA_real_)
    out <- if (classification) gravity_class(finite) else gravity(finite)
    if (!is.finite(out)) out <- if (classification) mode_class(finite) else mean(finite)
    as.numeric(out)
  }

  .sanitize_predictions <- function(predicted, fallback_y) {
    predicted <- as.numeric(predicted)
    bad <- !is.finite(predicted)
    if (any(bad)) {
      replacement <- .central_value(predicted[!bad], is_class)
      if (!is.finite(replacement)) replacement <- .central_value(fallback_y, is_class)
      if (!is.finite(replacement)) {
        stop("NNS.reg returned no finite predictions.", call. = FALSE)
      }
      predicted[bad] <- replacement
    }
    if (is_class) {
      predicted <- pmin(pmax(predicted, min(y)), max(y))
      predicted <- ifelse(predicted %% 1 < 0.5, floor(predicted), ceiling(predicted))
    }
    predicted
  }

  .has_all_classes <- function(train_y) {
    !is_class || identical(sort(unique(train_y)), sort(unique(y)))
  }

  .random_validation_index <- function() {
    size <- max(1L, min(n_obs - 1L, as.integer(round(cv_fraction * n_obs))))
    for (attempt in seq_len(200L)) {
      idx <- sort(sample.int(n_obs, size = size, replace = FALSE))
      if (.has_all_classes(y[-idx])) return(idx)
    }
    stop("Unable to create a validation split retaining every class in training. Reduce CV.size or provide more observations per class.",
         call. = FALSE)
  }

  validation_index <- if (!is.null(ts.test)) {
    idx <- seq.int(n_obs - ts.test + 1L, n_obs)
    if (!.has_all_classes(y[-idx])) {
      stop("The chronological training prefix does not contain every response class.",
           call. = FALSE)
    }
    idx
  } else {
    .random_validation_index()
  }

  .balance_training <- function(train_x, train_y) {
    if (!balance) return(list(x = train_x, y = train_y))

    groups <- split(seq_along(train_y), train_y)
    if (length(groups) < 2L || any(lengths(groups) == 0L)) {
      stop("Balancing requires at least two non-empty classes in the fitting split.",
           call. = FALSE)
    }

    smallest <- min(lengths(groups))
    largest <- max(lengths(groups))

    down_idx <- unlist(lapply(groups, function(g) sample(g, smallest, replace = FALSE)),
                       use.names = FALSE)
    up_idx <- unlist(lapply(groups, function(g) sample(g, largest, replace = TRUE)),
                     use.names = FALSE)

    idx <- sort(c(down_idx, up_idx))
    list(x = train_x[idx, , drop = FALSE], y = train_y[idx])
  }

  .fit_subset <- function(feature_index, train_index, test_index) {
    train_x <- x[train_index, feature_index, drop = FALSE]
    train_y <- y[train_index]
    test_x <- x[test_index, feature_index, drop = FALSE]

    balanced <- .balance_training(train_x, train_y)

    fit <- suppressWarnings(
      NNS.reg(
        balanced$x,
        balanced$y,
        point.est = test_x,
        dim.red.method = if (ncol(balanced$x) > 1L) "equal" else NULL,
        plot = FALSE,
        residual.plot = FALSE,
        order = depth,
        ncores = 1,
        type = type,
        point.only = TRUE
      )
    )

    .sanitize_predictions(fit$Point.est, balanced$y)
  }

  # ---------------------------------------------------------------------------
  # Generate learner feature subsets
  # ---------------------------------------------------------------------------

  total_sets <- if (n_features <= 30L) 2^n_features - 1 else Inf
  exhaustive <- is.finite(total_sets) && total_sets <= learner.trials

  if (exhaustive) {
    test.features <- unlist(
      lapply(seq_len(n_features), function(k) {
        as.list(as.data.frame(utils::combn(n_features, k)))
      }),
      recursive = FALSE
    )
  } else {
    target_trials <- if (is.finite(total_sets)) {
      min(learner.trials, as.integer(total_sets))
    } else {
      learner.trials
    }

    test.features <- vector("list", target_trials)
    seen <- new.env(hash = TRUE, parent = emptyenv())
    filled <- 0L
    attempts <- 0L
    max_attempts <- max(1000L, target_trials * 200L)

    while (filled < target_trials && attempts < max_attempts) {
      attempts <- attempts + 1L
      k <- sample.int(n_features, 1L)
      candidate <- sort(sample.int(n_features, k, replace = FALSE))
      key <- paste(candidate, collapse = ",")
      if (!exists(key, envir = seen, inherits = FALSE)) {
        filled <- filled + 1L
        test.features[[filled]] <- candidate
        assign(key, TRUE, envir = seen)
      }
    }

    if (filled < target_trials) {
      test.features <- test.features[seq_len(filled)]
      warning("Fewer unique feature subsets were generated than requested.", call. = FALSE)
    }
  }

  learner_count <- length(test.features)
  results <- rep(NA_real_, learner_count)
  train_index <- setdiff(seq_len(n_obs), validation_index)
  actual <- y[validation_index]

  for (i in seq_len(learner_count)) {
    if (status) {
      message("Current Threshold Iterations Remaining = ",
              learner_count - i, " ", "\r", appendLF = FALSE)
    }

    predicted <- .fit_subset(test.features[[i]], train_index, validation_index)
    results[i] <- .score(predicted, actual)
  }

  finite_results <- which(is.finite(results))
  if (!length(finite_results)) {
    stop("No learner trial produced a finite objective value.", call. = FALSE)
  }

  supplied_threshold <- !is.null(threshold)
  if (!supplied_threshold) {
    if (extreme) {
      threshold <- if (objective == "max") {
        max(results[finite_results])
      } else {
        min(results[finite_results])
      }
    } else {
      threshold <- as.numeric(stats::quantile(
        results[finite_results],
        probs = if (objective == "max") 0.75 else 0.25,
        names = FALSE,
        type = 2
      ))
    }
  }

  if (status) {
    message(
      paste0("\nLearner Accuracy Threshold = ",
             format(threshold, digits = 4, nsmall = 2), "           "),
      appendLF = TRUE
    )
  }

  passes_learner <- is.finite(results) & if (objective == "max") {
    results >= threshold
  } else {
    results <= threshold
  }

  reduced.test.features <- test.features[passes_learner]

  if (!length(reduced.test.features)) {
    if (supplied_threshold) {
      if (objective == "min") {
        stop("No learner subset met [threshold]; increase the threshold.", call. = FALSE)
      } else {
        stop("No learner subset met [threshold]; reduce the threshold.", call. = FALSE)
      }
    }

    best_index <- if (objective == "min") {
      finite_results[which.min(results[finite_results])]
    } else {
      finite_results[which.max(results[finite_results])]
    }
    reduced.test.features <- list(test.features[[best_index]])
  }

  # Preserve repeated survivor counts. Features with no survivor count retain a very
  # small probability so every predictor can still be explored during epochs.
  feature_count <- tabulate(unlist(reduced.test.features), nbins = n_features)
  feature_prob <- as.numeric(feature_count)
  feature_prob <- feature_prob + max(1, sum(feature_prob)) * 1e-12
  feature_prob <- feature_prob / sum(feature_prob)

  # ---------------------------------------------------------------------------
  # Weighted epoch stage
  # ---------------------------------------------------------------------------

  keeper.features <- list()

  if (!exhaustive && epochs > 0L) {
    keeper.features <- vector("list", epochs)

    for (j in seq_len(epochs)) {
      if (status) {
        message("% of epochs = ",
                format(j / epochs, digits = 3, nsmall = 2),
                "     ", "\r", appendLF = FALSE)
      }

      k <- sample.int(n_features, 1L)
      features_j <- sort(sample.int(
        n_features,
        size = k,
        replace = FALSE,
        prob = feature_prob
      ))

      predicted <- .fit_subset(features_j, train_index, validation_index)
      new_result <- .score(predicted, actual)

      passes <- is.finite(new_result) && if (objective == "max") {
        new_result >= threshold
      } else {
        new_result <= threshold
      }

      keeper.features[[j]] <- if (passes) features_j else NULL
    }

    keeper.features <- keeper.features[!vapply(keeper.features, is.null, logical(1L))]
  } else {
    keeper.features <- reduced.test.features
  }

  if (!length(keeper.features)) {
    if (supplied_threshold) {
      if (objective == "min") {
        stop("No epoch subset met [threshold]; increase the threshold.", call. = FALSE)
      } else {
        stop("No epoch subset met [threshold]; reduce the threshold.", call. = FALSE)
      }
    }

    best_index <- if (objective == "min") {
      finite_results[which.min(results[finite_results])]
    } else {
      finite_results[which.max(results[finite_results])]
    }
    keeper.features <- list(test.features[[best_index]])
  }

  plot.table <- table(factor(
    unlist(keeper.features),
    levels = seq_len(n_features)
  ))
  names(plot.table) <- names(x)
  plot.table <- sort(plot.table[plot.table > 0L], decreasing = TRUE)

  if (features.only) {
    return(.NNS.out(list(
      feature.weights = plot.table / sum(plot.table),
      feature.frequency = plot.table
    )))
  }

  if (status) message("\nGenerating Final Estimate", "\r", appendLF = TRUE)

  # ---------------------------------------------------------------------------
  # Build a training-fitted numeric design and frequency-weighted X*.
  # Categorical features are one-hot encoded with training levels only. Each active
  # dummy receives its original feature's weight, so a categorical feature's total
  # per-row contribution remains equal to that feature weight.
  # ---------------------------------------------------------------------------

  .numeric_design <- function(train, test) {
    train_blocks <- list()
    test_blocks <- list()
    source_feature <- character()

    for (j in seq_along(train)) {
      nm <- names(train)[j]
      tr <- train[[j]]
      te <- test[[j]]

      if (is.factor(tr)) {
        lev <- levels(tr)
        if (!length(lev)) {
          stop(sprintf("Predictor [%s] has no observed levels.", nm), call. = FALSE)
        }
        tr_block <- vapply(lev, function(level) as.numeric(as.character(tr) == level),
                           numeric(length(tr)))
        te_block <- vapply(lev, function(level) as.numeric(as.character(te) == level),
                           numeric(length(te)))
        if (is.null(dim(tr_block))) tr_block <- matrix(tr_block, ncol = 1L)
        if (is.null(dim(te_block))) te_block <- matrix(te_block, ncol = 1L)
        colnames(tr_block) <- paste0(nm, "__", make.names(lev, unique = TRUE))
        colnames(te_block) <- colnames(tr_block)
        source_feature <- c(source_feature, rep(nm, length(lev)))
      } else {
        tr_block <- matrix(as.numeric(tr), ncol = 1L,
                           dimnames = list(NULL, nm))
        te_block <- matrix(as.numeric(te), ncol = 1L,
                           dimnames = list(NULL, nm))
        source_feature <- c(source_feature, nm)
      }

      train_blocks[[j]] <- tr_block
      test_blocks[[j]] <- te_block
    }

    train_matrix <- do.call(cbind, train_blocks)
    test_matrix <- do.call(cbind, test_blocks)
    storage.mode(train_matrix) <- "double"
    storage.mode(test_matrix) <- "double"

    list(train = train_matrix, test = test_matrix, source = source_feature)
  }

  design <- .numeric_design(x, z)

  train_min <- apply(design$train, 2L, min)
  train_max <- apply(design$train, 2L, max)
  train_range <- train_max - train_min
  train_range[!is.finite(train_range) | train_range == 0] <- 1

  train_norm <- sweep(design$train, 2L, train_min, "-")
  train_norm <- sweep(train_norm, 2L, train_range, "/")
  test_norm <- sweep(design$test, 2L, train_min, "-")
  test_norm <- sweep(test_norm, 2L, train_range, "/")

  feature_weights <- as.numeric(plot.table / sum(plot.table))
  names(feature_weights) <- names(plot.table)
  coef_design <- feature_weights[design$source]
  coef_design[is.na(coef_design)] <- 0

  if (!any(coef_design > 0)) {
    stop("No positive feature weights were available for the final estimate.",
         call. = FALSE)
  }

  xstar_train <- as.numeric(train_norm %*% coef_design)
  xstar_test <- as.numeric(test_norm %*% coef_design)

  if (any(!is.finite(xstar_train)) || any(!is.finite(xstar_test))) {
    stop("The final synthetic predictor contains non-finite values.", call. = FALSE)
  }

  xstar_frame <- function(v) {
    data.frame(xstar = v, xstar2 = v, check.names = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Select n.best locally. This avoids balancing before CV and avoids inheriting
  # the reversed ts.test split currently present in NNS.stack.
  # ---------------------------------------------------------------------------

  .final_validation_splits <- function() {
    if (!is.null(ts.test)) {
      return(list(validation_index))
    }

    splits <- vector("list", 5L)
    for (b in seq_len(5L)) splits[[b]] <- .random_validation_index()
    splits
  }

  final_splits <- .final_validation_splits()
  minimum_train_size <- min(vapply(final_splits, function(idx) n_obs - length(idx), integer(1L)))
  k_small <- max(1L, floor(sqrt(minimum_train_size)))
  k_candidates <- unique(c(seq_len(k_small), minimum_train_size))
  k_scores <- rep(NA_real_, length(k_candidates))

  for (ki in seq_along(k_candidates)) {
    k_value <- k_candidates[ki]
    split_scores <- rep(NA_real_, length(final_splits))

    for (b in seq_along(final_splits)) {
      test_idx <- final_splits[[b]]
      train_idx <- setdiff(seq_len(n_obs), test_idx)

      fold_train_x <- xstar_frame(xstar_train[train_idx])
      fold_train_y <- y[train_idx]
      fold_test_x <- xstar_frame(xstar_train[test_idx])

      fold_balanced <- .balance_training(fold_train_x, fold_train_y)

      fold_fit <- suppressWarnings(
        NNS.reg(
          fold_balanced$x,
          fold_balanced$y,
          point.est = fold_test_x,
          plot = FALSE,
          residual.plot = FALSE,
          n.best = min(k_value, nrow(fold_balanced$x)),
          order = depth,
          ncores = 1,
          type = type,
          factor.2.dummy = FALSE,
          dist = "L2",
          point.only = TRUE
        )
      )

      fold_pred <- .sanitize_predictions(fold_fit$Point.est, fold_balanced$y)
      split_scores[b] <- .score(fold_pred, y[test_idx])
    }

    finite_split_scores <- split_scores[is.finite(split_scores)]
    if (length(finite_split_scores)) k_scores[ki] <- mean(finite_split_scores)

    if (status) {
      message(
        sprintf(
          "Current NNS.reg(..., n.best = %d) | mean eval(obj.fn) = %s | Iterations Remaining = %d",
          k_value,
          if (is.finite(k_scores[ki])) format(k_scores[ki], digits = 6) else "NA",
          length(k_candidates) - ki
        )
      )
    }
  }

  finite_k <- which(is.finite(k_scores))
  if (!length(finite_k)) {
    stop("No n.best candidate produced a finite objective value.", call. = FALSE)
  }

  best_k_index <- if (objective == "min") {
    finite_k[which.min(k_scores[finite_k])]
  } else {
    finite_k[which.max(k_scores[finite_k])]
  }
  best_k <- k_candidates[best_k_index]

  final_train <- xstar_frame(xstar_train)
  final_test <- xstar_frame(xstar_test)
  final_balanced <- .balance_training(final_train, y)

  final_fit <- suppressWarnings(
    NNS.reg(
      final_balanced$x,
      final_balanced$y,
      point.est = final_test,
      plot = FALSE,
      residual.plot = FALSE,
      n.best = min(best_k, nrow(final_balanced$x)),
      order = depth,
      ncores = 1,
      type = type,
      factor.2.dummy = FALSE,
      dist = "L2",
      point.only = FALSE,
      confidence.interval = pred.int
    )
  )

  estimates_code <- .sanitize_predictions(final_fit$Point.est, final_balanced$y)
  pred_int_out <- final_fit$pred.int

  if (is_class) {
    estimates_code <- pmin(pmax(estimates_code, 1L), length(class_values))
    estimates_code <- as.integer(round(estimates_code))

    if (response_was_numeric) {
      estimates <- class_values[estimates_code]
    } else {
      estimates <- estimates_code
    }

    if (!is.null(pred_int_out)) {
      pred_int_out <- as.data.frame(pred_int_out)
      pred_int_out[] <- lapply(pred_int_out, function(v) {
        code <- pmin(pmax(as.integer(round(v)), 1L), length(class_values))
        if (response_was_numeric) class_values[code] else code
      })
    }
  } else {
    estimates <- estimates_code
  }

  if (feature.importance) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    top <- tail(sort(plot.table, decreasing = FALSE),
                min(length(plot.table), 10L))
    label_margin <- max(graphics::strwidth(names(top), "inch") + 0.4,
                        na.rm = TRUE)
    graphics::par(mai = c(1.0, label_margin, 0.8, 0.5))
    graphics::barplot(
      top,
      horiz = TRUE,
      col = "steelblue",
      main = "Feature Frequency in Final Estimate",
      xlab = "Frequency",
      las = 1
    )
  }

  .NNS.out(list(
    results = estimates,
    pred.int = pred_int_out,
    feature.weights = plot.table / sum(plot.table),
    feature.frequency = plot.table
  ))
}
