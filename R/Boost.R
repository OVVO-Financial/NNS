#' NNS Boost
#'
#' Ensemble method for classification using the NNS multivariate regression \link{NNS.reg} as the base learner instead of trees.
#'
#' @param IVs.train a matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimensions to \code{(IVs.train)}.
#' @param IVs.test a matrix or data frame of variables of numeric or factor data types with compatible dimensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.  Columns are matched to \code{(IVs.train)} by name when the two share the same predictor names, and positionally when they share no names at all (as in \code{cbind(test.x_1, test.x_2)} against \code{cbind(x_1, x_2)}).  Names that only partly overlap the training predictors are ambiguous and return an error.
#' @param type \code{NULL} (default).  To perform a classification of discrete integer classes from factor target variable \code{(DV.train)} with a base category of 1, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.
#' @param depth options: (integer, NULL, "max"); \code{(depth = NULL)}(default) Specifies the \code{order} parameter in the \link{NNS.reg} routine, assigning a number of splits in the regressors, analogous to tree depth.
#' @param learner.trials integer; 100 (default) Sets the number of trials to obtain an accuracy \code{threshold} level.  If the number of all possible feature combinations is less than selected value, the minimum of the two values will be used.
#' @param epochs integer; \code{2*length(DV.train)} (default) Number of repeated holdout re-evaluations of the learner-trial feature subsets that pass the accuracy threshold.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size.  Defaults to a random value between 0.2 and 0.33 for a random sampling of the training set.
#' @param balance logical; \code{FALSE} (default) Uses both up and down sampling to balance the classes.  \code{type="CLASS"} required.
#' @param ts.test integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.
#' @param threshold numeric [0, 1]; \code{NULL} (default) Probability supplied to \link{LPM.VaR} over the learner-trial objective distribution to determine the objective cutoff for keeping feature combinations.  Defaults to 0.80 when \code{objective = "max"} and 0.20 when \code{objective = "min"}.  It is not a literal objective-score cutoff.
#' @param obj.fn expression;
#' \code{expression( sum((predicted - actual)^2) )} (default) Sum of squared errors is the default objective function.  Any \code{expression(...)} using the specific terms \code{predicted} and \code{actual} can be used.  Automatically selects an accuracy measure when \code{(type = "CLASS")}.
#' @param objective options: ("min", "max") \code{"max"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param extreme logical; \code{FALSE} (default) Sets the \link{LPM.VaR} probability to 1 (0) for maximization (minimization) \code{objective}, i.e. the most extreme learner-trial objective value becomes the cutoff.  Overrides \code{threshold}.
#' @param features.only logical; \code{FALSE} (default) Returns only the final feature loadings along with the final feature frequencies.
#' @param feature.importance logical; \code{TRUE} (default) Draws a two-panel diagnostic: the learner-trial objective distribution with its \link{LPM.VaR} cutoff, and the frequency of features used in the final estimate.
#' @param pred.int numeric [0,1]; \code{NULL} (default) Returns the associated prediction intervals for the final estimate.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param seed Optional integer random seed used for reproducible resampling, fold construction, and stochastic fitting steps. If `NULL`, the current random-number-generator state is used.
#' @param dist options:(NULL, "NNS", "L1", "L2", "FACTOR") the method of distance calculation passed to delegated \link{NNS.reg} and \link{NNS.stack} calls. \code{dist = NULL} is the default and selects the native blended NNS distance; \code{dist = "NNS"} is an explicit alias for the default.
#' @param folds integer; 5 (default) Number of cross-validation \code{folds} passed to the final \link{NNS.stack} call.
#'
#' @return Returns a vector of fitted values for the dependent variable test set \code{$results}, prediction intervals \code{$pred.int}, the final feature loadings \code{$feature.weights}, final feature frequencies \code{$feature.frequency}, and (for classification) the class labels \code{$class.levels}.  Classification results are numeric: a factor or character \code{DV.train} yields integer class codes with a base category of 1 (label recoverable as \code{class.levels[results]}), and a numeric \code{DV.train} yields its original numeric class values.
#'
#' @note
#' \itemize{
#' \item{} Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].  The response variable base category should be 1 for classification problems.
#'
#' \item{} Incorporate any objective function from external packages (such as \code{Metrics::mape}) via \code{NNS.boost(..., obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")}
#'}
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"  \doi{10.2139/ssrn.2864711}
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  a <- NNS.boost(iris[1:140, 1:4], iris[1:140, 5],
#'  IVs.test = iris[141:150, 1:4],
#'  epochs = 100, learner.trials = 100,
#'  type = "CLASS", depth = NULL, balance = TRUE)
#'
#'  ## Test accuracy
#'  mean(a$results == as.numeric(iris[141:150, 5]))
#'
#'  ## Recover the labels
#'  a$class.levels[a$results]
#'  }
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
                      seed = 123L,
                      dist = NULL,
                      folds = 5) {
  dist <- .nns_reg_validate_dist(dist)
  # ---------------------------------------------------------------------------
  
  # Local validation and coercion helpers
  
  # ---------------------------------------------------------------------------
  
  
  .scalar_logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("[%s] must be TRUE or FALSE.", name), call. = FALSE)
    }
    x
  }
  
  .scalar_integer <- function(x,
                              name,
                              minimum = 0L,
                              allow_null = FALSE) {
    if (allow_null && is.null(x))
      return(NULL)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x)) {
      stop(sprintf("[%s] must be an integer >= %d.", name, minimum),
           call. = FALSE)
    }
    as.integer(x)
  }
  
  .as_train_frame <- function(x) {
    if (any(class(x) %in% c("tbl", "data.table")))
      x <- as.data.frame(x)
    if (is.null(dim(x)))
      x <- data.frame(X1 = x, check.names = FALSE)
    x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) < 1L)
      stop("[IVs.train] must contain at least one predictor.", call. = FALSE)
    if (is.null(names(x)) || any(names(x) == "")) {
      names(x) <- paste0("X", seq_len(ncol(x)))
    }
    # De-duplicate repeated predictor names with make.unique() exactly as
    
    # NNS.reg's .nns_reg_as_frame() does, so the cbind(x, x) dimension trick
    
    # works in NNS.boost instead of erroring. c("x", "x") becomes
    
    # c("x", "x.1").
    
    names(x) <- make.unique(names(x), sep = ".")
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
        if (!is.null(supplied) && all(nzchar(supplied))) {
          # Normalize duplicate names exactly as the training frame does, then
          # align a named test row by the training predictor names only when
          # the two name sets describe the same predictors.
          supplied <- make.unique(supplied, sep = ".")
          ordering <- .nns_match_predictor_names(supplied, train_names, "IVs.test")
          if (!is.null(ordering)) x <- x[ordering]
        }
        x <- as.data.frame(as.list(x),
                           check.names = FALSE,
                           stringsAsFactors = FALSE)
        names(x) <- train_names
      } else {
        stop(
          "A vector [IVs.test] must contain one complete test row, unless the training data have one predictor.",
          call. = FALSE
        )
      }
    } else {
      x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
    if (ncol(x) != p) {
      stop("[IVs.test] must have the same number of predictors as [IVs.train].",
           call. = FALSE)
    }
    
    if (had_column_names) {
      # Normalize duplicate names with make.unique() identically to the
      # training frame (and to NNS.reg), so cbind(x, x) test input aligns
      # with the c("x", "x.1") training columns rather than erroring.
      supplied <- make.unique(names(x), sep = ".")
      ordering <- .nns_match_predictor_names(supplied, train_names, "IVs.test")
      if (!is.null(ordering)) x <- x[, ordering, drop = FALSE]
    }
    names(x) <- train_names

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
        lev <- if (is.factor(tr))
          levels(droplevels(tr))
        else
          sort(unique(tr_chr))
        unseen <- setdiff(unique(te_chr), lev)
        if (length(unseen)) {
          stop(
            sprintf(
              "Test predictor [%s] contains unseen level(s): %s.",
              nm,
              paste(unseen, collapse = ", ")
            ),
            call. = FALSE
          )
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
    if (anyNA(x))
      stop(sprintf("[%s] contains missing values.", name), call. = FALSE)
    for (j in seq_along(x)) {
      if (is.numeric(x[[j]]) && any(!is.finite(x[[j]]))) {
        stop(
          sprintf(
            "[%s] predictor [%s] contains non-finite values.",
            name,
            names(x)[j]
          ),
          call. = FALSE
        )
      }
    }
  }
  
  .restore_rng <- local({
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old <- if (existed)
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    else
      NULL
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
  
  if (is.null(obj.fn) ||
      !(is.expression(obj.fn) || is.call(obj.fn))) {
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
        stop("[depth] must be NULL, a positive integer, or 'max'.",
             call. = FALSE)
      }
      depth <- "max"
    } else {
      depth <- .scalar_integer(depth, "depth", minimum = 1L)
    }
  }
  
  learner.trials <- .scalar_integer(learner.trials, "learner.trials", minimum = 1L)
  epochs <- .scalar_integer(epochs, "epochs", minimum = 0L, allow_null = TRUE)
  folds <- .scalar_integer(folds, "folds", minimum = 1L)
  
  if (!is.null(CV.size)) {
    if (!is.numeric(CV.size) ||
        length(CV.size) != 1L || !is.finite(CV.size) ||
        CV.size <= 0 || CV.size >= 1) {
      stop("[CV.size] must be a finite scalar strictly between 0 and 1.",
           call. = FALSE)
    }
    CV.size <- as.numeric(CV.size)
  }
  
  ts.test <- .scalar_integer(ts.test,
                             "ts.test",
                             minimum = 1L,
                             allow_null = TRUE)
  
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) ||
        length(threshold) != 1L || !is.finite(threshold)) {
      stop("[threshold] must be a finite numeric scalar or NULL.",
           call. = FALSE)
    }
    threshold <- as.numeric(threshold)
  }
  
  if (!is.null(pred.int)) {
    if (!is.numeric(pred.int) ||
        length(pred.int) != 1L || !is.finite(pred.int) ||
        pred.int <= 0 || pred.int >= 1) {
      stop("[pred.int] must be a finite scalar strictly between 0 and 1.",
           call. = FALSE)
    }
    pred.int <- as.numeric(pred.int)
  }
  
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
    stop("NNS.boost requires at least four training observations.",
         call. = FALSE)
  }
  if (anyNA(DV.train))
    stop("[DV.train] contains missing values.", call. = FALSE)
  if (is.numeric(DV.train) && any(!is.finite(DV.train))) {
    stop("[DV.train] contains non-finite values.", call. = FALSE)
  }
  
  response_was_numeric <- is.numeric(DV.train) ||
    is.integer(DV.train)
  auto_class <- is.factor(DV.train) ||
    is.character(DV.train) || is.logical(DV.train)
  if (auto_class && is.null(type))
    type <- "class"
  is_class <- identical(type, "class")
  
  original_response <- DV.train
  class_values <- NULL
  
  if (is_class) {
    if (response_was_numeric) {
      class_values <- sort(unique(as.numeric(DV.train)))
      y <- match(as.numeric(DV.train), class_values)
    } else {
      class_factor <- if (is.factor(DV.train))
        droplevels(DV.train)
      else
        factor(DV.train)
      class_values <- levels(class_factor)
      y <- as.integer(class_factor)
    }
    y <- as.numeric(y)
    if (length(unique(y)) < 2L) {
      stop("Classification requires at least two response classes.",
           call. = FALSE)
    }
    
    if (identical(obj.fn, expression(sum((
      predicted - actual
    )^2)))) {
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
  
  if (is.null(epochs))
    epochs <- as.integer(2L * n_obs)
  cv_fraction <- if (is.null(CV.size))
    stats::runif(1L, 0.2, 1 / 3)
  else
    CV.size
  
  # ---------------------------------------------------------------------------
  
  # Scoring, splitting, balancing, and prediction helpers
  
  # ---------------------------------------------------------------------------
  
  
  .score <- function(predicted, actual) {
    if (length(predicted) != length(actual)) {
      stop(
        "The objective received predicted and actual vectors of different lengths.",
        call. = FALSE
      )
    }
    value <- eval(
      obj.fn,
      envir = list(predicted = predicted, actual = actual),
      enclos = parent.frame()
    )
    if (!is.numeric(value) || length(value) != 1L) {
      stop("[obj.fn] must return one numeric value.", call. = FALSE)
    }
    value <- as.numeric(value)
    if (!is.finite(value))
      return(NA_real_)
    value
  }
  
  .central_value <- function(v, classification = FALSE) {
    finite <- v[is.finite(v)]
    if (!length(finite))
      return(NA_real_)
    out <- if (classification)
      gravity_class(finite)
    else
      gravity(finite)
    if (!is.finite(out))
      out <- if (classification)
        mode_class(finite)
    else
      mean(finite)
    as.numeric(out)
  }
  
  .sanitize_predictions <- function(predicted, fallback_y) {
    predicted <- as.numeric(predicted)
    bad <- !is.finite(predicted)
    if (any(bad)) {
      replacement <- .central_value(predicted[!bad], is_class)
      if (!is.finite(replacement))
        replacement <- .central_value(fallback_y, is_class)
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
    size <- max(1L, min(n_obs - 1L, as.integer(round(
      cv_fraction * n_obs
    ))))
    for (attempt in seq_len(200L)) {
      idx <- sort(sample.int(n_obs, size = size, replace = FALSE))
      if (.has_all_classes(y[-idx]))
        return(idx)
    }
    stop(
      "Unable to create a validation split retaining every class in training. Reduce CV.size or provide more observations per class.",
      call. = FALSE
    )
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
    if (!balance)
      return(list(x = train_x, y = train_y))
    
    groups <- split(seq_along(train_y), train_y)
    if (length(groups) < 2L || any(lengths(groups) == 0L)) {
      stop("Balancing requires at least two non-empty classes in the fitting split.",
           call. = FALSE)
    }
    
    smallest <- min(lengths(groups))
    largest <- max(lengths(groups))
    
    down_idx <- unlist(lapply(groups, function(g)
      sample(g, smallest, replace = FALSE)), use.names = FALSE)
    up_idx <- unlist(lapply(groups, function(g)
      sample(g, largest, replace = TRUE)), use.names = FALSE)
    
    idx <- sort(c(down_idx, up_idx))
    list(x = train_x[idx, , drop = FALSE], y = train_y[idx])
  }
  
  .fit_subset <- function(feature_index, train_index, test_index) {
    train_x <- x[train_index, feature_index, drop = FALSE]
    train_y <- y[train_index]
    test_x <- x[test_index, feature_index, drop = FALSE]
    
    balanced <- .balance_training(train_x, train_y)
    
    # Every learner trial must remain a genuine NNS.reg base learner on the
    # sampled feature subset. Do not collapse multivariate subsets to an
    # equal-weight synthetic X* before scoring them.
    fit <- suppressWarnings(
      NNS.reg(
        balanced$x,
        balanced$y,
        point.est = test_x,
        dim.red.method = NULL,
        plot = FALSE,
        residual.plot = FALSE,
        order = depth,
        ncores = 1,
        type = type,
        point.only = TRUE,
        dist = dist
      )
    )
    
    .sanitize_predictions(fit$Point.est, balanced$y)
  }
  
  # ---------------------------------------------------------------------------
  
  # Generate learner feature subsets
  
  # ---------------------------------------------------------------------------
  
  
  total_sets <- if (n_features <= 30L)
    2^n_features - 1
  else
    Inf
  exhaustive <- is.finite(total_sets) &&
    total_sets <= learner.trials
  
  if (exhaustive) {
    test.features <- unlist(lapply(seq_len(n_features), function(k) {
      as.list(as.data.frame(utils::combn(n_features, k)))
    }), recursive = FALSE)
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
      warning("Fewer unique feature subsets were generated than requested.",
              call. = FALSE)
    }
  }
  
  learner_count <- length(test.features)
  learner.results <- rep(NA_real_, learner_count)
  
  for (i in seq_len(learner_count)) {
    if (status) {
      message("Current Threshold Iterations Remaining = ",
              learner_count - i,
              " ",
              "\r",
              appendLF = FALSE)
    }
    
    # Cross-sectional learner trials draw a fresh validation holdout per trial
    # (the same resampling geometry as the epoch stage). Time-series trials
    # keep the chronological terminal block.
    trial_validation_index <- if (is.null(ts.test)) {
      .random_validation_index()
    } else {
      validation_index
    }
    trial_train_index <- setdiff(seq_len(n_obs), trial_validation_index)
    
    predicted <- .fit_subset(test.features[[i]],
                             trial_train_index,
                             trial_validation_index)
    learner.results[i] <- .score(predicted, y[trial_validation_index])
  }
  
  finite_results <- which(is.finite(learner.results))
  if (!length(finite_results)) {
    stop("No learner trial produced a finite objective value.", call. = FALSE)
  }
  
  # The public [threshold] argument is a probability supplied to LPM.VaR over
  # the learner-trial objective distribution; the objective-score cutoff is the
  # distinct value LPM.VaR returns. Neither variable overwrites the other.
  threshold_info <- .nns_boost_threshold(
    threshold = threshold,
    objective = objective,
    extreme = extreme,
    learner.results = learner.results[finite_results]
  )
  threshold.probability <- threshold_info$probability
  learner.threshold <- threshold_info$cutoff
  
  if (status) {
    message(
      sprintf(
        "\nLearner threshold probability = %.2f; objective cutoff = %.6f",
        threshold.probability,
        learner.threshold
      )
    )
  }
  
  passes_learner <- is.finite(learner.results) & if (objective == "max") {
    learner.results >= learner.threshold
  } else {
    learner.results <= learner.threshold
  }
  
  reduced.test.features <- test.features[passes_learner]
  
  if (!length(reduced.test.features)) {
    # Defensive: LPM.VaR returns a value inside the observed range, so the
    # best trial always passes; guard anyway for degenerate distributions.
    best_index <- if (objective == "min") {
      finite_results[which.min(learner.results[finite_results])]
    } else {
      finite_results[which.max(learner.results[finite_results])]
    }
    reduced.test.features <- list(test.features[[best_index]])
  }
  
  # ---------------------------------------------------------------------------
  
  # Epoch stability stage
  
  # The learner trials establish the objective threshold and identify the
  # survivor feature combinations. Epochs then repeatedly re-test only those
  # survivor combinations. For ordinary cross-sectional data, each epoch uses
  # a fresh holdout so the final feature frequencies measure out-of-sample
  # stability rather than repeated scoring on one fixed validation sample.
  #
  # Exhaustive learner trials do not disable this stage: complete enumeration
  # answers which combinations passed once, while epochs answer which of those
  # passing combinations continue to pass under repeated holdouts.
  
  # ---------------------------------------------------------------------------
  
  
  keeper.features <- list()
  
  if (epochs > 0L) {
    keeper.features <- vector("list", epochs)
    survivor_count <- length(reduced.test.features)
    
    # Give every survivor approximately equal re-test exposure. Randomizing the
    # balanced schedule avoids ordering effects without allowing exposure count
    # alone to masquerade as feature importance.
    epoch_survivor_id <- rep(seq_len(survivor_count), length.out = epochs)
    if (length(epoch_survivor_id) > 1L) {
      epoch_survivor_id <- sample(epoch_survivor_id,
                                  length(epoch_survivor_id),
                                  replace = FALSE)
    }
    
    # For time-series data, construct expanding-window chronological holdouts of
    # the same size as ts.test. This makes repeated epochs genuine stability
    # checks rather than repetitions of the same terminal block.
    chronological_splits <- NULL
    epoch_split_id <- NULL
    if (!is.null(ts.test)) {
      possible_blocks <- floor((n_obs - 1L) / ts.test)
      block_starts <- n_obs - seq_len(possible_blocks) * ts.test + 1L
      chronological_splits <- lapply(block_starts, function(start) {
        validation <- seq.int(start, start + ts.test - 1L)
        training <- seq_len(start - 1L)
        list(train = training, validation = validation)
      })
      valid_split <- vapply(chronological_splits, function(split) {
        length(split$train) >= 3L && .has_all_classes(y[split$train])
      }, logical(1L))
      chronological_splits <- chronological_splits[valid_split]
      if (!length(chronological_splits)) {
        stop("No chronological epoch split retained enough training observations and every response class.",
             call. = FALSE)
      }
      epoch_split_id <- rep(seq_along(chronological_splits),
                            length.out = epochs)
      if (length(epoch_split_id) > 1L) {
        epoch_split_id <- sample(epoch_split_id,
                                 length(epoch_split_id),
                                 replace = FALSE)
      }
    }
    
    for (j in seq_len(epochs)) {
      if (status) {
        message("% of epochs = ",
                format(j / epochs, digits = 3, nsmall = 2),
                "     ",
                "\r",
                appendLF = FALSE)
      }
      
      features_j <- reduced.test.features[[epoch_survivor_id[j]]]
      
      # Cross-sectional epochs draw a fresh random holdout. Time-series epochs
      # cycle through expanding-window chronological holdouts.
      if (is.null(ts.test)) {
        epoch_validation_index <- .random_validation_index()
        epoch_train_index <- setdiff(seq_len(n_obs), epoch_validation_index)
      } else {
        epoch_split <- chronological_splits[[epoch_split_id[j]]]
        epoch_train_index <- epoch_split$train
        epoch_validation_index <- epoch_split$validation
      }
      epoch_actual <- y[epoch_validation_index]
      
      predicted <- .fit_subset(features_j,
                               epoch_train_index,
                               epoch_validation_index)
      new_result <- .score(predicted, epoch_actual)
      
      passes <- is.finite(new_result) && if (objective == "max") {
        new_result >= learner.threshold
      } else {
        new_result <= learner.threshold
      }
      
      keeper.features[[j]] <- if (passes)
        features_j
      else
        NULL
    }
    
    keeper.features <- keeper.features[!vapply(keeper.features, is.null, logical(1L))]
  } else {
    keeper.features <- reduced.test.features
  }
  
  if (!length(keeper.features)) {
    warning(
      "No feature combination re-passed the objective cutoff during epochs; using the best learner-trial combination. Consider a lower [threshold] probability.",
      call. = FALSE
    )
    best_index <- if (objective == "min") {
      finite_results[which.min(learner.results[finite_results])]
    } else {
      finite_results[which.max(learner.results[finite_results])]
    }
    keeper.features <- list(test.features[[best_index]])
  }
  
  plot.table <- table(factor(unlist(keeper.features), levels = seq_len(n_features)))
  names(plot.table) <- names(x)
  plot.table <- sort(plot.table[plot.table > 0L], decreasing = TRUE)
  
  # Both diagnostic panels are drawn once the final feature frequencies exist,
  # so they are complete before an early features.only return and cannot be
  # disturbed by the final NNS.stack call (whose NNS.reg fits never plot).
  if (feature.importance) {
    .nns_boost_plot_diagnostics(
      learner.results = learner.results[finite_results],
      threshold.probability = threshold.probability,
      learner.threshold = learner.threshold,
      feature.frequency = plot.table
    )
  }
  
  if (features.only) {
    return(.NNS.out(
      list(
        feature.weights = plot.table / sum(plot.table),
        feature.frequency = plot.table
      )
    ))
  }
  
  if (status)
    message("\nGenerating Final Estimate", "\r", appendLF = TRUE)
  
  # ---------------------------------------------------------------------------
  
  # Final estimate: replicate the original keeper predictors by their relative
  # epoch frequencies and fit a genuine multivariate Method 1 NNS.stack.
  #
  # The historical scaling rule converts positive keeper-feature counts into
  # integer replication factors, so a more stable feature contributes
  # proportionally more columns to the final design. No synthetic scalar X* is
  # constructed and the frequencies are never passed through dim.red.method:
  # method = 1 with stack = FALSE keeps the final estimator a multivariate
  # NNS.reg whose n.best is selected by NNS.stack's cross-validation.
  
  # ---------------------------------------------------------------------------
  
  
  feature.frequency <- plot.table
  relative.frequency <- as.numeric(feature.frequency) / min(as.numeric(feature.frequency))
  replication.count <- pmax(1L, as.integer(round(relative.frequency)))
  
  replicated_names <- rep(names(feature.frequency), times = replication.count)
  replicated.train <- x[, replicated_names, drop = FALSE]
  replicated.test <- z[, replicated_names, drop = FALSE]
  names(replicated.train) <- make.unique(names(replicated.train), sep = ".")
  names(replicated.test) <- names(replicated.train)
  
  final.stack <- suppressWarnings(
    NNS.stack(
      IVs.train = replicated.train,
      DV.train = y,
      IVs.test = replicated.test,
      type = type,
      obj.fn = obj.fn,
      objective = objective,
      optimize.threshold = FALSE,
      dist = dist,
      CV.size = CV.size,
      balance = FALSE,
      ts.test = ts.test,
      folds = folds,
      order = depth,
      method = 1L,
      stack = FALSE,
      pred.int = pred.int,
      status = status,
      ncores = 1,
      seed = seed
    )
  )
  
  results <- final.stack$reg
  pred.int.output <- final.stack$reg.pred.int
  
  estimates_code <- .sanitize_predictions(results, y)
  pred_int_out <- pred.int.output
  
  if (is_class) {
    estimates_code <- pmin(pmax(estimates_code, 1L), length(class_values))
    estimates_code <- as.integer(round(estimates_code))
    
    # Classification results keep the historical numeric coding: integer
    # class codes with a base category of 1 (numeric responses recover their
    # original numeric class values). $class.levels supplies the label for
    # each code, so labels are always recoverable via class.levels[results].
    if (response_was_numeric) {
      estimates <- class_values[estimates_code]
    } else {
      estimates <- estimates_code
    }
    
    if (!is.null(pred_int_out)) {
      pred_int_out <- as.data.frame(pred_int_out)
      pred_int_out[] <- lapply(pred_int_out, function(v) {
        code <- pmin(pmax(as.integer(round(v)), 1L), length(class_values))
        if (response_was_numeric)
          class_values[code]
        else
          code
      })
    }
  } else {
    estimates <- estimates_code
  }
  
  .NNS.out(
    list(
      results = estimates,
      pred.int = pred_int_out,
      feature.weights = plot.table / sum(plot.table),
      feature.frequency = plot.table,
      class.levels = if (is_class) class_values else NULL
    )
  )
}


#' Learner threshold from the trial objective distribution
#'
#' Maps the public \code{threshold} probability (with \code{objective} and
#' \code{extreme} defaults) onto the objective-score cutoff via
#' \link{LPM.VaR}. The probability and the cutoff remain distinct values.
#'
#' @keywords internal
#' @noRd
.nns_boost_threshold <- function(threshold, objective, extreme, learner.results) {
  threshold.probability <- threshold
  
  if (is.null(threshold.probability)) {
    threshold.probability <- if (objective == "max") 0.80 else 0.20
  }
  
  if (extreme) {
    threshold.probability <- if (objective == "max") 1 else 0
  }
  
  learner.threshold <- as.numeric(
    LPM.VaR(
      percentile = threshold.probability,
      degree = 1,
      x = learner.results
    )
  )
  
  list(probability = threshold.probability, cutoff = learner.threshold)
}


#' Panel 1: learner-trial objective distribution with its LPM.VaR cutoff
#'
#' @keywords internal
#' @noRd
.nns_boost_plot_learner_distribution <- function(learner.results,
                                                 threshold.probability,
                                                 learner.threshold) {
  graphics::hist(
    learner.results,
    main = "Distribution of Learner Trials Objective Function",
    xlab = "Objective Function",
    col = "steelblue"
  )
  
  graphics::abline(
    v = learner.threshold,
    col = "red",
    lty = 2,
    lwd = 2
  )
  
  graphics::mtext(
    sprintf(
      "LPM.VaR(p = %.2f) = %.4f",
      threshold.probability,
      learner.threshold
    ),
    side = 3,
    col = "red"
  )
  
  invisible(NULL)
}


#' Panel 2: horizontal feature-frequency bar plot
#'
#' @keywords internal
#' @noRd
.nns_boost_plot_feature_frequency <- function(feature.frequency) {
  sorted <- sort(feature.frequency, decreasing = FALSE)

  # Widen this panel's left margin so horizontal las = 1 labels are not
  # clipped by the default margin. The change is scoped to this helper (and
  # margins are per-panel, so the surrounding mfrow layout is untouched).
  old_mai <- graphics::par("mai")
  on.exit(graphics::par(mai = old_mai), add = TRUE)
  label_margin <- max(graphics::strwidth(names(sorted), units = "inches") + 0.4,
                      na.rm = TRUE)
  graphics::par(mai = c(old_mai[1L],
                        max(old_mai[2L], label_margin),
                        old_mai[3L],
                        old_mai[4L]))

  graphics::barplot(
    sorted,
    horiz = TRUE,
    col = "steelblue",
    main = "Feature Frequency in Final Estimate",
    xlab = "Frequency",
    las = 1
  )

  invisible(NULL)
}


#' Two-panel NNS.boost diagnostic on one graphics device
#'
#' Saves the caller's graphics parameters, draws the learner-trial
#' distribution and the feature-frequency panels under
#' \code{par(mfrow = c(2, 1))}, and restores the original parameters only
#' after both panels are complete.
#'
#' @keywords internal
#' @noRd
.nns_boost_plot_diagnostics <- function(learner.results,
                                        threshold.probability,
                                        learner.threshold,
                                        feature.frequency) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::par(mfrow = c(2, 1))
  
  .nns_boost_plot_learner_distribution(
    learner.results = learner.results,
    threshold.probability = threshold.probability,
    learner.threshold = learner.threshold
  )
  
  .nns_boost_plot_feature_frequency(feature.frequency = feature.frequency)
  
  invisible(NULL)
}
