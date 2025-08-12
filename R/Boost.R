#' NNS Boost
#'
#' Ensemble method for classification using the NNS multivariate regression \link{NNS.reg} as the base learner instead of trees.
#'
#' @param IVs.train a matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimensions to \code{(IVs.train)}.
#' @param IVs.test a matrix or data frame of variables of numeric or factor data types with compatible dimensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.
#' @param type \code{NULL} (default).  To perform a classification of discrete integer classes from factor target variable \code{(DV.train)} with a base category of 1, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.
#' @param depth options: (integer, NULL, "max"); \code{(depth = NULL)}(default) Specifies the \code{order} parameter in the \link{NNS.reg} routine, assigning a number of splits in the regressors, analogous to tree depth.
#' @param learner.trials integer; 100 (default) Sets the number of trials to obtain an accuracy \code{threshold} level.  If the number of all possible feature combinations is less than selected value, the minimum of the two values will be used.
#' @param epochs integer; \code{2*length(DV.train)} (default) Total number of feature combinations to run.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size.  Defaults to a random value between 0.2 and 0.33 for a random sampling of the training set.
#' @param balance logical; \code{FALSE} (default) Uses both up and down sampling to balance the classes.  \code{type="CLASS"} required.
#' @param ts.test integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.
#' @param folds integer; 5 (default) Sets the number of \code{folds} in the \link{NNS.stack} procedure for optimal \code{n.best} parameter.
#' @param threshold numeric; \code{NULL} (default) Sets the \code{obj.fn} threshold to keep feature combinations.
#' @param obj.fn expression;
#' \code{expression( sum((predicted - actual)^2) )} (default) Sum of squared errors is the default objective function.  Any \code{expression(...)} using the specific terms \code{predicted} and \code{actual} can be used.  Automatically selects an accuracy measure when \code{(type = "CLASS")}.
#' @param objective options: ("min", "max") \code{"max"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param extreme logical; \code{FALSE} (default) Uses the maximum (minimum) \code{threshold} obtained from the \code{learner.trials}, rather than the upper (lower) quintile level for maximization (minimization) \code{objective}.
#' @param features.only logical; \code{FALSE} (default) Returns only the final feature loadings along with the final feature frequencies.
#' @param feature.importance logical; \code{TRUE} (default) Plots the frequency of features used in the final estimate.
#' @param pred.int numeric [0,1]; \code{NULL} (default) Returns the associated prediction intervals for the final estimate.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#'
#' @return Returns a vector of fitted values for the dependent variable test set \code{$results}, prediction intervals \code{$pred.int}, and the final feature loadings \code{$feature.weights}, along with final feature frequencies \code{$feature.frequency}.
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
#'  type = "CLASS", depth = NULL)
#'
#'  ## Test accuracy
#'  mean(a$results == as.numeric(iris[141:150, 5]))
#'  }
#'
#' @export


NNS.boost <- function(
    IVs.train,
    DV.train,
    IVs.test = NULL,
    type = NULL,
    depth = NULL,
    learner.trials = 100L,
    epochs = NULL,
    CV.size = NULL,
    balance = FALSE,
    ts.test = NULL,
    folds = 5L,
    threshold = NULL,
    obj.fn = expression(sum((predicted - actual)^2)),
    objective = "min",
    extreme = FALSE,
    features.only = FALSE,
    feature.importance = TRUE,
    pred.int = NULL,
    status = TRUE,
    ncores = NULL,
    seed = 123L
) {
  # --- Basic checks -----------------------------------------------------------
  if (is.null(obj.fn)) stop("Please provide an objective function.")
  if (sum(is.na(cbind(IVs.train, DV.train))) > 0) stop("You have some missing values, please address.")
  if (length(DV.train) != nrow(IVs.train)) stop("DV.train length must match nrow(IVs.train).")
  if (!is.null(IVs.test) && !is.null(colnames(IVs.train)) &&
      !setequal(colnames(IVs.test), colnames(IVs.train)))
    stop("IVs.test must have the same columns as IVs.train (names may be in any order).")
  
  objective <- match.arg(tolower(objective), c("min","max"))
  
  if (balance && is.null(type)) warning("type = 'CLASS' selected due to balance = TRUE.")
  if (balance) type <- "CLASS"
  
  if (!is.null(type)) {
    type <- toupper(type)
    if (type == "CLASS" && identical(obj.fn, expression(sum((predicted - actual)^2)))) {
      obj.fn <- expression(mean(predicted == as.numeric(actual)))
      objective <- "max"
    }
  }
  
  if (!is.null(type) && min(as.numeric(as.factor(DV.train))) == 0)
    warning("Base response variable category should be 1, not 0.")
  
  # --- Coerce storage once ----------------------------------------------------
  if (inherits(IVs.train, c("tbl","data.table"))) IVs.train <- as.data.frame(IVs.train)
  if (inherits(DV.train,  c("tbl","data.table"))) DV.train  <- as.vector(unlist(DV.train))
  
  if (is.null(colnames(IVs.train))) {
    colnames(IVs.train) <- paste0("X", seq_len(ncol(IVs.train)))
    if (!is.null(IVs.test)) colnames(IVs.test) <- colnames(IVs.train)
  }
  
  features <- sort(colnames(IVs.train))
  IVs.train <- IVs.train[, features, drop = FALSE]
  
  if (is.null(IVs.test)) {
    IVs.test <- IVs.train
  } else {
    if (inherits(IVs.test, c("tbl","data.table"))) IVs.test <- as.data.frame(IVs.test)
    IVs.test <- IVs.test[, features, drop = FALSE]
  }
  
  # --- Balance block -----------------------
  if (balance) {
    set.seed(seed)  # your up/downSample use sample(); make it reproducible
    y_fac <- factor(DV.train)
    
    dwn <- downSample(IVs.train, y_fac, list = FALSE, yname = "DV.train")
    up  <- upSample  (IVs.train, y_fac, list = FALSE, yname = "DV.train")
    
    training <- rbind(dwn, up)
    
    IVs.train <- training[, features, drop = FALSE]
    DV.train  <- as.integer(factor(training[["DV.train"]]))
    
    if (anyNA(IVs.train) || anyNA(DV.train)) {
      stop("balance step produced missing values; verify predictors are numeric.")
    }
  }
  
  # --- Data prep --------------------------------------------------------------
  x <- data.table::as.data.table(IVs.train)
  y <- DV.train
  z <- data.table::as.data.table(IVs.test)
  stopifnot(nrow(x) == length(y))   # fail fast if anything drifted
  
  n <- ncol(x)
  N <- nrow(x)
  
  if (is.null(epochs)) epochs <- 2L * length(y)
  dist <- if (!is.null(ts.test)) "DTW" else "L2"
  
  # --- Representative samples -------------------
  # 5-number summary for y (length 5)
  rep.y <- stats::fivenum(as.numeric(y))
  # 5 x p matrix for X, then 5-row data.frame with column names = features
  rep.x <- apply(as.matrix(x), 2, function(v) stats::fivenum(as.numeric(v)))
  rep.x <- as.data.frame(rep.x)
  colnames(rep.x) <- features  # rep.x: 5 rows, p columns
  
  # --- Deterministic mode -----------------------------------------------------
  sets <- sum(choose(n, 1:n))
  deterministic <- ((sets < length(y)) || n <= 10)
  if (deterministic) {
    learner.trials <- sets
    combn_list <- lapply(1:n, function(k) combn(n, k, simplify = FALSE))
    deterministic.sets <- unlist(combn_list, recursive = FALSE)
  } else {
    if (is.null(learner.trials) || learner.trials < 1) learner.trials <- length(y)
  }
  
  pick_indices <- function(i, l = length(y), frac = 0.25, ts_len = NULL) {
    if (!is.null(ts_len)) return(seq.int(l - ts_len + 1L, l))
    k <- max(1L, as.integer((if (is.null(CV.size)) stats::runif(1, .2, 1/3) else CV.size) * l))
    if (i <= l/4) as.integer(seq.int(i, l, length.out = k)) else sample.int(l, k, replace = FALSE)
  }
  
  set.seed(seed)
  if (is.null(threshold)) {
    if (!extreme) epochs <- NULL
    results <- numeric(learner.trials)
    test.features <- vector("list", learner.trials)
    
    for (i in seq_len(learner.trials)) {
      idx <- pick_indices(i, l = length(y), ts_len = if (!is.null(ts.test)) ts.test else NULL)
      keep <- setdiff(seq_len(N), idx)
      
      new_iv_train <- data.table::copy(x[keep])
      new_dv_train <- y[keep]
      
      # anchors: 5-number summaries on the current training subset (5 x p)
      anchors <- apply(as.matrix(new_iv_train), 2, function(v) stats::fivenum(as.numeric(v)))
      anchors <- as.data.frame(anchors)
      colnames(anchors) <- features
      
      # stack anchors (5 rows) + raw rows
      new_iv_train <- data.table::rbindlist(list(anchors, new_iv_train), use.names = TRUE, fill = FALSE)
      new_dv_train <- c(rep.y, new_dv_train)
      
      actual  <- as.numeric(y[idx])
      new_iv_test <- x[idx]
      
      if (status) message("Current Threshold Iterations Remaining = ", learner.trials + 1L - i, " \r", appendLF = FALSE)
      
      feats <- if (deterministic) deterministic.sets[[i]] else sort(sample.int(n, sample.int(n - 1L, 1L) + 1L, replace = FALSE))
      test.features[[i]] <- feats
      
      Xtr <- as.matrix(new_iv_train)[, feats, drop = FALSE]
      Xte <- as.matrix(new_iv_test)[ , feats, drop = FALSE]
      
      predicted <- NNS.reg(
        x = Xtr,
        y = new_dv_train,
        point.est = Xte,
        dim.red.method = "equal",
        plot = FALSE, plot.regions = FALSE, residual.plot = FALSE,
        order = depth, ncores = ncores, type = if (!is.null(type)) type else NULL, dist = dist, smooth = TRUE
      )$Point.est
      
      if (anyNA(predicted)) predicted[is.na(predicted)] <- NNS.gravity(na.omit(predicted))
      if (!is.null(type)) {
        predicted <- pmin(predicted, max(as.numeric(y)))
        predicted <- pmax(predicted, min(as.numeric(y)))
      }
      
      results[i] <- eval(obj.fn)
    }
    
    if (extreme) {
      threshold <- if (objective == "max") max(results) else min(results)
    } else {
      q <- stats::fivenum(results)
      threshold <- if (objective == "max") q[4L] else q[2L]
    }
    
  } else {
    results <- as.numeric(threshold)
  }
  
  # --- Dual plot setup (panel 1 now, panel 2 later) --------------------------
  if (feature.importance) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(2, 1))
    
    hist(results, main = "Distribution of Learner Trials Objective Function",
         xlab = "Objective Function", col = "steelblue")
    abline(v = threshold, col = 'red', lty = 2, lwd = 2)
    mtext(round(threshold, 3), side = 1, col = "red", at = threshold)
    if (extreme) {
      mtext(if (objective == "max") "Threshold >" else "< Threshold", side = 3, col = "red", at = threshold)
    }
  }
  
  # --- Feature filtering ------------------------------------------------------
  pass_idx <- if (extreme) {
    if (objective == "max") which.max(results) else which.min(results)
  } else {
    if (objective == "max") which(results >= threshold) else which(results <= threshold)
  }
  reduced.test.features <- if (length(pass_idx)) test.features[pass_idx] else list()
  
  if (length(reduced.test.features)) {
    freq_vec <- tabulate(unlist(reduced.test.features), nbins = n)
    names(freq_vec) <- features
  } else {
    if (is.null(threshold)) stop("No feature sets meet the threshold; adjust `threshold` or increase trials.")
    best_idx <- if (objective == "max") which.max(results) else which.min(results)
    freq_vec <- tabulate(unlist(test.features[[best_idx]]), nbins = n)
    names(freq_vec) <- features
  }
  
  # --- Epoch refinement -------------------------------------------------------
  keeper.features <- list()
  if (!deterministic && !is.null(epochs) && epochs > 0L) {
    for (j in seq_len(epochs)) {
      idx <- pick_indices(j, l = length(y), ts_len = if (!is.null(ts.test)) (2L * ts.test) else NULL)
      keep <- setdiff(seq_len(N), idx)
      
      new_iv_train <- data.table::copy(x[keep])
      new_dv_train <- y[keep]
      
      anchors <- apply(as.matrix(new_iv_train), 2, function(v) stats::fivenum(as.numeric(v)))
      anchors <- as.data.frame(anchors)
      colnames(anchors) <- features
      
      new_iv_train <- data.table::rbindlist(list(anchors, new_iv_train), use.names = TRUE, fill = FALSE)
      new_dv_train <- c(rep.y, new_dv_train)
      
      actual  <- as.numeric(y[idx])
      new_iv_test <- x[idx]
      
      if (status) {
        message(sprintf("%% of epochs = %.2f  \r", j/epochs), appendLF = FALSE)
        if (j == epochs) { message(sprintf("%% of epochs %d = 1.00  \r", j), appendLF = FALSE); flush.console() }
      }
      
      base_pool <- which(freq_vec > 0L)
      extra <- sample.int(n, sample.int(n, 1L), replace = FALSE)
      feats <- sort(unique(c(base_pool, extra)))
      
      Xtr <- as.matrix(new_iv_train)[, feats, drop = FALSE]
      Xte <- as.matrix(new_iv_test)[ , feats, drop = FALSE]
      
      predicted <- NNS.reg(
        x = Xtr,
        y = new_dv_train,
        point.est = Xte,
        dim.red.method = "equal",
        plot = FALSE, plot.regions = FALSE, residual.plot = FALSE,
        order = depth, ncores = ncores, type = if (!is.null(type)) type else NULL, dist = dist, smooth = TRUE
      )$Point.est
      
      if (anyNA(predicted)) predicted[is.na(predicted)] <- NNS.gravity(na.omit(predicted))
      if (!is.null(type)) {
        predicted <- pmin(predicted, max(as.numeric(y)))
        predicted <- pmax(predicted, min(as.numeric(y)))
      }
      
      score <- eval(obj.fn)
      keep_it <- if (objective == "max") (ifelse(is.na(score), 0.99 * threshold, score) >= threshold)
      else                      (ifelse(is.na(score), 1.01 * threshold, score) <= threshold)
      if (keep_it) keeper.features[[length(keeper.features) + 1L]] <- feats
    }
  } else {
    keeper.features <- reduced.test.features
  }
  
  keeper.features <- keeper.features[lengths(keeper.features) > 0L]
  if (!length(keeper.features)) {
    best_idx <- if (objective == "max") which.max(results) else which.min(results)
    keeper.features <- list(test.features[[best_idx]])
  }
  
  final_counts <- tabulate(unlist(keeper.features), nbins = n)
  names(final_counts) <- features
  final_counts <- final_counts[final_counts > 0L]
  final_counts <- sort(final_counts, decreasing = TRUE)
  
  if (features.only) {
    fw <- final_counts / sum(final_counts)
    return(list(feature.weights = fw, feature.frequency = final_counts))
  }
  
  # Add representative rows before final fit
  if (!is.null(rep.y)) {
    x <- data.table::as.data.table(rbind(rep.x, x))  # rep.x: 5 rows x p
    y <- c(rep.y, y)                                  # rep.y: length 5
  }
  
  if (status) message("Generating Final Estimate\r", appendLF = TRUE)
  
  keep_cols <- intersect(colnames(x), names(final_counts))
  X_train <- as.matrix(x)[, keep_cols, drop = FALSE]
  X_test  <- as.matrix(z)[, keep_cols, drop = FALSE]
  
  model <- NNS.stack(
    IVs.train = X_train,
    DV.train  = y,
    IVs.test  = X_test,
    order = depth,
    dim.red.method = "all",
    ncores = ncores,
    stack = FALSE,
    status = status,
    type = if (!is.null(type)) type else NULL,
    dist = dist,
    folds = folds,
    pred.int = pred.int
  )
  
  estimates <- model$stack
  if (anyNA(estimates)) {
    fill <- NNS.mode(na.omit(estimates))
    estimates[is.na(estimates)] <- fill
  }
  if (!is.null(type)) {
    estimates <- pmin(estimates, max(as.numeric(y)))
    estimates <- pmax(estimates, min(as.numeric(y)))
    estimates <- ifelse(estimates %% 1 < .5, floor(estimates), ceiling(estimates))
  }
  
  # --- Dual plot panel 2 (feature barplot) -----------------------------------
  if (feature.importance) {
    linch <- max(strwidth(names(final_counts), "inch") + 0.4, na.rm = TRUE)
    par(mai = c(1.0, linch, 0.8, 0.5))
    
    k <- min(length(final_counts), 10L)
    barplot(
      sort(final_counts, decreasing = FALSE)[seq_len(k)],
      horiz = TRUE, col = "steelblue",
      main = "Feature Frequency in Final Estimate",
      xlab = "Frequency", las = 1
    )
  }
  
  list(
    results = estimates,
    pred.int = model$pred.int,
    feature.weights = final_counts / sum(final_counts),
    feature.frequency = final_counts
  )
}
