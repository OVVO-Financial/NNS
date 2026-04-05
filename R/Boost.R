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
#'  type = "CLASS", depth = NULL, balance = TRUE)
#'
#'  ## Test accuracy
#'  mean(a$results == as.numeric(iris[141:150, 5]))
#'  }
#'
#' @export



# NNS Boost (balanced + robust)
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
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = "min",
                      extreme = FALSE,
                      features.only = FALSE,
                      feature.importance = TRUE,
                      pred.int = NULL,
                      status = TRUE){
  
  .core <- function() {
    
    if (anyNA(cbind(IVs.train, DV.train))) stop("You have some missing values, please address.")
    if (is.null(obj.fn)) stop("Please provide an objective function")
    
    if (balance && is.null(type)) warning("type = 'CLASS' selected due to balance = TRUE.")
    if (balance) type <- "CLASS"
    
    if (!is.null(type) && min(as.numeric(as.factor(DV.train))) == 0)
      warning("Base response variable category should be 1, not 0.")
    
    if (any(class(IVs.train) %in% c("tbl","data.table"))) IVs.train <- as.data.frame(IVs.train)
    if (any(class(DV.train) %in% c("tbl","data.table"))) DV.train <- as.vector(unlist(DV.train))
    
    if (!is.null(type)) {
      type <- tolower(type)
      if (type == "class" && identical(obj.fn, expression( sum((predicted - actual)^2) ))) {
        obj.fn <- expression(mean(predicted == as.numeric(actual)))
        objective <- "max"
      }
    }
    
    objective <- tolower(objective)
    
    if (is.null(colnames(IVs.train))) {
      colnames(IVs.train) <- paste0("X", seq_len(ncol(IVs.train)))
    }
    features <- colnames(IVs.train)
    
    IVs.train <- IVs.train[, sort(features), drop = FALSE]
    transform  <- data.matrix(cbind(DV.train, IVs.train))
    IVs.train  <- transform[, -1, drop = FALSE]
    colnames(IVs.train) <- sort(features)
    DV.train   <- transform[,  1]
    
    if (is.null(IVs.test)) {
      IVs.test <- IVs.train
    } else {
      if (any(class(IVs.test) %in% c("tbl","data.table"))) IVs.test <- as.data.frame(IVs.test)
      colnames(IVs.test) <- colnames(IVs.train)
    }
    
    if (balance) {
      y_train  <- as.factor(DV.train)
      ycol <- "Class"
      training_1 <- downSample(IVs.train, y_train, list = FALSE, yname = ycol)
      training_2 <- upSample(IVs.train,   y_train, list = FALSE, yname = ycol)
      training <- rbind.data.frame(training_1, training_2)
      IVs.train <- training[, setdiff(names(training), ycol), drop = FALSE]
      DV.train  <- as.numeric(as.factor(training[[ycol]]))
      colnames(IVs.test) <- colnames(IVs.train)  
    }
    
    x <- data.table::data.table(IVs.train)
    y <- DV.train
    z <- data.table::data.table(IVs.test)
    
    
    n <- ncol(x)
    if (is.null(epochs)) epochs <- 2*length(y)
    dist <- if (!is.null(ts.test)) "DTW" else "L2"
    
    old.threshold <- 0
    
    sets <- sum(choose(n, 1:n))
    deterministic <- FALSE
    if ((sets < length(y)) || n <= 10) {
      deterministic <- TRUE
      learner.trials <- sets
      combn_vec <- Vectorize(Rfast::comb_n, vectorize.args = "k")
      deterministic.sets <- unlist(lapply(combn_vec(n, 1:n), function(df) as.list(as.data.frame(df))), recursive = FALSE)
    }
    
    if (is.null(threshold)) {
      new.CV.size <- if (is.null(CV.size)) round(runif(1, .2, 1/3), 3) else CV.size
      old.threshold <- 1
      if (is.null(learner.trials)) learner.trials <- length(y)
      
      results <- numeric(learner.trials)
      test.features <- vector(mode = "list", learner.trials)
      
      for (i in 1:learner.trials) {
        set.seed(123 + i)
        l <- length(y)
        if (i <= l/4) new.index <- as.integer(seq(i, length(y), length.out = as.integer(new.CV.size * length(y))))
        else          new.index <- sample(l, as.integer(new.CV.size * l), replace = FALSE)
        if (!is.null(ts.test)) new.index <- 1:(length(y) - ts.test)
        new.index <- unlist(new.index)
        
        new.iv.train <- cbind(y[-new.index], x[-new.index,])
        new.iv.train <- new.iv.train[, lapply(.SD, as.double)]
        new.iv.train <- new.iv.train[, lapply(.SD, function(z) fivenum(as.numeric(z)))]
        new.dv.train <- unlist(new.iv.train[,1])
        new.iv.train <- as.data.frame(new.iv.train)
        new.iv.train <- new.iv.train[, unlist(colnames(new.iv.train) %in% colnames(IVs.train)), drop = FALSE]
        new.iv.train <- data.table::rbindlist(list(new.iv.train, x[-new.index,]), use.names = FALSE)
        new.dv.train <- c(new.dv.train, y[-new.index])
        
        colnames(new.iv.train) <- c(colnames(IVs.train))
        
        actual <- as.numeric(y[new.index])
        new.iv.test <- x[new.index,]
        
        if (status) message("Current Threshold Iterations Remaining = ", learner.trials + 1 - i, " ", "\r", appendLF = FALSE)
        
        if (deterministic) test.features[[i]] <- deterministic.sets[[i]] else test.features[[i]] <- sort(sample(n, sample(2:n, 1), replace = FALSE))
        
        learning.IVs <- as.data.frame(new.iv.train)[ , as.integer(unlist(test.features[[i]])), drop = FALSE]
        point.IVs    <- as.data.frame(new.iv.test)[  , as.integer(unlist(test.features[[i]])), drop = FALSE]
        
        predicted <- NNS.reg(learning.IVs,
                             new.dv.train,
                             point.est = point.IVs,
                             dim.red.method = "equal",
                             plot = FALSE, order = depth,
                             ncores = 1, type = type)$Point.est
        
        predicted[is.na(predicted)] <- gravity(na.omit(predicted))
        
        if (!is.null(type)) {
          predicted <- pmin(predicted, max(as.numeric(y)))
          predicted <- pmax(predicted, min(as.numeric(y)))
        }
        
        results[i] <- eval(obj.fn)
      }
    } else {
      results <- threshold
    }
    
    if (extreme) {
      threshold <- if (objective == "max") max(results) else min(results)
    } else {
      threshold <- if (objective == "max") fivenum(results)[4] else fivenum(results)[2]
    }
    
    if (status) {
      message(paste0("\nLearner Accuracy Threshold = ", format(threshold, digits = 3, nsmall = 2), "           "), appendLF = TRUE)
    }
    
    if (extreme) {
      reduced.test.features <- if (objective == "max") test.features[which.max(results)] else test.features[which.min(results)]
    } else {
      reduced.test.features <- if (objective == "max") test.features[which(results >= threshold)] else test.features[which(results <= threshold)]
    }
    
    # Build a weighted feature sampling pool from the surviving learner-trial sets.
    # scale_factor_rf gives each feature index a count proportional to how often
    # it appeared across the surviving sets; feature.pool is a flat index vector
    # used for weighted random sampling inside the epoch loop.
    rf <- data.table::data.table(table(as.character(reduced.test.features)))
    rf$N <- rf$N / sum(rf$N)
    rf_reduced <- apply(rf, 1, function(x) eval(parse(text = x[1])))
    scale_factor_rf <- table(unlist(rf_reduced)) / min(table(unlist(rf_reduced)))
    feature.pool <- as.numeric(rep(names(scale_factor_rf),
                                   ifelse(scale_factor_rf %% 1 < .5,
                                          floor(scale_factor_rf),
                                          ceiling(scale_factor_rf))))
    # reduced.test.features remains a list for set-level operations downstream
    
    keeper.features <- list()
    if (deterministic) epochs <- NULL
    
    if (!is.null(epochs) && !deterministic) {
      new.CV.size <- if (is.null(CV.size)) round(runif(1, .2, 1/3), 3) else CV.size
      for (j in 1:epochs) {
        set.seed(123 * j)
        l <- length(y)
        if (j <= l/4) new.index <- as.integer(seq(j, length(y), length.out = as.integer(new.CV.size * length(y))))
        else          new.index <- sample(l, as.integer(new.CV.size * l), replace = FALSE)
        if (!is.null(ts.test)) new.index <- length(y) - (2 * ts.test):0
        new.index <- unlist(new.index)
        
        new.iv.train <- cbind(y[-new.index], x[-new.index, ])
        new.iv.train <- new.iv.train[, lapply(.SD, as.double)]
        new.iv.train <- new.iv.train[, lapply(.SD, function(z) fivenum(as.numeric(z)))]
        new.dv.train <- unlist(new.iv.train[, 1])
        new.iv.train <- as.data.frame(new.iv.train)
        new.iv.train <- new.iv.train[, unlist(colnames(new.iv.train) %in% colnames(IVs.train)), drop = FALSE]
        new.iv.train <- data.table::rbindlist(list(new.iv.train, x[-new.index, ]), use.names = FALSE)
        new.dv.train <- c(new.dv.train, y[-new.index])
        colnames(new.iv.train) <- colnames(IVs.train)
        
        actual <- as.numeric(y[new.index])
        new.iv.test <- x[new.index, ]
        
        if (status) message("% of epochs = ", format(j / epochs, digits = 3, nsmall = 2), "     ", "\r", appendLF = FALSE)
        
        # Each epoch: draw a random number of features (1..n) from the weighted
        # pool so that both the COUNT and the COMBINATION vary across epochs.
        # feature.pool has high-frequency features repeated more often, so
        # sample(..., replace = FALSE) naturally over-represents them while
        # still allowing any combination of size 1..n to appear.
        # unique() + sort() prevents duplicate columns (corrupts L2 distance).
        features_j <- if (deterministic) {
          unlist(deterministic.sets[[j]])
        } else {
          k <- sample(seq_len(n), 1L)                          # random feature count 1..n
          sort(unique(sample(feature.pool, k, replace = TRUE))) # weighted draw, unique indices
        }
        
        learning_IVs_epoch <- as.data.frame(new.iv.train)[ , as.integer(features_j), drop = FALSE]
        point_est_epoch    <- as.data.frame(new.iv.test)[  , as.integer(features_j), drop = FALSE]
        
        predicted <- NNS.reg(learning_IVs_epoch,
                             new.dv.train,
                             point.est = point_est_epoch,
                             dim.red.method = "equal",
                             plot = FALSE, residual.plot = FALSE, order = depth,
                             ncores = 1, type = type)$Point.est
        
        predicted[is.na(predicted)] <- gravity(na.omit(predicted))
        if (!is.null(type)) {
          predicted <- pmin(pmax(predicted, min(as.numeric(y))), max(as.numeric(y)))
        }
        
        new.results <- eval(obj.fn)
        passes <- if (objective == "max") {
          if (is.na(new.results)) new.results <- .99 * threshold
          new.results >= threshold
        } else {
          if (is.na(new.results)) new.results <- 1.01 * threshold
          new.results <= threshold
        }
        keeper.features[[j]] <- if (passes) features_j else NULL
      }
    } else {
      keeper.features <- reduced.test.features
    }
    
    keeper.features <- keeper.features[!sapply(keeper.features, is.null)]
    
    # Fallback: if no epoch passed the threshold, use the single best learner-trial
    if (length(keeper.features) == 0) {
      if (old.threshold == 0) {
        if (objective == "min") stop("Please increase [threshold].") else stop("Please reduce [threshold].")
      }
      best_feat <- if (objective == "min") test.features[[which.min(results)]] else test.features[[which.max(results)]]
      keeper.features <- list(best_feat)
    }
    
    plot.table <- table(unlist(keeper.features))
    names(plot.table) <- colnames(IVs.train)[as.numeric(names(plot.table))]
    if (features.only || feature.importance) plot.table <- plot.table[rev(order(plot.table))]
    if (features.only) {
      return(list("feature.weights"   = plot.table / sum(plot.table),
                  "feature.frequency" = plot.table))
    }
    
    if (status) message("\nGenerating Final Estimate", "\r", appendLF = TRUE)
    
    # Build a frequency-weighted synthetic predictor X* from all features, where
    # each column is weighted by how often it survived the threshold filter.
    # NNS.reg with dim.red.method = coef_aligned computes X* internally and
    # exposes it via $x.star (training rows). $Point.est is the regression
    # output (predicted Y), NOT the X* projection of test rows, so we compute
    # the test-set X* explicitly using the same joint-normalisation that
    # NNS.reg applies internally:
    #   norm.x  <- apply(rbind(test, train), 2, rescale)
    #   X*      <- norm.x %*% coef / sum(abs(coef) > 0)
    # X* is then duplicated into cbind(xstar, xstar) so that NNS.stack
    # method = 1 can cross-validate n.best on a two-column design matrix.
    # The duplicate column satisfies the multivariate path requirement
    # without adding new information. The same obj.fn and objective carried
    # through the boost loop govern n.best selection.
    freq_weights <- as.numeric(plot.table / sum(plot.table))  # normalised frequencies
    names(freq_weights) <- names(plot.table)
    # align to column order of IVs.train (x)
    coef_aligned <- freq_weights[colnames(x)]
    coef_aligned[is.na(coef_aligned)] <- 0
    
    xstar_fit <- suppressWarnings(
      NNS.reg(as.data.frame(x), y,
              dim.red.method = coef_aligned,
              plot           = FALSE,
              residual.plot  = FALSE,
              order          = depth,
              ncores         = 1,
              type           = NULL,
              point.only     = FALSE)
    )
    
    xstar_train <- as.numeric(unlist(xstar_fit$x.star))
    xstar_train[is.na(xstar_train)] <- gravity(na.omit(xstar_train))
    
    # Replicate NNS.reg joint-normalisation to project test rows onto X*
    x_mat  <- data.matrix(as.data.frame(x))
    z_mat  <- data.matrix(as.data.frame(z))
    joint  <- rbind(z_mat, x_mat)
    joint_norm <- apply(joint, 2, function(col) {
      rng <- max(col) - min(col)
      (col - min(col)) / ifelse(rng == 0, 1, rng)
    })
    xn <- sum(abs(coef_aligned) > 0)
    if (xn == 0) xn <- 1L
    xstar_test <- as.numeric(
      joint_norm[seq_len(nrow(z_mat)), , drop = FALSE] %*% coef_aligned / xn
    )
    xstar_test[is.na(xstar_test)] <- gravity(na.omit(xstar_test))
    
    IVs.xstar.train <- data.frame(xstar = xstar_train, xstar2 = xstar_train)
    IVs.xstar.test  <- data.frame(xstar = xstar_test,  xstar2 = xstar_test)
    
    final_fit <- suppressWarnings(
      NNS.stack(IVs.train  = IVs.xstar.train,
                DV.train   = y,
                IVs.test   = IVs.xstar.test,
                method     = 1,
                obj.fn     = obj.fn,
                objective  = objective,
                type       = type,
                pred.int   = pred.int,
                status     = status)
    )
    
    estimates <- final_fit$stack
    if (is.null(estimates)) estimates <- final_fit$reg
    estimates[is.na(estimates)] <- gravity(na.omit(estimates))
    
    if (!is.null(type)) {
      estimates <- pmin(pmax(estimates, min(as.numeric(y))), max(as.numeric(y)))
      estimates <- ifelse(estimates %% 1 < .5, floor(estimates), ceiling(estimates))
    }
    
    if (feature.importance) {
      linch <- max(strwidth(names(plot.table), "inch") + 0.4, na.rm = TRUE)
      par(mai = c(1.0, linch, 0.8, 0.5))
      if (length(plot.table) != 1) {
        barplot(sort(plot.table, decreasing = FALSE)[1:min(n, 10)], horiz = TRUE,
                col = 'steelblue', main = "Feature Frequency in Final Estimate",
                xlab = "Frequency", las = 1)
      } else {
        barplot(sort(plot.table, decreasing = FALSE), horiz = TRUE,
                col = 'steelblue', main = "Feature Frequency in Final Estimate",
                xlab = "Frequency", las = 1)
      }
      par(mfrow = c(1,1))
    }
    
    return(list("results"           = estimates,
                "pred.int"          = final_fit$pred.int,
                "feature.weights"   = plot.table / sum(plot.table),
                "feature.frequency" = plot.table,
                "n.best"            = final_fit$NNS.reg.n.best))
  } # end .core
  
  out <- tryCatch(
    suppressWarnings(.core()),
    error = function(e) {
      if (isTRUE(balance)) {
        warning("[retry] First attempt failed; retrying with balance = FALSE")
        return(NNS.boost(IVs.train = IVs.train, DV.train = DV.train, IVs.test = IVs.test,
                         type = type, depth = depth, learner.trials = learner.trials,
                         epochs = epochs, CV.size = CV.size, balance = FALSE, ts.test = ts.test,
                         threshold = threshold, obj.fn = obj.fn,
                         objective = objective, extreme = extreme, features.only = features.only,
                         feature.importance = feature.importance, pred.int = pred.int, status = status))
      }
      stop(e)
    }
  )
  
  out
}