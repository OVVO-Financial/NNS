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
                      folds = 5,
                      threshold = NULL,
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = "min",
                      extreme = FALSE,
                      features.only = FALSE,
                      feature.importance = TRUE,
                      pred.int = NULL,
                      status = TRUE){
  
  # Check for missing values
  if(sum(is.na(cbind(IVs.train, DV.train)) > 0)) 
     stop("You have some missing values, please address.")
     
     if(is.null(obj.fn)) 
       stop("Please provide an objective function")
     
     # Handle balancing
     if(balance && is.null(type)) {
       warning("type = 'CLASS' selected due to balance = TRUE.")
       type <- "CLASS"
     }
     
     if(!is.null(type) && min(as.numeric(as.factor(DV.train))) == 0) 
       warning("Base response variable category should be 1, not 0.")
     
     # Convert to data frames
     if(any(class(IVs.train) %in% c("tbl", "data.table"))) 
       IVs.train <- as.data.frame(IVs.train)
     
     if(any(class(DV.train) %in% c("tbl", "data.table"))) 
       DV.train <- as.vector(unlist(DV.train))
     
     # Set classification objective function
     if(!is.null(type)){
       type <- tolower(type)
       if(type == "class" && identical(obj.fn, expression( sum((predicted - actual)^2) ))){
         obj.fn <- expression(mean( predicted == as.numeric(actual)))
         objective <- "max"
       }
     }
     
     objective <- tolower(objective)
     
     # Handle column names
     if(is.null(colnames(IVs.train))){
       colnames(IVs.train) <- paste0("X", 1:ncol(IVs.train))
     }
     
     if(!is.null(IVs.test) && is.null(colnames(IVs.test))) {
       colnames(IVs.test) <- colnames(IVs.train)
     }
     
     features <- colnames(IVs.train)
     IVs.train <- IVs.train[, sort(features), drop = FALSE]
     
     # Handle test set
     if(is.null(IVs.test)){
       IVs.test <- IVs.train
     } else {
       if(any(class(IVs.test) %in% c("tbl", "data.table"))) 
         IVs.test <- as.data.frame(IVs.test)
       IVs.test <- IVs.test[, sort(features), drop = FALSE]
     }
     
     # Balance classes if requested
     if(balance){
       if (!requireNamespace("caret", quietly = TRUE)) {
         stop("Package 'caret' needed for balancing. Please install it.")
       }
       
       DV.factor <- as.factor(DV.train)
       balanced_data <- caret::downSample(IVs.train, DV.factor)
       balanced_data_up <- caret::upSample(IVs.train, DV.factor)
       
       IVs.train <- rbind(balanced_data[, -ncol(balanced_data)], 
                          balanced_data_up[, -ncol(balanced_data_up)])
       DV.train <- c(as.character(balanced_data$Class), 
                     as.character(balanced_data_up$Class))
       
       # Convert to numeric for classification
       if(!is.null(type)) {
         DV.train <- as.numeric(as.factor(DV.train))
       }
     }
     
     x <- IVs.train
     y <- DV.train
     z <- IVs.test
     
     n <- ncol(x)
     n_train <- nrow(x)
     
     # Set default epochs
     if(is.null(epochs)) epochs <- 2 * n_train
     
     # Set distance metric
     dist <- if(!is.null(ts.test)) "DTW" else "L2"
     
     old.threshold <- 0
     sets <- sum(choose(n, 1:n))
     deterministic <- FALSE 
     
     # Handle deterministic feature sets
     if((sets < n_train) || n <= 10){
       deterministic <- TRUE
       learner.trials <- min(sets, learner.trials)
       
       deterministic.sets <- list()
       for(k in 1:n) {
         deterministic.sets <- c(deterministic.sets, combn(n, k, simplify = FALSE))
       }
     }
     
     # Threshold learning phase
     if(is.null(threshold)){
       if(!extreme) epochs <- NULL
       if(is.null(CV.size)) new.CV.size <- round(runif(1, 0.2, 1/3), 3) else new.CV.size <- CV.size
       
       old.threshold <- 1
       if(is.null(learner.trials)) learner.trials <- n_train
       
       results <- numeric(learner.trials)
       test.features <- vector("list", learner.trials)
       
       for(i in 1:learner.trials){
         set.seed(123 + i)
         
         # Create test index
         if(!is.null(ts.test)) {
           new.index <- 1:(n_train - ts.test)
         } else if(i <= n_train/4) {
           new.index <- round(seq(i, n_train, length.out = round(new.CV.size * n_train)))
         } else {
           new.index <- sample(n_train, round(new.CV.size * n_train), replace = FALSE)
         }
         
         # Create training and test sets
         train.idx <- setdiff(1:n_train, new.index)
         test.idx <- new.index
         
         # Handle feature selection
         if(deterministic && i <= length(deterministic.sets)) {
           features <- deterministic.sets[[i]]
         } else {
           k <- sample(2:min(n, 10), 1)  # Limit feature combinations
           features <- sort(sample(n, k, replace = FALSE))
         }
         test.features[[i]] <- features
         
         # Get subset data
         x.train <- x[train.idx, features, drop = FALSE]
         y.train <- y[train.idx]
         x.test <- x[test.idx, features, drop = FALSE]
         actual <- y[test.idx]
         
         if(status && i %% 10 == 0) {
           message("Threshold trials: ", i, "/", learner.trials, "\r", appendLF = FALSE)
         }
         
         # Run regression
         model <- tryCatch({
           NNS.reg(x.train, y.train, point.est = x.test, 
                   dim.red.method = "equal", plot = FALSE, order = depth,
                   ncores = 1, type = type, smooth = TRUE)
         }, error = function(e) {
           list(Point.est = rep(mean(y.train), length(actual)))
         })
         
         predicted <- model$Point.est
         
         # Handle NAs and out-of-range predictions
         predicted[is.na(predicted)] <- mean(y.train, na.rm = TRUE)
         if(!is.null(type)) {
           predicted <- pmin(pmax(predicted, min(y, na.rm = TRUE)), max(y, na.rm = TRUE))
         }
         
         # Evaluate objective function
         results[i] <- eval(obj.fn)
       }
     } else {
       results <- threshold
     }
     
     # Calculate threshold
     if(extreme){
       threshold <- if(objective == "max") max(results, na.rm = TRUE) else min(results, na.rm = TRUE)
     } else {
       q <- quantile(results, probs = c(0.25, 0.75), na.rm = TRUE)
       threshold <- if(objective == "max") q[2] else q[1]
     }
     
     # Plot feature importance
     if(feature.importance && is.null(threshold)){
       par(mfrow = c(1, 1))
       hist(results, main = "Distribution of Learner Trials Objective Function",
            xlab = "Objective Function", col = "steelblue")
       abline(v = threshold, col = 'red', lty = 2, lwd = 2)
       mtext(round(threshold, 2), side = 1, col = "red", at = threshold)
       legend_text <- if(objective == 'max') "Threshold >" else "< Threshold"
       mtext(legend_text, side = 3, col = "red", at = threshold)
     }
     
     if(status) message("Learner Accuracy Threshold = ", round(threshold, 3))
     
     # Select feature sets based on threshold
     if(extreme){
       idx <- if(objective == "max") which.max(results) else which.min(results)
       reduced.test.features <- list(test.features[[idx]])
     } else {
       if(objective == "max") {
         reduced.test.features <- test.features[results >= threshold]
       } else {
         reduced.test.features <- test.features[results <= threshold]
       }
     }
     
     # Calculate feature frequencies
     feature.freq <- table(unlist(reduced.test.features))
     feature.freq <- feature.freq / sum(feature.freq)
     
     # Epochs training phase
     keeper.features <- list()
     
     if(!is.null(epochs) && !deterministic && length(reduced.test.features) > 0){
       if(is.null(CV.size)) new.CV.size <- round(runif(1, 0.2, 1/3), 3) else new.CV.size <- CV.size
       
       for(j in 1:min(epochs, 1000)) {  # Limit epochs
         set.seed(123 * j)
         
         # Create test index
         if(!is.null(ts.test)) {
           test.idx <- (n_train - 2*ts.test + 1):n_train
         } else if(j <= n_train/4) {
           test.idx <- round(seq(j, n_train, length.out = round(new.CV.size * n_train)))
         } else {
           test.idx <- sample(n_train, round(new.CV.size * n_train), replace = FALSE)
         }
         
         train.idx <- setdiff(1:n_train, test.idx)
         
         if(status && j %% 10 == 0) {
           message("Epochs: ", j, "/", epochs, "\r", appendLF = FALSE)
         }
         
         # Select features
         if(length(reduced.test.features) > 1) {
           base.features <- reduced.test.features[[sample(length(reduced.test.features), 1)]]
         } else {
           base.features <- reduced.test.features[[1]]
         }
         
         # Add random features
         extra.size <- sample(0:min(3, n-length(base.features)), 1)  # Limit extra features
         if(extra.size > 0) {
           extra.features <- sample(setdiff(1:n, base.features), extra.size)
           features <- sort(c(base.features, extra.features))
         } else {
           features <- base.features
         }
         
         # Get subset data
         x.train <- x[train.idx, features, drop = FALSE]
         y.train <- y[train.idx]
         x.test <- x[test.idx, features, drop = FALSE]
         actual <- y[test.idx]
         
         # Run regression
         model <- tryCatch({
           NNS.reg(x.train, y.train, point.est = x.test, 
                   dim.red.method = "equal", plot = FALSE, order = depth,
                   ncores = 1, type = type, smooth = TRUE)
         }, error = function(e) {
           list(Point.est = rep(mean(y.train), length(actual)))
         })
         
         predicted <- model$Point.est
         
         # Handle NAs and out-of-range predictions
         predicted[is.na(predicted)] <- mean(y.train, na.rm = TRUE)
         if(!is.null(type)) {
           predicted <- pmin(pmax(predicted, min(y, na.rm = TRUE)), max(y, na.rm = TRUE))
         }
         
         # Evaluate objective function
         new.results <- eval(obj.fn)
         
         # Keep features if they meet threshold
         if((objective == "max" && new.results >= threshold) ||
            (objective == "min" && new.results <= threshold)) {
           keeper.features[[j]] <- features
         }
       }
     } else {
       keeper.features <- reduced.test.features
     }
     
     # Process keeper features
     keeper.features <- keeper.features[!sapply(keeper.features, is.null)]
     if(length(keeper.features) == 0) keeper.features <- reduced.test.features
     
     # Calculate feature frequencies
     plot.table <- table(unlist(keeper.features))
     names(plot.table) <- colnames(IVs.train)[as.numeric(names(plot.table))]
     plot.table <- sort(plot.table, decreasing = TRUE)
     
     # Return early if only features requested
     if(features.only){
       return(list(feature.weights = plot.table/sum(plot.table),
                   feature.frequency = plot.table))
     }
     
     # Final prediction
     if(status) message("Generating final estimate...")
     
     # Get all unique features used
     all.features <- unique(unlist(keeper.features))
     x.final <- x[, all.features, drop = FALSE]
     z.final <- z[, all.features, drop = FALSE]
     
     # Use NNS.stack for final prediction
     model <- NNS.stack(IVs.train = x.final, 
                        DV.train = y,
                        IVs.test = z.final,
                        order = depth, 
                        dim.red.method = "all",
                        ncores = 1,
                        stack = FALSE, 
                        status = status,
                        type = type, 
                        dist = dist, 
                        folds = folds,
                        pred.int = pred.int)
     
     estimates <- model$stack
     
     # Process estimates
     if(!is.null(type)) {
       # For classification, round to nearest class
       class.levels <- sort(unique(y))
       estimates <- sapply(estimates, function(p) {
         class.levels[which.min(abs(p - class.levels))]
       })
     }
     
     # Final feature importance plot
     if(feature.importance) {
       important.features <- head(plot.table, min(10, length(plot.table)))
       barplot(important.features, horiz = TRUE, las = 1,
               main = "Top Feature Frequencies", 
               xlab = "Frequency", col = "steelblue")
     }
     
     # Return results
     return(list(results = estimates,
                 pred.int = model$pred.int,
                 feature.weights = plot.table/sum(plot.table),
                 feature.frequency = plot.table))
}