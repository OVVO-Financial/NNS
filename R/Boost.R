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
  
  .core <- function() {
    
    if(sum(is.na(cbind(IVs.train,DV.train))) > 0) stop("You have some missing values, please address.")
    if(is.null(obj.fn)) stop("Please provide an objective function")
    
    if(balance && is.null(type)) warning("type = 'CLASS' selected due to balance = TRUE.")
    if(balance) type <- "CLASS"
    
    if(!is.null(type) && min(as.numeric(as.factor(DV.train)))==0) warning("Base response variable category should be 1, not 0.")
    
    if(any(class(IVs.train)%in%c("tbl","data.table"))) IVs.train <- as.data.frame(IVs.train)
    if(any(class(DV.train)%in%c("tbl","data.table"))) DV.train <- as.vector(unlist(DV.train))
    
    if(!is.null(type)){
      type <- tolower(type)
      if(type == "class" && identical(obj.fn,expression( sum((predicted - actual)^2) ))){
        obj.fn <- expression(mean( predicted == as.numeric(actual)))
        objective <- "max"
      }
    }
    
    objective <- tolower(objective)
    
    if(is.null(colnames(IVs.train))){
      colnames.list <- lapply(1 : dim(IVs.train)[2], function(i) paste0("X", i))
      colnames(IVs.test) <- colnames(IVs.train) <- as.character(colnames.list)
    }
    
    features <- colnames(IVs.train)
    IVs.train <- IVs.train[ ,sort(features)]
    
    transform <- data.matrix(cbind(DV.train, IVs.train))
    
    IVs.train <- transform[,-1]
    colnames(IVs.train) <- sort(features)
    
    DV.train <- transform[,1]
    
    if(is.null(IVs.test)){
      IVs.test <- IVs.train
    } else {
      if(any(class(IVs.test)%in%c("tbl","data.table"))) IVs.test <- as.data.frame(IVs.test)
    }
    
    if(balance){
      DV.train <- as.numeric(as.factor(DV.train))
      
      y_train <- as.factor(as.character(DV.train))
      
      training_1 <- do.call(cbind, downSample(IVs.train, y_train, list = TRUE))
      training_2 <- do.call(cbind, upSample(IVs.train, y_train, list = TRUE))
      
      training <- rbind.data.frame(training_1, training_2)
      
      colnames(training) <- c(colnames(IVs.train), names(DV.train))
      
      IVs.train <- training[, -ncol(training)]
      DV.train <- as.numeric(as.character(training[,ncol(training)]))
    }
    
    x <- data.table::data.table(IVs.train)
    y <- DV.train
    z <- IVs.test
    
    ### Representative samples
    yx <- cbind(y, x)
    rep.x <- yx[,lapply(.SD, function(z) fivenum(as.numeric(z)))]
    rm(yx)
    
    rep.y <- unlist(rep.x[,1])
    rep.x <- rep.x[,-1]
    
    rep.x <- as.data.frame(rep.x)
    
    n <- ncol(x)
    
    if(is.null(epochs)) epochs <- 2*length(y)
    
    if(!is.null(ts.test)) dist <- "DTW" else dist <- "L2"
    
    estimates <- list()
    fold <- list()
    
    old.threshold <- 0
    
    sets <- sum(choose(n, 1:n))
    deterministic <- FALSE 
    if((sets < length(y)) || n <= 10){
      deterministic <- TRUE
      learner.trials <- sets
      combn_vec <- Vectorize(Rfast::comb_n, vectorize.args = "k")
      deterministic.sets <- unlist(lapply(combn_vec(n, 1:n), function(df) as.list(as.data.frame(df))), recursive = FALSE)
    }
    
    # Add test loop for highest threshold ...
    if(is.null(threshold)){
      if(!extreme) epochs <- NULL
      if(is.null(CV.size)) new.CV.size <- round(runif(1, .2, 1/3), 3) else new.CV.size <- CV.size
      
      old.threshold <- 1
      
      if(is.null(learner.trials)){learner.trials <- length(y)}
      
      results <- numeric(learner.trials)
      test.features <- vector(mode = "list", learner.trials)
      
      for(i in 1:learner.trials){
        set.seed(123 + i)
        
        l <- length(y)
        if(i<=l/4) new.index <- as.integer(seq(i, length(y), length.out = as.integer(new.CV.size * length(y)))) else {
          new.index <- sample(l, as.integer(new.CV.size * l), replace = FALSE)
        }
        
        if(!is.null(ts.test)) new.index <- 1:(length(y) - ts.test)
        
        new.index <- unlist(new.index)
        
        new.iv.train <- cbind(y[-new.index], x[-new.index,]) 
        new.iv.train <- new.iv.train[,lapply(.SD, as.double)]
        
        new.iv.train <- new.iv.train[,lapply(.SD, function(z) fivenum(as.numeric(z)))]
        
        new.dv.train <- unlist(new.iv.train[,1])
        new.iv.train <- as.data.frame(new.iv.train)
        new.iv.train <- new.iv.train[,unlist(colnames(new.iv.train)%in%colnames(IVs.train))]
        
        new.iv.train <- data.table::rbindlist(list(new.iv.train, x[-new.index,]), use.names = FALSE)
        new.dv.train <- c(new.dv.train, y[-new.index])
        
        colnames(new.iv.train) <- c( features)
        
        actual <- as.numeric(y[new.index])
        new.iv.test <- x[new.index,]
        
        if(status) message("Current Threshold Iterations Remaining = " ,learner.trials+1-i," ","\r",appendLF=FALSE)
        
        if(deterministic) test.features[[i]] <- deterministic.sets[[i]] else test.features[[i]] <- sort(sample(n, sample(2:n, 1), replace = FALSE))
        
        learning.IVs <- new.iv.train[,.SD, .SDcols = unlist(test.features[i])]
        
        #If estimate is > threshold, store 'features'
        predicted <- NNS.reg(learning.IVs,
                             new.dv.train,
                             point.est = new.iv.test[, .SD, .SDcols=unlist(test.features[[i]])],
                             dim.red.method = "equal",
                             plot = FALSE, order = depth,
                             ncores = 1, type = type)$Point.est
        
        predicted[is.na(predicted)] <- gravity(na.omit(predicted))
        
        # Do not predict a new unseen class
        if(!is.null(type)){
          predicted <- pmin(predicted, max(as.numeric(y)))
          predicted <- pmax(predicted, min(as.numeric(y)))
        }
        
        results[i] <- eval(obj.fn)
        
      } # i in learner.trials
    } else {
      results <- threshold
    } # NULL threshold
    
    if(extreme){
      if(objective=="max") threshold <- max(results) else threshold <- min(results)
    } else {
      if(objective=="max") threshold <- fivenum(results)[4] else threshold <- fivenum(results)[2]
    }
    
    if(feature.importance){
      par(mfrow = c(2,1))
      par(mai = c(1.0,.5,0.8,0.5))
      hist(results, main = "Distribution of Learner Trials Objective Function",
           xlab = "Objective Function", col = "steelblue")
      abline(v = threshold, col = 'red', lty = 2, lwd = 2)
      mtext(round(threshold, 2), side = 1, col = "red", at = threshold)
      if(extreme){
        if(objective=='max') mtext("Threshold >", side = 3, col = "red", at = threshold, adj = 1) else mtext("< Threshold", side = 3, col = "red", at = threshold, adj = 0)
      } else {
        if(objective=='max') mtext("Threshold >", side = 3, col = "red", at = threshold) else mtext("< Threshold", side = 3, col = "red", at = threshold)
      }
    }
    
    if(status){
      message(paste0("Learner Accuracy Threshold = ", format(threshold, digits = 3, nsmall = 2),"           "), appendLF = TRUE)
      # Clear message line
      message("                                       ", "\r", appendLF = FALSE)
    }
    
    if(extreme){
      if(objective=="max") reduced.test.features <- test.features[which.max(results)] else reduced.test.features <- test.features[which.min(results)]
    } else {
      if(objective=="max") reduced.test.features <- test.features[which(results>=threshold)] else reduced.test.features <- test.features[which(results<=threshold)]
    }
    
    rf <- data.table::data.table(table(as.character(reduced.test.features)))
    rf$N <- rf$N / sum(rf$N)
    
    rf_reduced <- apply(rf, 1, function(x) eval(parse(text=x[1])))
    
    scale_factor_rf <- table(unlist(rf_reduced))/min(table(unlist(rf_reduced)))
    
    reduced.test.features <- as.numeric(rep(names(scale_factor_rf), ifelse(scale_factor_rf%%1 < .5, floor(scale_factor_rf), ceiling(scale_factor_rf))))
    
    keeper.features <- list()
    
    if(deterministic) epochs <- NULL
    
    if(!is.null(epochs) && !deterministic){
      
      if(is.null(CV.size)) new.CV.size <- round(runif(1, .2, 1/3), 3) else new.CV.size <- CV.size
      
      for(j in 1:epochs){
        set.seed(123 * j)
        
        l <- length(y)
        if(j<=l/4) new.index <- as.integer(seq(j, length(y), length.out = as.integer(new.CV.size * length(y)))) else {
          new.index <- sample(l, as.integer(new.CV.size * l), replace = FALSE)
        }
        if(!is.null(ts.test)) new.index <- length(y) - (2*ts.test):0
        
        new.index <- unlist(new.index)
        
        new.iv.train <- cbind(y[-new.index], x[-new.index, ])
        new.iv.train <- new.iv.train[, lapply(.SD, as.double)]
        
        new.iv.train <- new.iv.train[,lapply(.SD,fivenum), by = .(y[-new.index])]
        
        new.dv.train <- unlist(new.iv.train[, 1])
        new.iv.train <- as.data.frame(new.iv.train)
        new.iv.train <- new.iv.train[,unlist(colnames(new.iv.train)%in%colnames(IVs.train))]
        
        new.iv.train <- data.table::rbindlist(list(new.iv.train, x[-new.index,]), use.names = FALSE)
        new.dv.train <- c(new.dv.train, y[-new.index])
        
        actual <- as.numeric(y[new.index])
        new.iv.test <- x[new.index,]
        
        if(status){
          message("% of epochs = ", format(j/epochs,digits =  3,nsmall = 2),"     ","\r",appendLF=FALSE)
          if(j == epochs){
            message("% of epochs ",j," = 1.000     ","\r",appendLF = FALSE)
            flush.console()
          }
        }
        
        if(deterministic) features <- unlist(deterministic.sets[[j]]) else features <- sort(c(unlist(reduced.test.features), sample(c(1:n), sample(1:n, 1), replace = FALSE)))
        
        if(length(features) == 1) point.est.values <- unlist(new.iv.test[, as.numeric(features)]) else point.est.values <- new.iv.test[, as.numeric(features)]
        
        #If estimate is > threshold, store 'features'
        predicted <- NNS.reg(new.iv.train[, as.numeric(features)],
                             new.dv.train, point.est = point.est.values,
                             dim.red.method = "equal",
                             plot = FALSE, residual.plot = FALSE, order = depth,
                             ncores = 1, type = type)$Point.est
        
        predicted[is.na(predicted)] <- gravity(na.omit(predicted))
        # Do not predict a new unseen class
        if(!is.null(type)){
          predicted <- pmin(predicted,max(as.numeric(y)))
          predicted <- pmax(predicted,min(as.numeric(y)))
        }
        
        new.results <- eval(obj.fn)
        
        if(objective=="max"){
          if(is.na(new.results)) new.results <- .99*threshold
          if(new.results>=threshold) keeper.features[[j]] <- features else keeper.features[[j]] <- NULL
        } else {
          if(is.na(new.results)) new.results <- 1.01*threshold
          if(new.results<=threshold) keeper.features[[j]] <- features else keeper.features[[j]] <- NULL
        }
      }
    } else { # !is.null(epochs)
      keeper.features <- reduced.test.features
    }
    
    keeper.features <- keeper.features[!sapply(keeper.features, is.null)]
    if(length(keeper.features)==0){
      if(old.threshold==0){
        if(objective=="min") stop("Please increase [threshold].") else stop("Please reduce [threshold].")
      } else {
        keeper.features <- test.features[which.max(results)]
      }
    }
    
    plot.table <- table(unlist(keeper.features))
    
    names(plot.table) <- colnames(IVs.train)[eval(as.numeric(names(plot.table)))]
    
    if(features.only || feature.importance) plot.table <- plot.table[rev(order(plot.table))]
    
    if(features.only){
      return(list("feature.weights" = plot.table/sum(plot.table),
                  "feature.frequency" = plot.table))
    }
    
    if(!is.null(rep.y)){
      x <- rbind(rep.x, (x))
      y <- c(rep.y, y)
    }
    
    kf <- data.table::data.table(table(as.character(keeper.features)))
    kf$N <- kf$N / sum(kf$N)
    
    kf_reduced <- apply(kf, 1, function(x) eval(parse(text=x[1])))
    
    scale_factor <- table(unlist(kf_reduced))/min(table(unlist(kf_reduced)))
    
    final_scale <- as.numeric(rep(names(scale_factor), ifelse(scale_factor%%1 < .5, floor(scale_factor), ceiling(scale_factor))))
    
    if(status) message("Generating Final Estimate" ,"\r", appendLF = TRUE)
    
    model <- NNS.stack(x[, keeper.features],
                       y,
                       IVs.test = z[, keeper.features],
                       order = depth, dim.red.method = "all",
                       ncores = 1,
                       stack = FALSE, status = status,
                       type = type, dist = dist, folds = folds,
                       pred.int = pred.int)
    
    estimates <- model$stack
    
    estimates[is.na(unlist(estimates))] <- ifelse(!is.null(type), mode_class(unlist(na.omit(estimates))), mode(unlist(na.omit(estimates))))
    
    if(!is.null(type)){
      estimates <- pmin(estimates, max(as.numeric(y)))
      estimates <- pmax(estimates, min(as.numeric(y)))
    }
    
    if(feature.importance){
      linch <-  max(strwidth(names(plot.table), "inch") + 0.4, na.rm = TRUE)
      par(mai=c(1.0, linch, 0.8, 0.5))
      
      if(length(plot.table)!=1){
        barplot(sort(plot.table, decreasing = FALSE)[1:min(n, 10)],
                horiz = TRUE,
                col='steelblue',
                main="Feature Frequency in Final Estimate",
                xlab = "Frequency",las=1)
      } else {
        barplot(sort(plot.table,decreasing = FALSE),
                horiz = TRUE,
                col='steelblue',
                main="Feature Frequency in Final Estimate",
                xlab = "Frequency", las = 1)
      }
      par(mfrow=c(1,1))
    }
    
    if(!is.null(type)) estimates <- ifelse(estimates%%1 < .5, floor(estimates), ceiling(estimates))
    
    return(list("results" = estimates,
                "pred.int" = model$pred.int,
                "feature.weights" = plot.table/sum(plot.table),
                "feature.frequency" = plot.table))
  } # end .core
  
  # First attempt: suppress all warnings; on error, retry once with balance=FALSE
  out <- tryCatch(
    suppressWarnings(.core()),
    error = function(e) {
      if (isTRUE(balance)) {
        warning("[retry] First attempt failed; retrying with balance = FALSE")
        return(Recall(IVs.train = IVs.train,
                      DV.train = DV.train,
                      IVs.test = IVs.test,
                      type = type,
                      depth = depth,
                      learner.trials = learner.trials,
                      epochs = epochs,
                      CV.size = CV.size,
                      balance = FALSE,          # single retry without balance
                      ts.test = ts.test,
                      folds = folds,
                      threshold = threshold,
                      obj.fn = obj.fn,
                      objective = objective,
                      extreme = extreme,
                      features.only = features.only,
                      feature.importance = feature.importance,
                      pred.int = pred.int,
                      status = status))
      }
      stop(e)
    }
  )
  
  out
}
