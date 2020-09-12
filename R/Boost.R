#' NNS Boost
#'
#' Ensemble method for classification using the predictions of the NNS multivariate regression \link{NNS.reg} collected from uncorrelated feature combinations.
#'
#' @param IVs.train a matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimensions to \code{(IVs.train)}.
#' @param IVs.test a matrix or data frame of variables of numeric or factor data types with compatible dimensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.
#' @param type \code{NULL} (default).  To perform a classification of discrete integer classes from factor target variable \code{(DV.train)}, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.
#' @param representative.sample logical; \code{FALSE} (default) Reduces observations of \code{IVs.train} to a set of representative observations per regressor.
#' @param depth options: (integer, NULL, "max"); Specifies the \code{order} parameter in the \link{NNS.reg} routine, assigning a number of splits in the regressors.  \code{(depth = "max")}(default) will be significantly faster, but increase the variance of results, which is suggested for mixed continuous and discrete (unordered, ordered) data.
#' @param n.best integer; \code{NULL} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at \code{sqrt(# of regressors)}. Analogous to \code{k} in a \code{k Nearest Neighbors} algorithm.  If \code{NULL}, determines the optimal clusters via the \link{NNS.stack} procedure.
#' @param learner.trials integer; \code{NULL} (default) Sets the number of trials to obtain an accuracy \code{threshold} level.  \code{(learner.trials = 100)} is the default setting.
#' @param epochs integer; \code{2*length(DV.train)} (default) Total number of feature combinations to run.
#' @param CV.size numeric [0, 1]; \code{(CV.size = .25)} (default) Sets the cross-validation size.  Defaults to 0.25 for a 25 percent random sampling of the training set.
#' @param balance logical; \code{FALSE} (default) Uses both up and down sampling from \code{caret} to balance the classes.  \code{type="CLASS"} required.
#' @param ts.test integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.
#' @param folds integer; 5 (default) Sets the number of \code{folds} in the \link{NNS.stack} procedure for optimal \code{n.best} parameter.
#' @param threshold numeric; \code{NULL} (default) Sets the \code{obj.fn} threshold to keep feature combinations.
#' @param obj.fn expression;
#' \code{expression( sum((predicted - actual)^2) )} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.  Automatically selects an accuracy measure when \code{(type = "CLASS")}.
#' @param objective options: ("min", "max") \code{"max"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param extreme logical; \code{FALSE} (default) Uses the maximum (minimum) \code{threshold} obtained from the \code{learner.trials}, rather than the upper (lower) quintile level for maximization (minimization) \code{objective}.
#' @param feature.importance logical; \code{TRUE} (default) Plots the frequency of features used in the final estimate.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns a vector of fitted values for the dependent variable test set \code{$results}, and the final feature loadings \code{$feature.weights}.
#'
#' @note Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#' \url{https://www.ssrn.com/abstract=2864711}
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  a <- NNS.boost(iris[1:140, 1:4], iris[1:140, 5],
#'  IVs.test = iris[141:150, 1:4],
#'  epochs = 100, learner.trials = 100,
#'  type = "CLASS")
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
                      representative.sample = FALSE,
                      depth = "max",
                      n.best = NULL,
                      learner.trials = 100,
                      epochs = NULL,
                      CV.size = .25,
                      balance = FALSE,
                      ts.test = NULL,
                      folds = 5,
                      threshold = NULL,
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = "min",
                      extreme = FALSE,
                      feature.importance = TRUE,
                      status = TRUE,
                      ncores = NULL){

  if(is.null(obj.fn)){ stop("Please provide an objective function")}

  if(!is.null(type)){
    type <- tolower(type)
    if(type=="class"){
      obj.fn <- expression( mean(predicted == as.numeric(actual)) )
      objective <- "max"
    }
  }

  objective <- tolower(objective)



  # Parallel process...
  if (is.null(ncores)) {
    cores <- detectCores()
    num_cores <- cores - 1
  } else {
    cores <- detectCores()
    num_cores <- ncores
  }

  if((num_cores)>cores){ stop(paste0("Please ensure total number of cores [ncores] is less than ", cores))}
  if(balance && type!="class") stop("Please select 'type='CLASS'' when 'balance=TRUE'.")

  if(balance && type=="class"){
    if(1%in%DV.train){
      DV.train <- as.numeric(as.character(DV.train))
    } else {
      DV.train <- as.numeric(as.character(as.numeric(DV.train)))
    }

    y_train <- as.factor(as.character(DV.train))

    training_1 <- do.call(cbind, caret::downSample(IVs.train, y_train, list = TRUE))
    training_2 <- do.call(cbind, caret::upSample(IVs.train, y_train, list = TRUE))

    training <- rbind.data.frame(training_1, training_2)

    colnames(training) <- c(colnames(IVs.train), names(DV.train))

    IVs.train <- training[, -ncol(training)]
    DV.train <- as.numeric(as.character(training[,ncol(training)]))


    if(is.null(IVs.test)) {IVs.test <- IVs.train}
  }

  if(is.null(IVs.test)) {IVs.test <- IVs.train}

  x <- IVs.train
  y <- DV.train
  z <- IVs.test


  ### Representative samples
  rep.x <- data.table::data.table(x)

  fivenum.x <- rep.x[,lapply(.SD, function(z) fivenum(as.numeric(z))), by = .(y)]
  mode.x <- rep.x[,lapply(.SD, function(z) mode(as.numeric(z))), by = .(y)]
  mean.x <- rep.x[,lapply(.SD, function(z) mean(as.numeric(z))), by = .(y)]

  rep.x <- data.table::rbindlist(list(fivenum.x, mode.x, mean.x), use.names = FALSE)

  rep.y <- unlist(rep.x[,1])
  rep.x <- rep.x[,-1]

  rep.x <- as.data.frame(rep.x)

  n <- ncol(x)

  if(is.null(epochs)){
    epochs <- 2*length(y)
  }

  if(!is.null(ts.test)){
    dist <- "DTW"
  } else {
    dist <- "L2"
  }

  estimates <- list()
  fold <- list()

  old.threshold <- 0

  if(is.null(n.best)){
    if(status){
      message("Currently determining optimal [n.best] clusters...","\r",appendLF=TRUE)
    }

    n.best <- NNS.stack(x, y, folds = 1, status = status,
                        method = 1, order = depth,
                        obj.fn = obj.fn, ts.test = ts.test,
                        objective = objective,
                        ncores = ncores, type = type)$NNS.reg.n.best

    if(status){
      message("Currently determining learning threshold...","\r",appendLF=TRUE)
    }
  }

  # Add test loop for highest threshold ...
  if(is.null(threshold)){
    if(!extreme){
      epochs <- NULL
    }
    old.threshold <- 1

    if(is.null(learner.trials)){learner.trials <- length(y)}

    learner.trials <- min(sum(choose(n, 1:n)), learner.trials)

    results <- numeric(learner.trials)
    test.features <- list(learner.trials)

    for(i in 1:learner.trials){
      set.seed(123 + i)
      new.index <- sample(length(y), as.integer(CV.size*length(y)), replace = FALSE)

      if(i > 1){
        maxes <- as.vector(apply(x, 2, which.max))
        mins <- as.vector(apply(x, 2, which.min))
        new.index_half <- new.index.1[1:(length(new.index.1)/2)]
        new.index <- na.omit(unique(c(mins, maxes, new.index_half, new.index))[1:as.integer(CV.size*length(y))])
      }

      if(!is.null(ts.test)){
        new.index <- 1:(length(y) - ts.test)
      }

      new.index <- unlist(new.index)
      new.iv.train <- data.table::data.table(x[-new.index,])
      new.iv.train <- new.iv.train[,lapply(.SD, as.double)]

      fivenum.new.iv.train <- new.iv.train[,lapply(.SD,function(z) fivenum(as.numeric(z))), by = .(y[-new.index])]
      mode.new.iv.train <- new.iv.train[,lapply(.SD,function(z) mode(as.numeric(z))), by = .(y[-new.index])]
      mean.new.iv.train <- new.iv.train[,lapply(.SD,function(z) mean(as.numeric(z))), by = .(y[-new.index])]

      new.iv.train <- data.table::rbindlist(list(fivenum.new.iv.train,mode.new.iv.train,mean.new.iv.train), use.names = FALSE)

      new.iv.train <- as.data.frame(new.iv.train[,-1])
      new.dv.train <- unlist(new.iv.train[,1])

      if(!representative.sample){
        new.iv.train <- rbind(new.iv.train, data.matrix(x[-new.index,]))
        new.dv.train <- c(new.dv.train, y[-new.index])
      }

      actual <- as.numeric(y[new.index])
      new.iv.test <- x[new.index,]


      if(status){
        message("Current Threshold Iterations Remaining = " ,learner.trials+1-i," ","\r",appendLF=FALSE)
      }

      test.features[[i]] <- sort(sample(n, sample(2:n,1), replace = FALSE))

      #If estimate is > threshold, store 'features'
      predicted <- NNS.reg(new.iv.train[,test.features[[i]]],
                           new.dv.train,
                           point.est = new.iv.test[,test.features[[i]]],
                           plot = FALSE, residual.plot = FALSE, order = depth,
                           n.best = n.best, factor.2.dummy = FALSE,
                           ncores = 1, type = type)$Point.est

      # Do not predict a new unseen class
      if(!is.null(type)){
          predicted <- pmin(predicted, max(as.numeric(y)))
          predicted <- pmax(predicted, min(as.numeric(y)))
      }

      new.index.1 <- rev(order(abs(predicted - actual)))

      results[i] <- eval(obj.fn)

    } # i in learner.trials
  } else {
      results <- threshold
  } # NULL threshold


  if(extreme){
    if(objective=="max"){
      threshold <- max(results)
    } else {
      threshold <- min(results)
    }
  } else {
    if(objective=="max"){
      threshold <- fivenum(results)[4]
    } else {
      threshold <- fivenum(results)[2]
    }
  }

  if(feature.importance){
    original.par <- par(no.readonly = TRUE)
    par(mfrow = c(2,1))
    par(mai = c(1.0,.5,0.8,0.5))
    hist(results, main = "Distribution of Learner Trials Objective Function",
         xlab = "Objective Function", col = "steelblue")
    abline(v = threshold, col = 'red', lty = 2, lwd = 2)
    mtext(round(threshold, 2), side = 1, col = "red", at = threshold)
    if(extreme){
      if(objective=='max'){
        mtext("Threshold >", side = 3, col = "red", at = threshold, adj = 1)
      } else {
        mtext("< Threshold", side = 3, col = "red", at = threshold, adj = 0)
      }
    } else {
      if(objective=='max'){
        mtext("Threshold >", side = 3, col = "red", at = threshold)
      } else {
        mtext("< Threshold", side = 3, col = "red", at = threshold)
      }
    }
  }



  if(status){
    message(paste0("Learner Accuracy Threshold = ", format(threshold, digits = 3, nsmall = 2),"           "), appendLF = TRUE)

    # Clear message line
    message("                                       ", "\r", appendLF = FALSE)
  }


  keeper.features <- list()

  if(!is.null(epochs)){
    for(j in 1:epochs){
      set.seed(123 * j)
      new.index <- sample(length(y), as.integer(CV.size*length(y)), replace = FALSE)

      if(j > 1){
        maxes <- as.vector(apply(x, 2, which.max))
        mins <- as.vector(apply(x, 2, which.min))
        new.index_half <- new.index.1[1:(length(new.index.1)/2)]
        new.index <- na.omit(unique(c(mins, maxes, new.index_half, new.index))[1:as.integer(CV.size*length(y))])
      }

      if(!is.null(ts.test)){
        new.index <- length(y) - (2*ts.test):0
      }

      new.index <- unlist(new.index)
      new.iv.train <- data.table::data.table(x[-new.index, ])
      new.iv.train <- new.iv.train[, lapply(.SD,as.double)]

      fivenum.new.iv.train <- new.iv.train[,lapply(.SD,fivenum), by = .(y[-new.index])]
      mode.new.iv.train <- new.iv.train[,lapply(.SD,function(z) mode(as.numeric(z))), by = .(y[-new.index])]
      mean.new.iv.train <- new.iv.train[,lapply(.SD,function(z) mean(as.numeric(z))), by = .(y[-new.index])]


      names(fivenum.new.iv.train) <- c("y", colnames(new.iv.train))
      names(mean.new.iv.train) <- c("y", colnames(new.iv.train))
      names(mode.new.iv.train) <- c("y", colnames(new.iv.train))


      new.iv.train <- data.table::rbindlist(list(fivenum.new.iv.train, mode.new.iv.train, mean.new.iv.train), use.names = FALSE)

      new.iv.train <- as.data.frame(new.iv.train[, -1])
      new.dv.train <- unlist(new.iv.train[, 1])

      if(!representative.sample){
        new.iv.train <- rbind(new.iv.train, data.matrix(x[-new.index,]))
        new.dv.train <- c(new.dv.train, y[-new.index])
      }


      actual <- as.numeric(y[new.index])
      new.iv.test <- x[new.index,]

      if(status){
        message("% of epochs = ", format(j/epochs,digits =  3,nsmall = 2),"     ","\r",appendLF=FALSE)

        if(j == epochs){
          message("% of epochs ",j," = 1.000     ","\r",appendLF = FALSE)
          flush.console()
        }
      }

      features <- sort(sample(n, sample(2:n, 1), replace = FALSE))

      #If estimate is > threshold, store 'features'
      predicted <- NNS.reg(new.iv.train[, features],
                           new.dv.train, point.est = new.iv.test[, features],
                           plot = FALSE, residual.plot = FALSE, order = depth, n.best = n.best,
                           factor.2.dummy = FALSE, ncores = 1, type = type)$Point.est

      # Do not predict a new unseen class
      if(!is.null(type)){
          predicted <- pmin(predicted,max(as.numeric(y)))
          predicted <- pmax(predicted,min(as.numeric(y)))
      }

      new.index.1 <- rev(order(abs(predicted - actual)))

      new.results <- eval(obj.fn)

      if(objective=="max"){
        if(new.results>=threshold){
          keeper.features[[j]] <- features
        } else {
          keeper.features[[j]] <- NULL
        }
      } else {
        if(new.results<=threshold){
          keeper.features[[j]] <- features
        } else {
          keeper.features[[j]] <- NULL
        }
      }
    }
  } else { # !is.null(epochs)
    if(objective=="max"){
      keeper.features <- test.features[which(results>=threshold)]
    } else {
      keeper.features <- test.features[which(results<=threshold)]
    }
  }

  keeper.features <- keeper.features[!sapply(keeper.features, is.null)]
  if(length(keeper.features)==0){
    if(old.threshold==0){
      if(objective=="min"){
          stop("Please increase [threshold].")
      } else {
          stop("Please reduce [threshold].")
      }
    } else {
      keeper.features <- test.features[which.max(results)]
    }
  }

  x <- rbind(rep.x, data.matrix(x))
  y <- c(rep.y, y)

  kf <- data.table::data.table(table(as.character(keeper.features)))
  kf$N <- kf$N / sum(kf$N)



    for(i in 1:dim(kf)[1]){

      if(status){
        message("% of Final Estimate = ", format(i/dim(kf)[1], digits =  3, nsmall = 2),"     ","\r", appendLF = FALSE)
      }


      estimates[[i]] <- NNS.reg(data.matrix(x[, eval(parse(text=kf$V1[i]))]), y,
                                point.est = data.matrix(z[, eval(parse(text=kf$V1[i]))]),
                                plot = FALSE, residual.plot = FALSE, order = depth,
                                n.best = n.best,
                                factor.2.dummy = FALSE, ncores = ncores,
                                type = type, dist = dist)$Point.est/dim(kf)[1]

    }


  estimates <- Reduce("+", estimates)

  if(!is.null(type)){
      estimates <- pmin(estimates, max(as.numeric(y)))
      estimates <- pmax(estimates, min(as.numeric(y)))
  }

  plot.table <- table(unlist(keeper.features))
  names(plot.table) <- colnames(x[as.numeric(names(plot.table))])

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
      barplot(sort(plot.table,decreasing = FALSE)[1:min(n,10)],
              horiz = TRUE,
              col='steelblue',
              main="Feature Frequency in Final Estimate",
              xlab = "Frequency", las = 1)
    }
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    par(original.par)
  }
  gc()
  if(is.null(type)){
    return(list("results" = estimates,
                "feature.weights" = plot.table/sum(plot.table)))
  } else {
    estimates <- ifelse(estimates%%1 < 0.5, floor(estimates), ceiling(estimates))
    return(list("results" = estimates,
                "feature.weights" = plot.table/sum(plot.table)))
  }
}
