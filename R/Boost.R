#' NNS Boost
#'
#' Ensemble method for classification using the predictions of the NNS multivariate regression \link{NNS.reg} collected from uncorrelated feature combinations.
#'
#' @param IVs.train a matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param IVs.test a matrix or data frame of variables of numeric or factor data types with compatible dimsensions to \code{(IVs.train)}.
#' @param representative.sample logical; \code{TRUE} (default) Reduces observations of \code{IVs.train} to a set of representative observations per regressor.
#' @param depth integer; \code{NULL} (default) Specifies the \code{order} parameter in the \code{NNS.reg} routine, assigning a number of splits in the regressors.  \code{(depth = "max")} will be signifcantly faster, but increase the variance of results.
#' @param n.best integer; \code{3} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at \code{sqrt(# of regressors)}. Analogous to \code{k} in a \code{k Nearest Neighbors} algorithm.
#' @param learner.trials integer; \code{NULL} (default) Sets the number of trials to obtain an accuracy \code{threshold} level.  Number of observations in the training set is the default setting.
#' @param epochs integer; \code{2*length(DV.train)} (default) Total number of feature combinations to run.
#' @param CV.size numeric [0, 1]; \code{(CV.size = .25)} (default) Sets the cross-validation size.  Defaults to 0.25 for a 25 percent random sampling of the training set.
#' @param threshold numeric [0, 1]; \code{NULL} (default) Sets the \code{obj.fn} accuracy threshold to keep feature combinations.
#' @param obj.fn expression;
#' \code{expression(mean(round(predicted)==as.numeric(actual)))} (default) Mean accuracy is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"max"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param extreme logical; \code{FALSE} (default) Uses the maximum (minimum) \code{threshold} obtained from the \code{learner.trials}, rather than the upper (lower) quintile level for maximization (minimization) \code{objective}.
#' @param feature.importance logical; \code{TRUE} (default) Plots the frequency of features used in the final estimate.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized procedure. If NULL (default), the number of cores to be used is equal to half the number of cores of the machine.
#' @param subcores integer; value specifying the number of cores to be used in the parallelized procedure in the subroutine \link{NNS.reg}.  If NULL (default), the number of cores to be used is equal to half the number of cores of the machine - 1.
#'
#' @return Returns a vector of fitted values for the dependent variable test set \code{$results}, and the final feature loadings \code{$feature.weights}.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#' \url{https://ssrn.com/abstract=2864711}
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  a <- NNS.boost(iris[1:140, 1:4], iris[1:140, 5],
#'  IVs.test = iris[141:150, 1:4],
#'  epochs = 100, learner.trials = 100)
#'
#'  ## Test accuracy
#'  mean(round(a$results)==as.numeric(iris[141:150,5]))
#'  }
#'
#' @export


NNS.boost <- function(IVs.train,
                      DV.train,
                      IVs.test,
                      representative.sample = TRUE,
                      depth = NULL,
                      n.best = 3,
                      learner.trials = NULL,
                      epochs = NULL,
                      CV.size=.25,
                      threshold = NULL,
                      obj.fn = expression(mean(round(predicted)==as.numeric(actual))),
                      objective = "max",
                      extreme = FALSE,
                      feature.importance = TRUE,
                      status = TRUE,
                      ncores = NULL, subcores = NULL){



# Parallel process...
  if (is.null(ncores)) {
      cores <- detectCores()
      num_cores <- as.integer(cores / 2)
  } else {
      cores <- detectCores()
      num_cores <- ncores
  }

  if (is.null(subcores)) {
      subcores <- as.integer(cores / 2) - 1
  }

  if((num_cores+subcores)>cores){ stop(paste0("Please ensure total number of cores [ncores + subcores] is less than ", cores))}


  if(num_cores>1){
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
  } else { cl <- NULL }

  x <- IVs.train
  y <- DV.train
  z <- IVs.test


  if(!is.null(dim(x))){
      if(!is.numeric(x)){
          x <- sapply(x,factor_2_dummy)
    } else {
          x <- apply(x,2,as.double)
    }
    if(is.list(x)){
        x <- do.call(cbind,x)
        x <- apply(x,2,as.double)
    }

  } else {
      x <- factor_2_dummy(x)
      if(is.null(dim(x))){
          x <- as.double(x)
      } else {
          x <- apply(x,2,as.double)
      }
  }



  if(is.null(colnames(x))) {colnames(x) <- colnames(x, do.NULL = FALSE)}
  colnames(x) <- make.unique(colnames(x),sep = "_")

  if(!is.null(dim(z))){
      if(!is.numeric(z)){
          z <- sapply(z,factor_2_dummy)
      } else {
          z <- apply(z,2,as.double)
      }
      if(is.list(z)){
          z <- do.call(cbind,z)
          z <- apply(z,2,as.double)
      }
  } else {
      z <- factor_2_dummy(z)
      if(is.null(dim(z))){
          z <- as.double(z)
      } else {
          z <- apply(z,2,as.double)
      }
  }
  if(is.null(colnames(z))) {colnames(z) <- colnames(z, do.NULL = FALSE)}
  colnames(z) <- make.unique(colnames(z),sep = "_")

      y <- as.double(as.numeric(unlist(y)))


    ### Representative samples
      rep.x <- data.table(x)

      fivenum.x <- rep.x[,lapply(.SD,fivenum), by = .(y)]
      mode.x <- rep.x[,lapply(.SD,mode), by = .(y)]
      mean.x <- rep.x[,lapply(.SD,mean), by = .(y)]

      rep.x <- rbind(fivenum.x,mode.x,mean.x)
      rep.y <- unlist(rep.x[,1])
      rep.x <- rep.x[,-1]

      if(dim(t(t(x)))[2]!=dim(rep.x)[2]){
        Missing <- setdiff(colnames(x),colnames(rep.x))
        if(length(Missing)>0){
            rep.x[Missing] <- 0
            rep.x <- rep.x[colnames(x)]
        }
      }

      rep.x <- as.data.frame(rep.x)


  n <- ncol(x)

  if(is.null(epochs)){
      epochs <- 2*length(y)
  }

  estimates <- list()
  fold <- list()

  old.threshold <- 0


  # Add test loop for highest threshold ...
  if(is.null(threshold)){
      old.threshold <- 1
      test.features <- list()
      results <- numeric()

      if(is.null(learner.trials)){learner.trials <- length(y)}

      for(i in 1:learner.trials){
          set.seed(123 + i)
          new.index <- sample(length(y), as.integer(CV.size*length(y)), replace = FALSE)

          if(i > 1){
            new.index_half <- new.index.1[1:(length(y)/2)]
            new.index <- na.omit(unique(c(new.index_half,new.index))[1:length(y)])
          }


          new.iv.train <- data.table(x[-new.index,])
          new.iv.train <- new.iv.train[,lapply(.SD,as.double)]

          fivenum.new.iv.train <- new.iv.train[,lapply(.SD,fivenum), by = .(y[-new.index])]
          mode.new.iv.train <- new.iv.train[,lapply(.SD,mode), by = .(y[-new.index])]
          mean.new.iv.train <- new.iv.train[,lapply(.SD,mean), by = .(y[-new.index])]

          new.iv.train <- rbind(fivenum.new.iv.train,mode.new.iv.train,mean.new.iv.train)

          new.iv.train <- as.data.frame(new.iv.train[,-1])
          new.dv.train <- unlist(new.iv.train[,1])

          if(!representative.sample){
              new.iv.train <- rbind(new.iv.train,x[-new.index,])
              new.dv.train <- c(new.dv.train,y[-new.index])
          }

          actual <- y[new.index]
          new.iv.test <- x[new.index,]

          if(dim(new.iv.train)[2]!=dim(new.iv.test)[2]){
              Missing <- setdiff(colnames(new.iv.train),colnames(new.iv.test))
              if(length(Missing)>0){
                  new.iv.test[Missing] <- 0
                  new.iv.test <- new.iv.test[colnames(new.iv.train)]
              }
          }

          if(status){
              message("Current Threshold Iterations Remaining = " ,learner.trials+1-i," ","\r",appendLF=FALSE)
          }

          test.features[[i]] <- sort(sample(n,sample(2:n,1),replace = FALSE))

          #If estimate is > threshold, store 'features'
          predicted <- NNS.reg(new.iv.train[,test.features[[i]]],new.dv.train,point.est = new.iv.test[,test.features[[i]]],
                        plot=FALSE, residual.plot = FALSE, order=depth, n.best=n.best,
                        norm="std", factor.2.dummy = FALSE,ncores=subcores)$Point.est

          # Do not predict a new unseen class
          predicted <- pmin(predicted,max(as.numeric(y)))
          predicted <- pmax(predicted,min(as.numeric(y)))


          new.index.1 <- rev(order(abs(predicted - actual)))

      results[i] <- eval(obj.fn)

    } # i in learner.trials
  } # NULL thresholde

    if(feature.importance){
          original.par <- par(no.readonly = TRUE)
          par(mfrow=c(2,1))
          par(mai=c(1.0,.5,0.8,0.5))
          hist(results,main = "Distribution of Learner Trials Accuracy",
               xlab = "Accuracy",col = "steelblue")
    }

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

    if(status){
        message(paste0("Learner Accuracy Threshold = ", format(threshold,digits = 3,nsmall = 2),"           "),appendLF = TRUE)

        # Clear message line
        message("                                       ","\r",appendLF=FALSE)
    }


      keeper.features <- list()


      for(j in 1:epochs){
          set.seed(123 * j)
          new.index <- sample(1:length(y),as.integer(CV.size*length(y)),replace = FALSE)

      if(i > 1){
          new.index_half <- new.index.1[1:(length(y)/2)]
          new.index <- na.omit(unique(c(new.index_half,new.index))[1:length(y)])
      }


      new.iv.train <- data.table(x[-new.index,])
      new.iv.train <- new.iv.train[,lapply(.SD,as.double)]

      fivenum.new.iv.train <- new.iv.train[,lapply(.SD,fivenum), by = .(y[-new.index])]
      mode.new.iv.train <- new.iv.train[,lapply(.SD,mode), by = .(y[-new.index])]
      mean.new.iv.train <- new.iv.train[,lapply(.SD,mean), by = .(y[-new.index])]

      new.iv.train <- rbind(fivenum.new.iv.train,mode.new.iv.train,mean.new.iv.train)


        new.iv.train <- as.data.frame(new.iv.train[,-1])
        new.dv.train <- unlist(new.iv.train[,1])

        if(!representative.sample){
            new.iv.train <- rbind(new.iv.train,x[-new.index,])
            new.dv.train <- c(new.dv.train,y[-new.index])
        }



      actual <- as.numeric(y[new.index])
      new.iv.test <- x[new.index,]

          if(status){
              message("% of epochs = ", format(j/epochs,digits =  3,nsmall = 2),"     ","\r",appendLF=FALSE)

              if(j == epochs){
                  message("% of epochs ",j," = 1.000     ","\r",appendLF=FALSE)
                  flush.console()
              }
          }

          features <- sort(sample(n,sample(2:n,1),replace = FALSE))

          #If estimate is > threshold, store 'features'
          predicted <- NNS.reg(new.iv.train[,features],new.dv.train,point.est = new.iv.test[,features],
                            plot=FALSE, residual.plot = FALSE, order=depth, n.best=n.best,
                            norm="std", factor.2.dummy = FALSE, ncores=subcores)$Point.est

          # Do not predict a new unseen class
          predicted <- pmin(predicted,max(as.numeric(y)))
          predicted <- pmax(predicted,min(as.numeric(y)))

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

      keeper.features <- keeper.features[!sapply(keeper.features, is.null)]
      if(length(keeper.features)==0){
          if(old.threshold==0){
                stop("Please reduce [threshold].")
          } else {
                keeper.features <- test.features[which.max(results)]
          }
      }

      x <- rbind(rep.x,x)
      y <- c(rep.y,y)

      if(!is.null(cl)){
          clusterExport(cl,c("x","y"))
          if(status){
              message("Parallel process running, status unavailable...","\r",appendLF=FALSE)
          }


          estimates <- foreach(i = 1:length(keeper.features))%dopar%{

              NNS.reg(x[,keeper.features[[i]]],y,point.est = z[,keeper.features[[i]]],
                plot=FALSE, residual.plot = FALSE, order=depth, n.best=n.best,
                norm="std", factor.2.dummy = FALSE, ncores=subcores)$Point.est
          }
      } else {
          for(i in 1:length(keeper.features)){

              if(status){
                  message("% of Final Estimate = ", format(i/length(keeper.features),digits =  3,nsmall = 2),"     ","\r",appendLF=FALSE)
                  if(i == length(keeper.features)){
                      message("% of Final Estimate = 1.000             ","\r",appendLF=TRUE)
                      flush.console()
                  }
              }

              estimates[[i]] <- NNS.reg(x[,keeper.features[[i]]],y,point.est = z[,keeper.features[[i]]],
                          plot=FALSE, residual.plot = FALSE, order=depth, n.best=n.best,
                          norm="std", factor.2.dummy = FALSE, ncores=subcores)$Point.est
          }

      }

      if(!is.null(cl)){
          stopCluster(cl)
          registerDoSEQ()
      }


      estimates <- lapply(estimates, function(i) pmin(i,max(as.numeric(y))))
      estimates <- lapply(estimates, function(i) pmax(i,min(as.numeric(y))))


      plot.table <- table(unlist(keeper.features))
      names(plot.table) <- colnames(x[as.numeric(names(plot.table))])

      if(feature.importance){

          linch <-  max(strwidth(names(plot.table), "inch")+0.4, na.rm = TRUE)
          par(mai=c(1.0,linch,0.8,0.5))

          if(length(plot.table)!=1){
              barplot(sort(plot.table,decreasing = FALSE)[1:min(n,10)],
                    horiz = TRUE,
                    col='steelblue',
                    main="Feature Frequency in Final Estimate",
                    xlab = "Frequency",las=1)
          } else {
              barplot(sort(plot.table,decreasing = FALSE)[1:min(n,10)],
                    horiz = TRUE,
                    col='steelblue',
                    main="Feature Frequency in Final Estimate",
                    xlab = "Frequency",las=1)
      }
     par(mar=c(5.1, 4.1, 4.1, 2.1))
     par(original.par)
  }
  gc()
  return(list("results"=apply(do.call(cbind,estimates),1,mode),
              "feature.weights"=plot.table/sum(plot.table)))
}
