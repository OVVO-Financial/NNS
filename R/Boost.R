#' NNS Boost
#'
#' Ensemble method for classification using the predictions of the NNS multivariate regression \link{NNS.reg} collected from uncorrelated feature combinations.
#'
#' @param IVs.train a matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param IVs.test a matrix or data frame of variables of numeric or factor data types with compatible dimsensions to \code{(IVs.train)}.
#' @param representative.sample logical; \code{TRUE} (default) Reduces observations of \code{IVs.train} to a set of representative observations per regressor.
#' @param depth integer; \code{"max"} (default) Specifies the \code{order} parameter in the \code{NNS.reg} routine, assigning a number of splits in the regressors.
#' @param n.best integer; \code{3} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at \code{sqrt(# of regressors)}. Analogous to \code{k} in a \code{k Nearest Neighbors} algorithm.
#' @param learner.trials integer; \code{NULL} (default) Sets the number of trials to obtain an accuracy \code{threshold} level.  Number of observations in the training set is the default setting.
#' @param epochs integer; \code{2*length(DV.train)} (default) Total number of feature combinations to run.
#' @param CV.size numeric [0, 1]; \code{(CV.size = .25)} (default) Sets the cross-validation size.  Defaults to 0.25 for a 25 percent random sampling of the training set.
#' @param threshold numeric [0, 1]; \code{NULL} (default) Sets the \code{obj.fn} accuracy threshold to keep feature combinations.
#' @param extreme logical; \code{FALSE} (default) Uses the maximum \code{threshold} obtained from the \code{learner.trials}, rather than the upper quintile level.
#' @param feature.importance logical; \code{TRUE} (default) Plots the frequency of features used in the final estimate.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @return Returns a vector of fitted values for the dependent variable test set.
#'
#' @keywords classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#' \url{https://ssrn.com/abstract=2864711}
#' @note Currently only objective function is to maximize \code{mean(round(predicted)==as.numeric(actual))} since it is a classifier.  For a regression problem, simply use \code{NNS.reg}!
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  a = NNS.boost(iris[1:140, 1:4], iris[1:140, 5],
#'  IVs.test = iris[141:150, 1:4],
#'  epochs = 100, learner.trials = 100)
#'
#'  ## Test accuracy
#'  mean(round(a)==as.numeric(iris[141:150,5]))
#'  }
#'
#' @export


NNS.boost <- function(IVs.train,
                      DV.train,
                      IVs.test,
                      representative.sample = TRUE,
                      depth = "max",
                      n.best = 3,
                      learner.trials = NULL,
                      epochs = NULL,
                      CV.size=.25,
                      threshold = NULL,
                      extreme = FALSE,
                      feature.importance = TRUE,
                      status = TRUE,
                      ncores = NULL){


  obj.fn = expression(mean(round(predicted)==as.numeric(actual)))
  objective = "max"

  if (is.null(ncores)) {
    num_cores <- detectCores() - 1
  } else {
    num_cores <- ncores
  }


  mode = function(x){
    if(length(na.omit(x)) > 1){
      d <- density(na.omit(x))
      d$x[which.max(d$y)]
    } else {
      x
    }
  }

  x = data.frame(IVs.train); y = DV.train; z = data.frame(IVs.test)

    factor_2_dummy = function(x){
      if(class(x) == "factor"){
        n=length(unique(x))
        output = model.matrix(~x -1, x)[,-1]
      } else {
        output = x
      }
      output
    }

      x = apply(as.data.frame(x),2,factor_2_dummy)
      x = do.call(cbind, as.data.frame(x))
      x = as.data.frame(x)

      z = apply(as.data.frame(z),2,factor_2_dummy)
      z = do.call(cbind, as.data.frame(z))
      z = as.data.frame(z)



    ### Representative samples
      rep.x = data.matrix(x)

      rep.x = data.table(rep.x)

      fivenum.x = rep.x[,lapply(.SD,fivenum), by = .(y)]
      mode.x = rep.x[,lapply(.SD,mode), by = .(y)]
      mean.x = rep.x[,lapply(.SD,mean), by = .(y)]

      rep.x = rbind(fivenum.x,mode.x,mean.x)

      rep.y = unlist(rep.x[,1])
      rep.x = as.data.frame(rep.x[,-1])


  n = ncol(x)

  if(is.null(epochs)){
      epochs = 2*length(y)
  }

  estimates=list()
  fold = list()

  old.threshold = 0


  # Add test loop for highest threshold ...
  if(is.null(threshold)){
      old.threshold = 1
      test.features = list()
      results = numeric()

      if(is.null(learner.trials)){learner.trials = length(y)}

      for(i in 1:learner.trials){
          new.index = sample(length(y), as.integer(CV.size*length(y)), replace = FALSE)

          new.iv.train = data.table(x[-new.index,])

          fivenum.new.iv.train = new.iv.train[,lapply(.SD,fivenum), by = .(y[-new.index])]
          mode.new.iv.train = new.iv.train[,lapply(.SD,mode), by = .(y[-new.index])]
          mean.new.iv.train = new.iv.train[,lapply(.SD,mean), by = .(y[-new.index])]

          new.iv.train = rbind(fivenum.new.iv.train,mode.new.iv.train,mean.new.iv.train)

              new.iv.train = as.data.frame(new.iv.train[,-1])
              new.dv.train = unlist(new.iv.train[,1])

          if(!representative.sample){
              new.iv.train = rbind(new.iv.train,x[-new.index,])
              new.dv.train = c(new.dv.train,y[-new.index])
          }

        actual = y[new.index]
        new.iv.test = x[new.index,]

        if(status){
          message("Current Threshold Iterations Remaining = " ,learner.trials+1-i," ","\r",appendLF=FALSE)
        }

          test.features[[i]] = sort(sample(n,sample(2:n,1),replace = FALSE))

      #If estimate is > threshold, store 'features'
      predicted = NNS.reg(new.iv.train[,test.features[[i]]],new.dv.train,point.est = new.iv.test[,test.features[[i]]],plot=FALSE,residual.plot = FALSE,order=depth,n.best=n.best,norm="std")$Point.est

      predicted = pmin(predicted,max(as.numeric(y)))
      predicted = pmax(predicted,min(as.numeric(y)))

      results[i] = eval(obj.fn)

    } # i in learner.trials
  } # NULL thresholde

    if(feature.importance){
          original.par = par()
          par(mfrow=c(2,1))
          par(mai=c(1.0,.5,0.8,0.5))
          hist(results,main = "Distribution of Learner Trials Accuracy",
               xlab = "Accuracy",col = "steelblue")
    }

    if(extreme){
        threshold = max(results)
    } else {
        threshold = fivenum(results)[4]
    }

    if(status){
    message(paste0("Learner Accuracy Threshold = ", format(threshold,digits = 3,nsmall = 2),"           "),appendLF = TRUE)

    # Clear message line
    message("                                       ","\r",appendLF=FALSE)
    }


      keeper.features = list()


      for(j in 1:epochs){

      new.index = sample(1:length(y),as.integer(CV.size*length(y)),replace = FALSE)

      new.iv.train = data.table(x[-new.index,])

      fivenum.new.iv.train = new.iv.train[,lapply(.SD,fivenum), by = .(y[-new.index])]
      mode.new.iv.train = new.iv.train[,lapply(.SD,mode), by = .(y[-new.index])]
      mean.new.iv.train = new.iv.train[,lapply(.SD,mean), by = .(y[-new.index])]

      new.iv.train = rbind(fivenum.new.iv.train,mode.new.iv.train,mean.new.iv.train)


        new.iv.train = as.data.frame(new.iv.train[,-1])
        new.dv.train = unlist(new.iv.train[,1])

        if(!representative.sample){
        new.iv.train = rbind(new.iv.train,x[-new.index,])
        new.dv.train = c(new.dv.train,y[-new.index])
      }

      actual = y[new.index]
      new.iv.test = x[new.index,]

          if(status){
          message("% of epochs = ", format(j/epochs,digits =  3,nsmall = 2),"     ","\r",appendLF=FALSE)

          if(j == epochs){
              message("% of epochs ",j," = 1.000     ","\r",appendLF=FALSE)
              flush.console()
          }
          }

          features = sort(sample(n,sample(2:n,1),replace = FALSE))

          #If estimate is > threshold, store 'features'
          predicted = NNS.reg(new.iv.train[,features],new.dv.train,point.est = new.iv.test[,features],plot=FALSE,residual.plot = FALSE,order=depth,n.best=n.best,norm="std")$Point.est

          predicted = pmin(predicted,max(as.numeric(y)))
          predicted = pmax(predicted,min(as.numeric(y)))

          new.results = eval(obj.fn)

          if(new.results>=threshold){
              keeper.features[[j]]=features
          } else {
              keeper.features[[j]]=NULL
          }

      }

      keeper.features = keeper.features[!sapply(keeper.features, is.null)]
      if(length(keeper.features)==0){
          if(old.threshold==0){
                stop("Please reduce [threshold].")
          } else {
                keeper.features = test.features[which.max(results)]
          }
      }

      for(i in 1:length(keeper.features)){
      if(status){
      message(paste0("% of Final Estimate  = ", format(i/length(keeper.features),digits = 3,nsmall = 2),"     "),"\r",appendLF=FALSE)

      if(i == length(keeper.features)){
          message("% of Final Estimate  = 1.000     ","\r",appendLF=FALSE)
          flush.console()
      }
      }
    x = rbind(rep.x,x); y = c(rep.y,y)

      estimates[[i]]= NNS.reg(x[,keeper.features[[i]]],y,point.est = z[,keeper.features[[i]]],plot=FALSE,residual.plot = FALSE,order=depth,n.best=n.best,norm="std")$Point.est



  }

      estimates = lapply(estimates, function(i) pmin(i,max(as.numeric(y))))
      estimates = lapply(estimates, function(i) pmax(i,min(as.numeric(y))))

  if(feature.importance==TRUE){
      plot.table = table(unlist(keeper.features))
      names(plot.table)=names(IVs.train[as.numeric(names(plot.table))])

      linch <-  max(strwidth(names(plot.table), "inch")+0.4, na.rm = TRUE)
      par(mai=c(1.0,linch,0.8,0.5))

      if(length(plot.table)!=1){
          barplot(rev(sort(plot.table,decreasing = TRUE)),
                horiz = TRUE,
                col='steelblue',
                main="Feature Importance in Final Estimate",
                xlab = "Frequency",las=1)
      } else {
          barplot(plot.table,
                horiz = TRUE,
                col='steelblue',
                main="Feature Importance in Final Estimate",
                xlab = "Frequency",las=1)
      }

    par(original.par)
  }

  return(apply(do.call(cbind,estimates),1,mode))
}
