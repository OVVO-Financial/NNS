#' NNS Boost
#'
#' Ensemble method using the predictions of the NNS base models \link{NNS.reg} for various uncorrelated feature combinations.
#'
#' @param IVs.train a matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param IVs.test a matrix or data frame of variables of numeric or factor data types with compatible dimsensions to \code{(IVs.train)}.
#' @param depth integer; \code{NULL} (default) Specifies the \code{order} parameter in the \code{NNS.reg} routine, assigning a number of splits in the regressors.
#' @param epochs integer; 500 (default) Total number of feature combinations to run.
#' @param folds integer; 5 (default) Number of times to resample the training data.  Splits the \code{epochs} over the dataset evenly over each \code{folds}.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to 0.25 for a 25 percent random sampling of the training set under \code{(CV.size = NULL)}.
#' @param threshold numeric [0, 1]; \code{NULL} (default) Sets the \code{obj.fn} threshold to keep feature combinations.
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
#'  a = NNS.boost(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4], epochs = 100)
#'
#'  ## Test accuracy
#'  mean(round(a)==as.numeric(iris[141:150,5]))
#'  }
#'
#' @export


NNS.boost <- function(IVs.train,
                      DV.train,
                      IVs.test,
                      depth = NULL,
                      epochs=500,
                      folds=5,
                      CV.size=.2,
                      threshold=NULL,
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



  x=data.frame(IVs.train); y=DV.train; z=data.frame(IVs.test)

  n = ncol(x)

  epochs = min(epochs, sum(choose(n,1:n)), length(y))

  estimates=list()
  fold = list()

  # Test sample for threshold
  new.index = sample(length(x[,1]), as.integer(CV.size*length(x[,1])), replace = FALSE)

  actual = y[new.index]
  new.iv.test = x[new.index,]
  new.iv.train = x[-new.index,]
  new.dv.train = y[-new.index]


  # Add test loop for highest threshold
  if(is.null(threshold)){
    results = numeric()
    for(i in rep(seq(.99,0,-.01),each=1)){
      message("Current Threshold = ",i,"\r",appendLF=FALSE)
      features= sample(ncol(x),sample(2:ncol(x),1),replace = FALSE)

      #If estimate is > threshold, store 'features'
      predicted = NNS.reg(new.iv.train[,features],new.dv.train,point.est = new.iv.test[,features],plot=FALSE,residual.plot = FALSE,order=depth)$Point.est
      results[which(i%in%rep(seq(.99,0,-.01),each=1))] = eval(obj.fn)

      if(results[which(i%in%rep(seq(.99,0,-.01),each=1))]>=i){
        threshold = max(results,i)
        break
      }
    }
  }
  message("                           ","\r",appendLF=FALSE)
  fold = list()
  for(i in 1:folds){
    keeper.features = list()

    new.index = sample(1:length(x[,1]),as.integer(CV.size*length(x[,1])),replace = FALSE)

    #actual = y[new.index]
    new.iv.test = x[new.index,]
    new.iv.train = x[-new.index,]
    new.dv.train = y[-new.index]


    for(j in 1:as.integer(epochs/folds)){
      message("% of Fold ",i," = ", format(j/as.integer(epochs/folds),digits =  3,nsmall = 2),"     ","\r",appendLF=FALSE)
      if(j == as.integer(epochs/folds)){
        message("% of Fold ",i," = 1.000     ","\r",appendLF=FALSE)
        flush.console()
      }
      actual = y[new.index]
      features = sample(ncol(x),sample(2:ncol(x),1),replace = FALSE)

      if(i>1 && (list(features)%in%fold)){next}

      #If estimate is > threshold, store 'features'
      predicted = NNS.reg(new.iv.train[,features],new.dv.train,point.est = new.iv.test[,features],plot=FALSE,residual.plot = FALSE,order=depth)$Point.est

      results = eval(obj.fn)
      if(results>threshold){keeper.features[[j]]=features[order(features)]} else {NULL}
    }

    keeper.features = keeper.features[!sapply(keeper.features, is.null)]
    keeper.features = unique(keeper.features)
    fold[[i]]= keeper.features


  }

  fold = fold[!sapply(fold, is.null)]
  fold = unique(fold)

  if(length(fold)==0) stop("Please reduce [threshold]")

  final.features = do.call(c,fold)
  final.features = unique(final.features)

  if(length(final.features)==0) stop("Please reduce [threshold]")

  estimates = list()


  for(i in 1:length(final.features)){
    message(paste0("% of Final Estimate  = ", format(i/length(final.features),digits = 3,nsmall = 2),"     "),"\r",appendLF=FALSE)
    if(i == length(final.features)){
      message("% of Final Estimate  = 1.000     ","\r",appendLF=FALSE)
      flush.console()
    }

    estimates[[i]]= NNS.reg(x[,final.features[[i]]],y,point.est = z[,final.features[[i]]],plot=FALSE,residual.plot = FALSE,order=depth)$Point.est

  }


  return(apply(do.call(cbind,estimates),1,mode))
}
