#' NNS Stack
#'
#' Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
#'
#' @param IVs.train a vector, matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param IVs.test a vector, matrix or data frame of variables of numeric or factor data types with compatible dimsensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.
#' @param type \code{NULL} (default).  To perform a classification of disrete integer classes from factor target variable \code{(DV.train)}, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.   Like a logistic regression, this setting is not necessary for target variable of two classes e.g. [0, 1].
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param dist options:("L1", "L2", "DTW", "FACTOR") the method of distance calculation; Selects the distance calculation used. \code{dist = "L2"} (default) selects the Euclidean distance and \code{(dist = "L1")} seclects the Manhattan distance; \code{(dist = "DTW")} selects the dynamic time warping distance; \code{(dist = "FACTOR")} uses a frequency.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to 0.25 for a 25 percent random sampling of the training set under \code{(CV.size = NULL)}.
#' @param ts.test integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.
#' @param folds integer; \code{folds = 5} (default) Select the number of cross-validation folds.
#' @param order options: (integer, "max", NULL); \code{NULL} (default) Sets the order for \link{NNS.reg}, where \code{(order = "max")} is the k-nearest neighbors equivalent, which is suggested for mixed continuous and discrete (unordered, ordered) data.
#' @param norm options: ("std", "NNS", NULL); \code{NULL} (default) 3 settings offered: \code{NULL}, \code{"std"}, and \code{"NNS"}.  Selects the \code{norm} parameter in \link{NNS.reg}.
#' @param method numeric options: (1, 2); Select the NNS method to include in stack.  \code{(method = 1)} selects \link{NNS.reg}; \code{(method = 2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method = c(1, 2)}, which will reduce the dimension first, then find the optimal \code{n.best}.
#' @param stack logical; \code{TRUE} (default) Uses dimension reduction output in \code{n.best} optimization, otherwise performs both analyses independently.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all") method for determining synthetic X* coefficients.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for weights.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.reg}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns a vector of fitted values for the dependent variable test set for all models.
#' \itemize{
#' \item{\code{"NNS.reg.n.best"}} returns the optimum \code{"n.best"} paramater for the \link{NNS.reg} multivariate regression.  \code{"SSE.reg"} returns the SSE for the \link{NNS.reg} multivariate regression.
#' \item{\code{"OBJfn.reg"}} returns the \code{obj.fn} for the \link{NNS.reg} regression.
#' \item{\code{"NNS.dim.red.threshold"}} returns the optimum \code{"threshold"} from the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"OBJfn.dim.red"}} returns the \code{obj.fn} for the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"reg"}} returns \link{NNS.reg} output.
#' \item{\code{"dim.red"}} returns \link{NNS.reg} dimension reduction regression output.
#' \item{\code{"stack"}} returns the output of the stacked model.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#' \url{https://ssrn.com/abstract=2864711}
#'
#' @note
#' \itemize{
#' \item Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].
#'
#' \item Missing data should be handled prior as well using \link{na.omit} or \link{complete.cases} on the full dataset.
#' }
#'
#' If error received:
#'
#' \code{"Error in is.data.frame(x) : object 'RP' not found"}
#'
#' reduce the \code{CV.size}.
#'
#'
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4], type = "CLASS")
#'
#'  ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
#'  NNS.stack(iris[ , 1:4], iris[ , 5], type = "CLASS")
#'
#'  ## Selecting NNS.reg and dimension reduction techniques.
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], iris[141:150, 1:4], method = c(1, 2), type = "CLASS")}
#'
#' @export

NNS.stack <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      type = NULL,
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = "min",
                      dist = "L2",
                      CV.size = NULL,
                      ts.test = NULL,
                      folds = 5,
                      order = NULL,
                      norm = NULL,
                      method = c(1, 2),
                      stack = TRUE,
                      dim.red.method = "cor",
                      status = TRUE,
                      ncores = NULL){


  if(is.null(obj.fn)){ stop("Please provide an objective function")}

  if(is.null(dim(IVs.train)) || dim(IVs.train)[2]==1){
    IVs.train <- data.frame(IVs.train)
    method <- 1
    order <- NULL
  }



  if(!is.null(type)){
    type <- tolower(type)
    if(type=="class"){
      obj.fn <- expression(mean( predicted == as.numeric(actual)))
      objective <- "max"
    }
  }

  if(dist == "FACTOR"){
    method <- 1
  }

  objective <- tolower(objective)

  DV.train <- as.numeric(DV.train)

  n <- ncol(IVs.train)

  l <- floor(sqrt(length(IVs.train[ , 1])))

  if(is.null(IVs.test)){
    IVs.test <- IVs.train
  } else {
    if(is.null(dim(IVs.test)) || dim(IVs.test)[2]==1){
      IVs.test <- data.frame(IVs.test)
    }
  }


  if(is.null(CV.size)){
    if(is.null(IVs.test)){
      CV.size <- 0.25
    } else {
      CV.size <- mean(c(.2, min(length(IVs.test[ , 1]) / l, .5)))
    }
  }

  THRESHOLDS <- list(folds)
  best.k <- list(folds)
  best.nns.cv <- list(folds)
  best.nns.ord <- list(folds)

  if (is.null(ncores)) {
    num_cores <- detectCores() - 1
  } else {
    num_cores <- ncores
  }

  for(b in 1 : folds){
    if(status){
      message("Folds Remaining = " , folds-b," ","\r",appendLF=TRUE)
    }

    set.seed(123 * b)
    test.set <- sample(1 : length(IVs.train[ , 1]), as.integer(CV.size * length(IVs.train[ , 1])), replace = FALSE)

    if(b > 1){
      maxes <- as.vector(apply(IVs.train, 2, which.max))
      mins <- as.vector(apply(IVs.train, 2, which.min))
      test.set_half <- unique(c(rbind(test.set.1, test.set.2)))[1:(length(test.set)/2)]

      test.set <- unique(c(mins, maxes, test.set_half, sample(1 : length(IVs.train[ , 1]), replace = FALSE)))[1:length(test.set)]
      test.set <- na.omit(test.set)
    }

    if(!is.null(ts.test)){
      test.set <- 1:(length(DV.train) - ts.test)
      dist <- "DTW"
    }

    test.set <- unlist(test.set)

    CV.IVs.train <- data.frame(IVs.train[c(-test.set), ])
    CV.IVs.test <- data.frame(IVs.train[test.set, ])

    CV.DV.train <- DV.train[c(-test.set)]
    CV.DV.test <- DV.train[c(test.set)]

    training <- cbind(IVs.train[c(-test.set),], DV.train[c(-test.set)])
    training <- training[complete.cases(training),]

    CV.IVs.train <- data.frame(training[, -(ncol(training))])
    CV.DV.train <- training[, ncol(training)]


    # Dimension Reduction Regression Output
    if(2 %in% method && dim(IVs.train)[2]>1){

      actual <- CV.DV.test

      var.cutoffs_1 <- abs(round(NNS.reg(IVs.train, DV.train, dim.red.method = dim.red.method, plot = FALSE, residual.plot = FALSE, order=order, ncores = 1,
                                         type = type)$equation$Coefficient, digits = 2))

      var.cutoffs_2 <- abs(round(NNS.reg(CV.IVs.train, CV.DV.train, dim.red.method = dim.red.method, plot = FALSE, residual.plot = FALSE, order=order, ncores = 1,
                                         type = type)$equation$Coefficient, digits = 2))

      var.cutoffs <- (pmax(var.cutoffs_1, var.cutoffs_2) + pmin(var.cutoffs_1, var.cutoffs_2))/2

      var.cutoffs <- var.cutoffs[var.cutoffs < 1]

      var.cutoffs[is.na(var.cutoffs)] <- 0

      var.cutoffs <- rev(sort(unique(var.cutoffs)))[-1]

      if(is.null(var.cutoffs)) var.cutoffs <- 0

      if(n == 2) var.cutoffs <- c(var.cutoffs, 0)

      if(dist=="FACTOR"){
        var.cutoffs <- var.cutoffs[-1]
      }
      nns.ord <- numeric()

      for(i in 1:length(var.cutoffs)){
        if(status){
          message("Current NNS.reg(... , threshold = ", var.cutoffs[i] ," ) MAX Iterations Remaining = " ,length(var.cutoffs)-i," ","\r",appendLF=TRUE)
        }

        predicted <- NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, dim.red.method = dim.red.method, threshold = var.cutoffs[i], order = NULL, ncores = 1,
                             type = type, dist = dist)$Point.est

        nns.ord[i] <- eval(obj.fn)

        if(objective=="min"){
          best.threshold <- var.cutoffs[which.min(na.omit(nns.ord))]
          THRESHOLDS[[b]] <- best.threshold
          best.nns.ord[[b]] <- min(na.omit(nns.ord))
          if(i > 2 && is.na(nns.ord[i])) break
          if(i > 2 && (nns.ord[i] >= nns.ord[i-1])) break
        } else {
          best.threshold <- var.cutoffs[which.max(na.omit(nns.ord))]
          THRESHOLDS[[b]] <- best.threshold
          best.nns.ord[[b]] <- max(na.omit(nns.ord))
          if(i > 2 && is.na(nns.ord[i])) break
          if(i > 2  && (nns.ord[i] <= nns.ord[i-1])) break
        }
      }

      test.set.2 <- test.set[rev(order(abs(predicted - actual)))]

      relevant_vars <- colnames(IVs.train)

      if(b==folds){
        nns.ord.threshold <- as.numeric(names(sort(table(unlist(THRESHOLDS)), decreasing = TRUE)[1]))
        best.nns.ord <- mean(na.omit(unlist(best.nns.ord)))
        nns.method.2 <- NNS.reg(IVs.train, DV.train,point.est = IVs.test, dim.red.method = dim.red.method, plot = FALSE, order = order, threshold = nns.ord.threshold, ncores = 1,
                                type = type)

        rel_vars <- nns.method.2$equation
        rel_vars <- rel_vars[rel_vars$Coefficient>0,1][-.N]


        if(stack){
            relevant_vars <- colnames(IVs.train)%in%unlist(rel_vars)
        } else {
            relevant_vars <- colnames(IVs.train)%in%colnames(IVs.train)
        }

        if(all(relevant_vars=="FALSE")){
            relevant_vars <- colnames(IVs.train)%in%colnames(IVs.train)
        }

        if(!is.null(type) && !is.null(nns.method.2$Point.est)){
          nns.method.2 <- round(nns.method.2$Point.est)
        } else {
          nns.method.2 <- nns.method.2$Point.est
        }


      }

    } else {
      THRESHOLDS <- NA
      test.set.2 <- NULL
      nns.method.2 <- NA
      if(objective=='min'){best.nns.ord <- Inf} else {best.nns.ord <- -Inf}
      nns.ord.threshold <- NA
      relevant_vars <- colnames(IVs.train)%in%colnames(IVs.train)
    } # 2 %in% method



    if(1 %in% method){
      actual <- CV.DV.test
      nns.cv.1 <- numeric()
      if(length(relevant_vars)<=1 || is.logical(relevant_vars)){
          relevant_vars <- seq_len(dim(IVs.train)[2])
      }

      CV.IVs.train <- data.frame(CV.IVs.train[, relevant_vars])
      CV.IVs.test <- data.frame(CV.IVs.test[, relevant_vars])

      for(i in c(1:l, length(IVs.train[ , 1]))){
        index <- which(c(1:l, length(IVs.train[ , 1])) %in% i)
        if(status){
          message("Current NNS.reg(... , n.best = ", i ," ) MAX Iterations Remaining = " ,l-index+1," ","\r",appendLF=TRUE)
        }

        if(index==1){
          setup <- NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = i, order = order, ncores = 1,
                           type = type, factor.2.dummy = TRUE, dist = dist)
          predicted <- setup$Point.est
        } else {

          predicted <- list()


          if(!is.null(dim(CV.IVs.train))){
            if(dim(CV.IVs.train)[2]>1){
                cl <- makeCluster(num_cores)
                registerDoParallel(cl)

                CV.IVs.test <- do.call(cbind, lapply(data.frame(CV.IVs.test), factor_2_dummy_FR))

                predicted <- foreach(j = 1:nrow(CV.IVs.test), .packages=c("NNS","data.table","dtw"))%dopar%{
                    NNS.distance(setup$RPM, dist.estimate = as.vector(CV.IVs.test[j,]), type = dist, k = i)
                }

                stopCluster(cl)
                registerDoSEQ()
            } else {
              predicted <-  NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = i, order = order, ncores = 1,
                                    type = type, factor.2.dummy = TRUE, dist = dist)$Point.est
            }
          } else {
            predicted <-  NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = i, order = order, ncores = 1,
                                  type = type, factor.2.dummy = TRUE, dist = dist)$Point.est
          }
          predicted <- unlist(predicted)

          if(!is.null(type)){
            predicted <- round(predicted)
          }
        }

        nns.cv.1[index] <- eval(obj.fn)

        if(length(na.omit(nns.cv.1)) > 2){
          if(objective=='min' && nns.cv.1[index]>=nns.cv.1[index-1] && nns.cv.1[index]>=nns.cv.1[index-2]){ break }
          if(objective=='max' && nns.cv.1[index]<=nns.cv.1[index-1] && nns.cv.1[index]<=nns.cv.1[index-2]){ break }
        }
      }

      if(length(predicted > 0)){
        test.set.1 <- test.set[rev(order(abs(predicted - actual)))]
      } else {
        test.set.1 <- test.set
      }

      if(objective=='min'){
        ks <- c(1:l, length(IVs.train[ , 1]))[!is.na(nns.cv.1)]
        k <- ks[which.min(na.omit(nns.cv.1))]
        nns.cv.1 <- min(na.omit(nns.cv.1))
      } else {
        ks <- c(1:l, length(IVs.train[ , 1]))[!is.na(nns.cv.1)]
        k <- ks[which.max(na.omit(nns.cv.1))]
        nns.cv.1 <- max(na.omit(nns.cv.1))
      }


      best.k[[b]] <- k
      best.nns.cv[[b]] <- nns.cv.1

      if(b==folds){
        best.nns.cv <- mean(na.omit(unlist(best.nns.cv)))
        best.k <- round(fivenum(as.numeric(rep(names(table(unlist(best.k))), table(unlist(best.k)))))[4])
        nns.method.1 <- NNS.reg(IVs.train[ , relevant_vars], DV.train, point.est = IVs.test[, relevant_vars], plot = FALSE, n.best = best.k, order = order, ncores = 1,
                                type = type)$Point.est
        if(!is.null(type) && !is.null(nns.method.1)){
          nns.method.1 <- round(nns.method.1)
        }

      }


    } else {
      test.set.1 <- NULL
      best.k <- NA
      nns.method.1 <- NA
      if(objective=='min'){best.nns.cv <- Inf} else {best.nns.cv <- -Inf}
    }# 1 %in% method




  } # errors (b) loop


  ### Weights for combining NNS techniques
  best.nns.cv[best.nns.cv == 0] <- 1e-10
  best.nns.ord[best.nns.ord == 0] <- 1e-10

  if(objective=="min"){
    weights <- c(max(1e-10, 1 / best.nns.cv), max(1e-10, 1 / best.nns.ord))
  } else {
    weights <- c(max(1e-10, best.nns.cv), max(1e-10, best.nns.ord))
  }

  weights <- pmax(weights, c(0, 0))
  weights[!(c(1, 2) %in% method)] <- 0
  weights <- weights / sum(weights)

  if(identical(sort(method),c(1,2))){
    estimates <- (weights[1] * nns.method.1 + weights[2] * nns.method.2)
    if(!is.null(type)){
      estimates <- round(estimates)
    }
  } else {
    if(method==1){
      estimates <- nns.method.1
    } else {
      if(method==2) {
        estimates <- nns.method.2
      }
    }
  }

  gc()

  return(list(OBJfn.reg = best.nns.cv,
              NNS.reg.n.best = best.k,
              OBJfn.dim.red = best.nns.ord,
              NNS.dim.red.threshold = nns.ord.threshold,
              reg = nns.method.1,
              dim.red = nns.method.2,
              stack = estimates))

}
