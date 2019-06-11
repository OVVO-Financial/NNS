#' NNS Stack
#'
#' Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
#'
#' @param IVs.train a vector, matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param IVs.test a vector, matrix or data frame of variables of numeric or factor data types with compatible dimsensions to \code{(IVs.train)}.
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to 0.25 for a 25 percent random sampling of the training set under \code{(CV.size = NULL)}.
#' @param folds integer; \code{folds = 5} (default) Select the number of cross-validation folds.
#' @param order integer; \code{NULL} (default) Sets the order for \link{NNS.reg}, where \code{(order = 'max')} is the k-nearest neighbors equivalent.
#' @param norm options: ("std", "NNS", NULL); \code{NULL} (default) 3 settings offered: \code{NULL}, \code{"std"}, and \code{"NNS"}.  Selects the \code{norm} parameter in \link{NNS.reg}.
#' @param method numeric options: (1, 2); Select the NNS method to include in stack.  \code{(method = 1)} selects \link{NNS.reg}; \code{(method = 2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method = c(1, 2)}, including both NNS regression methods in the stack.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all") method for determining synthetic X* coefficients.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for weights.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.reg}. If NULL (default), the number of cores to be used is equal to half the number of cores of the machine - 1.
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
#' @note If character variables are used, transform them first to factors using \link{as.factor}, or \link{data.matrix} to ensure overall dataset is numeric.  A multifunction \link{sapply} can also be applied to the overall dataset: \code{data <- sapply(data,function(x){as.factor(x) ; as.numeric(x)})}.  Then run \code{NNS.stack} with transormed variables.
#'
#' Missing data should be handled prior as well using \link{na.omit} or \link{complete.cases} on the full dataset.
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
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4])}
#'
#'  ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
#'  \dontrun{
#'  NNS.stack(iris[ , 1:4], iris[ , 5])}
#'
#'  ## Selecting NNS.reg and dimension reduction techniques.
#'  \dontrun{
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], iris[141:150, 1:4], method = c(1, 2))}
#'
#' @export

NNS.stack <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = 'min',
                      CV.size = NULL,
                      folds = 5,
                      order = NULL,
                      norm = NULL,
                      method = c(1, 2),
                      dim.red.method = "cor",
                      ncores = NULL){


  DV.train <- as.numeric(DV.train)


  n <- ncol(IVs.train)

  l <- length(IVs.train[ , 1])

  if(is.null(CV.size)){
      if(is.null(IVs.test)){
          CV.size <- 0.25
      } else {
          CV.size <- mean(c(.2, min(length(IVs.test[ , 1]) / l, .5)))
      }
  }

  THRESHOLDS <- list()
  best.k <- list()
  best.nns.cv <- list()
  best.nns.ord <- list()

  for(b in 1 : folds){
      set.seed(123 * b)
      test.set <- sample(1 : l, as.integer(CV.size * l), replace = FALSE)

      if(b > 1){
          test.set_half <- unique(c(rbind(test.set.1,test.set.2)))[1:(length(test.set)/2)]
          test.set <- unique(c(test.set_half,sample(1 : l, replace = FALSE)))[1:length(test.set)]
      }

    CV.IVs.train <- IVs.train[c(-test.set), ]
    CV.IVs.test <- IVs.train[c(test.set), ]
    CV.DV.train <- DV.train[c(-test.set)]
    CV.DV.test <- DV.train[c(test.set)]

    training <- cbind(IVs.train[c(-test.set),],DV.train[c(-test.set)])
    training <- training[complete.cases(training),]

    CV.IVs.train <- training[,-(ncol(training))]
    CV.DV.train <- training[,ncol(training)]


    ### NORMALIZATION OF VARIABLES and SELECTION OF ORDER:
    if(!is.null(norm)){
        np <- nrow(CV.IVs.test)
        points.norm <- rbind(CV.IVs.test, CV.IVs.train)
        colnames(points.norm) <- colnames(CV.IVs.test)

        if(norm == "std"){
            CV.IVs.train <- apply(CV.IVs.train, 2, function(b) (b - min(b)) / (max(b) - min(b)))
            CV.IVs.test <- apply(points.norm, 2, function(b) (b - min(b)) / (max(b) - min(b)))[1 : np, ]
        }

        if(norm == "NNS"){
            CV.IVs.train <- NNS.norm(CV.IVs.train)
            CV.IVs.test <- NNS.norm(points.norm)[1 : np, ]
        }
    }

    if(1 %in% method){
        actual <- CV.DV.test
        nns.cv.1 <- numeric()

        if(objective=='min'){
            nns.cv.1[1] <- Inf
        } else {
            nns.cv.1[1] <- -Inf
        }

        for(i in 1:(2*n)){
            predicted <- NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = i, order=order, ncores = ncores)$Point.est

            nns.cv.1[i+1] <- eval(obj.fn)
            if(i > 1){
                if(objective=='min' && nns.cv.1[i]>nns.cv.1[i-1]){ break }
                if(objective=='max' && nns.cv.1[i]<nns.cv.1[i-1]){ break }
            }
        }

        test.set.1 <- test.set[rev(order(abs(predicted - actual)))]

        if(objective=='min'){
            k <- which.min(na.omit(nns.cv.1))-1
            nns.cv.1 <- min(na.omit((nns.cv.1)))
        } else {
            k <- which.max(na.omit(nns.cv.1))-1
            nns.cv.1 <- max(na.omit((nns.cv.1)))
        }


      best.k[[b]] <- k
      best.nns.cv[[b]] <- nns.cv.1

      if(b==folds){
          best.nns.cv <- mean(na.omit(unlist(best.nns.cv)))
          best.k <- as.numeric(names(sort(table(unlist(best.k)),decreasing = TRUE)[1]))

          nns.method.1 <- NNS.reg(IVs.train, DV.train, point.est = IVs.test, plot = FALSE, residual.plot = FALSE, n.best = best.k, order=order, ncores = ncores)$Point.est
      }

    } else {
        test.set.1 <- NULL
        best.k <- NA
        nns.method.1 <- NA
        if(objective=='min'){best.nns.cv <- Inf} else {best.nns.cv <- -Inf}
    }# 1 %in% method



    # Dimension Reduction Regression Output
    if(2 %in% method){
        actual <- CV.DV.test

        var.cutoffs <- abs(round(NNS.reg(CV.IVs.train, CV.DV.train, dim.red.method = dim.red.method, plot = FALSE, residual.plot = FALSE, order=order, ncores = ncores)$equation$Coefficient, digits = 2))

        var.cutoffs <- var.cutoffs - .005

        var.cutoffs <- var.cutoffs[var.cutoffs <= 1 & var.cutoffs >= 0]

        var.cutoffs <- rev(sort(unique(var.cutoffs)))

        nns.ord <- numeric()

        if(objective=='min'){nns.ord[1] <- Inf} else {nns.ord[1] <- -Inf}

        for(i in 2:length(var.cutoffs)){
            predicted <- NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, dim.red.method = dim.red.method, threshold = var.cutoffs[i], order=order, ncores = ncores)$Point.est

            nns.ord[i] <- eval(obj.fn)

            if(objective=='min'){
                best.threshold <- var.cutoffs[which.min(na.omit(nns.ord))]
                THRESHOLDS[[b]] <- best.threshold
                best.nns.ord[[b]] <- min(na.omit(nns.ord))
                if(nns.ord[i] > nns.ord[i-1]) break
                } else {
                    best.threshold <- var.cutoffs[which.max(na.omit(nns.ord))]
                    THRESHOLDS[[b]] <- best.threshold
                    best.nns.ord[[b]] <- max(na.omit(nns.ord))
                    if(nns.ord[i] < nns.ord[i-1]) break
                }
        }


        test.set.2 <- test.set[rev(order(abs(predicted - actual)))]

        if(b==folds){
            nns.ord.threshold <- as.numeric(names(sort(table(unlist(THRESHOLDS)),decreasing = TRUE)[1]))
            best.nns.ord <- mean(na.omit(unlist(best.nns.ord)))
            nns.method.2 <- NNS.reg(IVs.train, DV.train,point.est = IVs.test, dim.red.method = dim.red.method, plot = FALSE, order=order,
                               threshold = nns.ord.threshold, ncores = ncores)$Point.est
        }

    } else {
        THRESHOLDS <- NA
        test.set.2 <- NULL
        nns.method.2 <- NA
        if(objective=='min'){best.nns.ord <- Inf} else {best.nns.ord <- -Inf}
        nns.ord.threshold <- NA
    } # 2 %in% method


  } # errors (b) loop


  ### Weights for combining NNS techniques
  best.nns.cv[best.nns.cv == 0] <- 1e-10
  best.nns.ord[best.nns.ord == 0] <- 1e-10

  if(objective=="min"){
      weights <- c(max(1e-10, 1 / best.nns.cv), max(1e-10, 1 / best.nns.ord))
  } else {
      weights <- c(max(1e-10,best.nns.cv), max(1e-10, best.nns.ord))
  }

  weights <- pmax(weights, c(0, 0))
  weights[!(c(1, 2) %in% method)] <- 0
  weights <- weights / sum(weights)

  if(identical(sort(method),c(1,2))){
      estimates <- (weights[1] * nns.method.1 + weights[2] * nns.method.2)
  } else {
      if(method==1){
          estimates <- (nns.method.1)
      } else {
          if(method==2) {
              estimates <- (nns.method.2)
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
