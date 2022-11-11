#' NNS Stack
#'
#' Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
#'
#' @param IVs.train a vector, matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimensions to \code{(IVs.train)}.
#' @param IVs.test a vector, matrix or data frame of variables of numeric or factor data types with compatible dimensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.
#' @param type \code{NULL} (default).  To perform a classification of discrete integer classes from factor target variable \code{(DV.train)} with a base category of 1, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.   Like a logistic regression, this setting is not necessary for target variable of two classes e.g. [0, 1].
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param optimize.threshold logical; \code{TRUE} (default) Will optimize the probability threshold value for rounding in classification problems.  If \code{FALSE}, returns 0.5.
#' @param dist options:("L1", "L2", "DTW", "FACTOR") the method of distance calculation; Selects the distance calculation used. \code{dist = "L2"} (default) selects the Euclidean distance and \code{(dist = "L1")} selects the Manhattan distance; \code{(dist = "DTW")} selects the dynamic time warping distance; \code{(dist = "FACTOR")} uses a frequency.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to a random value between 0.2 and 0.33 for a random sampling of the training set.
#' @param balance logical; \code{FALSE} (default) Uses both up and down sampling from \code{caret} to balance the classes.  \code{type="CLASS"} required.
#' @param ts.test integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.
#' @param folds integer; \code{folds = 5} (default) Select the number of cross-validation folds.
#' @param order options: (integer, "max", NULL); \code{NULL} (default) Sets the order for \link{NNS.reg}, where \code{(order = "max")} is the k-nearest neighbors equivalent, which is suggested for mixed continuous and discrete (unordered, ordered) data.
#' @param norm options: ("std", "NNS", NULL); \code{NULL} (default) 3 settings offered: \code{NULL}, \code{"std"}, and \code{"NNS"}.  Selects the \code{norm} parameter in \link{NNS.reg}.
#' @param method numeric options: (1, 2); Select the NNS method to include in stack.  \code{(method = 1)} selects \link{NNS.reg}; \code{(method = 2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method = c(1, 2)}, which will reduce the dimension first, then find the optimal \code{n.best}.
#' @param stack logical; \code{TRUE} (default) Uses dimension reduction output in \code{n.best} optimization, otherwise performs both analyses independently.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "equal", "all") method for determining synthetic X* coefficients.  \code{(dim.red.method = "cor")} uses standard linear correlation for weights.  \code{(dim.red.method = "NNS.dep")} (default) uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.reg}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns a vector of fitted values for the dependent variable test set for all models.
#' \itemize{
#' \item{\code{"NNS.reg.n.best"}} returns the optimum \code{"n.best"} parameter for the \link{NNS.reg} multivariate regression.  \code{"SSE.reg"} returns the SSE for the \link{NNS.reg} multivariate regression.
#' \item{\code{"OBJfn.reg"}} returns the \code{obj.fn} for the \link{NNS.reg} regression.
#' \item{\code{"NNS.dim.red.threshold"}} returns the optimum \code{"threshold"} from the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"OBJfn.dim.red"}} returns the \code{obj.fn} for the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"probability.threshold"}} returns the optimum probability threshold for classification, else 0.5 when set to \code{FALSE}.
#' \item{\code{"reg"}} returns \link{NNS.reg} output.
#' \item{\code{"dim.red"}} returns \link{NNS.reg} dimension reduction regression output.
#' \item{\code{"stack"}} returns the output of the stacked model.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#' \url{https://www.ssrn.com/abstract=2864711}
#'
#' @note
#' \itemize{
#' \item Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].  The response variable base category should be 1 for multiple class problems.
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
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], iris[141:150, 1:4], method = c(1, 2), type = "CLASS")
#'  }
#' @export

NNS.stack <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      type = NULL,
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = "min",
                      optimize.threshold = TRUE,
                      dist = "L2",
                      CV.size = NULL,
                      balance = FALSE,
                      ts.test = NULL,
                      folds = 5,
                      order = NULL,
                      norm = NULL,
                      method = c(1, 2),
                      stack = TRUE,
                      dim.red.method = "cor",
                      status = TRUE,
                      ncores = NULL){
  
  if(sum(is.na(cbind(IVs.train,DV.train))) > 0) stop("You have some missing values, please address.")
  
  if(is.null(obj.fn)) stop("Please provide an objective function")
  
  if(balance && is.null(type)) warning("type = 'CLASS' selected due to balance = TRUE.")
  if(balance) type <- "CLASS"
  
  if(!is.null(type) && min(as.numeric(DV.train))==0) warning("Base response variable category should be 1, not 0.")
  
  if(any(class(IVs.train)%in%c("tbl","data.table"))) IVs.train <- as.data.frame(IVs.train)
  if(any(class(DV.train)%in%c("tbl","data.table"))) DV.train <- as.vector(unlist(DV.train))
  
  if(is.null(dim(IVs.train)) || dim(IVs.train)[2]==1){
    IVs.train <- data.frame(IVs.train)
    method <- 1
    order <- NULL
  }
  
  
  if(!is.null(type)){
    type <- tolower(type)
    if(type == "class" && identical(obj.fn,expression( sum((predicted - actual)^2) ))){
      obj.fn <- expression(mean( predicted == as.numeric(actual)))
      objective <- "max"
    }
  }
  
  objective <- tolower(objective)
  
  if(!is.null(type) && type=="class") DV.train <- as.numeric(factor(DV.train)) else DV.train <- as.numeric(DV.train)
  
  n <- ncol(IVs.train)
  
  l <- floor(sqrt(length(IVs.train[ , 1])))
  
  if(is.null(IVs.test)){
    IVs.test <- IVs.train
  } else {
    if(any(class(IVs.test)%in%c("tbl","data.table"))) IVs.test <- as.data.frame(IVs.test)
  }
  
  if(is.null(dim(IVs.test))) IVs.test <- data.frame(t(IVs.test)) else IVs.test <- data.frame(IVs.test)
  
  dist <- tolower(dist)
  
  THRESHOLDS <- vector(mode = "list", folds)
  best.k <- vector(mode = "list", folds)
  best.nns.cv <- vector(mode = "list", folds)
  best.nns.ord <- vector(mode = "list", folds)
  
  if(is.null(colnames(IVs.train))){
    colnames.list <- lapply(1 : dim(IVs.train)[2], function(i) paste0("X", i))
    colnames(IVs.test) <- colnames(IVs.train) <- as.character(colnames.list)
  }
  
  if(2 %in% method && dim(IVs.train)[2]>1){
    if(dim.red.method=="cor"){
      var.cutoffs_1 <- abs(round(cor(data.matrix(cbind(DV.train, IVs.train)), method = "spearman")[-1,1], digits = 2))
    } else {
      var.cutoffs_1 <- abs(round(suppressWarnings(NNS.reg(IVs.train, DV.train, dim.red.method = dim.red.method, plot = FALSE, residual.plot = FALSE, order=order, ncores = ncores,
                                                          type = type, point.only = TRUE)$equation$Coefficient[-(n+1)]), digits = 2))
    }
  }
  
  if(is.null(CV.size)) new.CV.size <- round(runif(1, .2, 1/3), 3) else new.CV.size <- CV.size
  
  for(b in 1 : folds){
    if(status) message("Folds Remaining = " , folds-b," ","\r",appendLF=TRUE)
    
    set.seed(123 * b)
    
    test.set <- as.integer(seq(b, length(unlist(IVs.train[ , 1])), length.out = as.integer(new.CV.size * length(unlist(IVs.train[ , 1])))))
    
    if(!is.null(ts.test)){
      test.set <- 1:(length(DV.train) - ts.test)
    }
    
    test.set <- unlist(test.set)
    
    CV.IVs.train <- data.frame(IVs.train[c(-test.set), ])
    
    if(dim(CV.IVs.train)[2]!=dim(IVs.train)[2]) CV.IVs.train <- t(CV.IVs.train)
    if(dim(CV.IVs.train)[2]!=dim(IVs.train)[2]) CV.IVs.train <- t(CV.IVs.train)
    
    
    CV.IVs.test <- data.frame(IVs.train[test.set, ])
    if(dim(CV.IVs.test)[2]!=dim(IVs.train)[2]) CV.IVs.test <- t(CV.IVs.test)
    if(dim(CV.IVs.test)[2]!=dim(IVs.train)[2]) CV.IVs.test <- t(CV.IVs.test)
    
    CV.DV.train <- DV.train[c(-test.set)]
    CV.DV.test <- DV.train[c(test.set)]
    
    training <- cbind(IVs.train[c(-test.set),], DV.train[c(-test.set)])
    training <- training[complete.cases(training),]
    
    if(balance){
      DV.train <- as.numeric(as.factor(DV.train))

      CV.DV.train <- DV.train[c(-test.set)]
      CV.DV.test <- DV.train[c(test.set)]
      
      y_train <- as.factor(CV.DV.train)
      training_1 <- do.call(cbind, caret::downSample(CV.IVs.train, y_train, list = TRUE))
      training_2 <- do.call(cbind, caret::upSample(CV.IVs.train, y_train, list = TRUE))
      
      training <- rbind.data.frame(training_1, training_2)
      
      colnames(training) <- c(colnames(CV.IVs.train), names(CV.DV.train))
    }
    
    
    CV.IVs.train <- data.frame(training[, -(ncol(training))])
    
    CV.DV.train <- as.numeric(as.character(training[, ncol(training)]))
    
    
    # Dimension Reduction Regression Output
    if(2 %in% method && dim(IVs.train)[2]>1){
      
      actual <- CV.DV.test
      if(dim.red.method=="cor"){
        
        var.cutoffs_2 <- abs(round(cor(data.matrix(cbind(CV.DV.train, CV.IVs.train)), method = "spearman")[-1,1], digits = 2))
      } else {
        
        var.cutoffs_2 <- abs(round(suppressWarnings(NNS.reg(CV.IVs.train, CV.DV.train, dim.red.method = dim.red.method, plot = FALSE, residual.plot = FALSE, order=order, ncores = ncores,
                                                            type = type, point.only = TRUE)$equation$Coefficient[-(n+1)]), digits = 2))
      }
      
      var.cutoffs <- c(pmin(var.cutoffs_1, (pmax(var.cutoffs_1, var.cutoffs_2) + pmin(var.cutoffs_1, var.cutoffs_2))/2))
      
      var.cutoffs <- var.cutoffs[var.cutoffs < 1 & var.cutoffs > 0]
      
      var.cutoffs[is.na(var.cutoffs)] <- 0
      
      var.cutoffs <- rev(sort(unique(var.cutoffs)))[-1]
      
      if(is.null(var.cutoffs)) var.cutoffs <- 0
      
      if(n == 2) var.cutoffs <- c(var.cutoffs, 0)
      
      if(dist=="factor") var.cutoffs <- var.cutoffs[-1]
      if(dim.red.method=="equal") var.cutoffs <- 0
      
      threshold_results_2 <- vector(mode = "list", length = length(var.cutoffs))
      nns.ord <- numeric(length(var.cutoffs))
      
      for(i in 1:length(var.cutoffs)){
        if(status){
          message("Current NNS.reg(... , threshold = ", var.cutoffs[i] ," ) MAX Iterations Remaining = " ,length(var.cutoffs)-i," ","\r",appendLF=TRUE)
        }
        
        predicted <- suppressWarnings(NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, dim.red.method = dim.red.method, threshold = var.cutoffs[i], order = NULL, ncores = ncores,
                                              type = NULL, dist = dist, point.only = TRUE)$Point.est)
        
        predicted[is.na(predicted)] <- gravity(na.omit(predicted))
        
        if(!is.null(type)){
          pred_matrix <- sapply(seq(.01, .99, .01), function(z) ifelse(predicted%%1<z, as.integer(floor(predicted)), as.integer(ceiling(predicted))))
          z <- apply(pred_matrix, 2, function(z) mean(z == as.numeric(actual)))
          threshold_results_2[[i]] <- seq(.01,.99, .01)[as.integer(median(which(z==max(z))))]
          
          predicted <- ifelse(predicted%%1 < threshold_results_2[[i]], floor(predicted), ceiling(predicted))
        }
        
        nns.ord[i] <- eval(obj.fn)
        
        if(objective=="min"){
          best.threshold <- gravity(var.cutoffs[na.omit(nns.ord)==min(na.omit(nns.ord))])
          THRESHOLDS[[b]] <- best.threshold
          best.nns.ord[[b]] <- min(na.omit(nns.ord))
          if(is.na(nns.ord[1])) nns.ord[1] <- Inf
          if(i > 2 && is.na(nns.ord[i])) break
          if(i > 2 && (nns.ord[i] >= nns.ord[i-1]) && (nns.ord[i] >= nns.ord[i-2])) break
        } else {
          best.threshold <- gravity(var.cutoffs[na.omit(nns.ord)==max(na.omit(nns.ord))])
          THRESHOLDS[[b]] <- best.threshold
          best.nns.ord[[b]] <- max(na.omit(nns.ord))
          if(is.na(nns.ord[1])) nns.ord[1] <- -Inf
          if(i > 2 && is.na(nns.ord[i])) break
          if(i > 2 && (nns.ord[i] <= nns.ord[i-1]) && (nns.ord[i] <= nns.ord[i-2])) break
        }
      }
      
      
      relevant_vars <- colnames(IVs.train)
      if(is.null(relevant_vars)) relevant_vars <- 1:dim(IVs.train)[2]
      
      if(b==folds){
        threshold.table <- sort(table(unlist(THRESHOLDS)), decreasing = TRUE)
        
        nns.ord.threshold <- gravity(as.numeric(names(threshold.table[threshold.table==max(threshold.table)])))
        
        nns.method.2 <- suppressWarnings(NNS.reg(IVs.train, DV.train, point.est = IVs.test, dim.red.method = dim.red.method, plot = FALSE, order = order, threshold = nns.ord.threshold, ncores = ncores,
                                                 type = NULL, point.only = TRUE))
        
        actual <- nns.method.2$Fitted.xy$y
        predicted <- nns.method.2$Fitted.xy$y.hat
        
        best.nns.ord <- eval(obj.fn)
        
        rel_vars <- nns.method.2$equation
        
        rel_vars <- which(rel_vars$Coefficient>0)
        rel_vars <- rel_vars[rel_vars <= dim(IVs.train)[2]]
        
        if(is.null(rel_vars) || length(rel_vars)==0) rel_vars <- 1:dim(IVs.train)[2]
        
        if(!stack) relevant_vars <- 1:dim(IVs.train)[2]
        
        if(all(relevant_vars=="FALSE")){
          relevant_vars <- 1:dim(IVs.train)[2]
        }
        
        if(!is.null(type) && !is.null(nns.method.2$Point.est)){
          threshold_results_2 <- mean(unlist(threshold_results_2))
          
          nns.method.2 <- ifelse(nns.method.2$Point.est%%1 < threshold_results_2, floor(nns.method.2$Point.est), ceiling(nns.method.2$Point.est))
          nns.method.2 <- pmin(nns.method.2, max(as.numeric(DV.train)))
          nns.method.2 <- pmax(nns.method.2, min(as.numeric(DV.train)))
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
      threshold_results_2 <- NA
      relevant_vars <- 1:dim(IVs.train)[2]
    } # 2 %in% method
    
    
    
    if(1 %in% method){
      actual <- CV.DV.test
      
      if(is.character(relevant_vars)) relevant_vars <- relevant_vars!=""
      
      if(is.logical(relevant_vars)){
        CV.IVs.train <- data.frame(CV.IVs.train[, relevant_vars])
        CV.IVs.test <- data.frame(CV.IVs.test[, relevant_vars])
      }
      
      
      if(dim(CV.IVs.train)[2]!=dim(IVs.train)[2]) CV.IVs.train <- t(CV.IVs.train)
      if(dim(CV.IVs.train)[2]!=dim(IVs.train)[2]) CV.IVs.train <- t(CV.IVs.train)
      
      if(dim(CV.IVs.test)[2]!=dim(IVs.train)[2]) CV.IVs.test <- t(CV.IVs.test)
      if(dim(CV.IVs.test)[2]!=dim(IVs.train)[2]) CV.IVs.test <- t(CV.IVs.test)
      
      threshold_results_1 <- vector(mode = "list", length(c(1:l, length(IVs.train[ , 1]))))
      nns.cv.1 <- numeric()
      
      for(i in c(1:l, length(IVs.train[ , 1]))){
        index <- which(c(1:l, length(IVs.train[ , 1])) == i)
        if(status){
          message("Current NNS.reg(... , n.best = ", i ," ) MAX Iterations Remaining = " ,l-index+1," ","\r",appendLF=TRUE)
        }
        
        if(index==1){
          setup <- suppressWarnings(NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = 1, order = order,
                                            type = type, factor.2.dummy = TRUE, dist = dist, ncores = ncores, point.only = TRUE))
          
          nns.id <- setup$Fitted.xy$NNS.ID
          original.DV <- setup$Fitted.xy$y
          
          predicted <- setup$Point.est
          
          predicted[is.na(predicted)] <- mean(predicted, na.rm = TRUE)
          pred_matrix <- sapply(seq(.01, .99, .01), function(z) ifelse(predicted%%1<z, as.integer(floor(predicted)), as.integer(ceiling(predicted))))
          
          threshold_results_1[index] <- seq(.01,.99, .01)[which.max(apply(pred_matrix, 2, function(z) mean(z == as.numeric(actual))))]
          
          predicted <- ifelse(predicted%%1 < threshold_results_1[index], floor(predicted), ceiling(predicted))
          
          RPM_CLASS <- apply(do.call(cbind, lapply(setup$RPM[,1:(dim(setup$RPM)[2]-1)], FUN = function(z) ifelse(z%%1 < .5, floor(z), ceiling(z)))), 2, as.integer)
        } else {
          
          
          if(!is.null(dim(CV.IVs.train))){
            if(dim(CV.IVs.train)[2]>1){
              CV.IVs.test.new <- data.table::data.table(do.call(cbind, lapply(data.frame(CV.IVs.test), factor_2_dummy_FR)))
              
              CV.IVs.test.new <- CV.IVs.test.new[, DISTANCES :=  NNS.distance(rpm = setup$RPM, rpm_class = RPM_CLASS, dist.estimate = .SD, type = dist, k = i, class = type)[1], by = 1:nrow(CV.IVs.test)]
              
              predicted <- as.numeric(unlist(CV.IVs.test.new$DISTANCES))
            } else {
              predicted <-  suppressWarnings(NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = i, order = order, ncores = ncores,
                                                     type = type, factor.2.dummy = TRUE, dist = dist, point.only = TRUE)$Point.est)
            }
            
            rm(CV.IVs.test.new)
          } else {
            predicted <-  suppressWarnings(NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, residual.plot = FALSE, n.best = i, order = order, ncores = ncores,
                                                   type = type, factor.2.dummy = TRUE, dist = dist, point.only = TRUE)$Point.est)
          }
          
          if(!is.null(type)){
            pred_matrix <- sapply(seq(.01, .99, .01), function(z) ifelse(predicted%%1<z, as.integer(floor(predicted)), as.integer(ceiling(predicted))))
            z <- apply(pred_matrix, 2, function(z) mean(z == as.numeric(actual)))
            threshold_results_1[[index]] <- seq(.01,.99, .01)[as.integer(median(which(z==max(z))))]
            
            predicted <- ifelse(predicted%%1 < threshold_results_1[[index]], floor(predicted), ceiling(predicted))
          }
          
          
        }
        
        nns.cv.1[index] <- eval(obj.fn)
        
        if(length(na.omit(nns.cv.1)) > 3){
          if(objective=="min") nns.cv.1[is.na(nns.cv.1)] <- max(na.omit(nns.cv.1)) else nns.cv.1[is.na(nns.cv.1)] <- min(na.omit(nns.cv.1))
          if(objective=='min' && nns.cv.1[index]>=nns.cv.1[index-1] && nns.cv.1[index]>=nns.cv.1[index-2]){ break }
          if(objective=='max' && nns.cv.1[index]<=nns.cv.1[index-1] && nns.cv.1[index]<=nns.cv.1[index-2]){ break }
        }
      }
      
      
      ks <- c(1:l, length(IVs.train[ , 1]))[!is.na(nns.cv.1)]
      
      if(objective=='min'){
        k <- ks[which.min(na.omit(nns.cv.1))]
        nns.cv.1 <- min(na.omit(nns.cv.1))
      } else {
        k <- ks[which.max(na.omit(nns.cv.1))]
        nns.cv.1 <- max(na.omit(nns.cv.1))
      }
      
      
      best.k[[b]] <- k
      best.nns.cv[[b]] <- if(!is.null(type)) min(max(nns.cv.1,0),1) else nns.cv.1
      
      if(b==folds){
        ks <- table(unlist(best.k))
        
        best.k <-  mode_class(as.numeric(rep(names(ks), as.numeric(unlist(ks)))))
        
        nns.method.1 <- suppressWarnings(NNS.reg(IVs.train[ , relevant_vars], DV.train, point.est = IVs.test[, relevant_vars], plot = FALSE, n.best = best.k, order = order, ncores = ncores,
                                                 type = NULL, point.only = FALSE))
        
        actual <- nns.method.1$Fitted.xy$y
        predicted <- nns.method.1$Fitted.xy$y.hat
        
        best.nns.cv <- eval(obj.fn)
        
        nns.method.1 <- nns.method.1$Point.est
        
        if(!is.null(type) && !is.null(nns.method.1)){
          threshold_results_1 <- mean(unlist(threshold_results_1))
          nns.method.1 <- ifelse(nns.method.1%%1 < threshold_results_1, floor(nns.method.1), ceiling(nns.method.1))
          nns.method.1 <- pmin(nns.method.1, max(as.numeric(DV.train)))
          nns.method.1 <- pmax(nns.method.1, min(as.numeric(DV.train)))
        }
        
      }
      
      
    } else {
      test.set.1 <- NULL
      best.k <- NA
      nns.method.1 <- NA
      threshold_results_1 <- NA
      if(objective=='min'){best.nns.cv <- Inf} else {best.nns.cv <- -Inf}
    }# 1 %in% method
    
    
    
    
  } # errors (b) loop
  
  
  ### Weights for combining NNS techniques
  best.nns.cv[best.nns.cv == 0] <- 1e-10
  best.nns.ord[best.nns.ord == 0] <- 1e-10
  
  if(objective=="min"){
    weights <- c(max(1e-10, 1 / best.nns.cv^2), max(1e-10, 1 / best.nns.ord^2))
  } else {
    weights <- c(max(1e-10, best.nns.cv^2), max(1e-10, best.nns.ord^2))
  }
  
  
  weights <- pmax(weights, c(0, 0))
  weights[!(c(1, 2) %in% method)] <- 0
  weights[is.nan(weights)] <- 0
  weights[is.infinite(weights)] <- 0
  
  if(sum(weights)>0)  weights <- weights / sum(weights) else weights <- c(.5, .5)
  
  if(!is.null(type)) probability.threshold <-  mean(c(threshold_results_1, threshold_results_2), na.rm = TRUE) else probability.threshold <- .5
  
  if(identical(sort(method),c(1,2))){
    if(sum(is.na(nns.method.1)>0)){
      na.1.index <- which(is.na(nns.method.1))
      nns.method.1[na.1.index] <- nns.method.2[na.1.index]
    }
    if(sum(is.na(nns.method.2)>0)){
      na.2.index <- which(is.na(nns.method.2))
      nns.method.2[na.2.index] <- nns.method.1[na.2.index]
    }
    
    estimates <- (weights[1] * nns.method.1 + weights[2] * nns.method.2)
    if(!is.null(type)){
      estimates <- ifelse(estimates%%1 < probability.threshold, floor(estimates), ceiling(estimates))
      estimates <- pmin(estimates, max(as.numeric(DV.train)))
      estimates <- pmax(estimates, min(as.numeric(DV.train)))
    }
  } else {
    if(method==1){
      estimates <- nns.method.1
    } else {
      if(method==2) estimates <- nns.method.2
    }
  }
  
  
  if(is.null(probability.threshold)) probability.threshold <- .5
  
  return(list(OBJfn.reg = best.nns.cv,
              NNS.reg.n.best = best.k,
              probability.threshold = probability.threshold,
              OBJfn.dim.red = best.nns.ord,
              NNS.dim.red.threshold = nns.ord.threshold,
              reg = nns.method.1,
              dim.red = nns.method.2,
              stack = estimates))
  
}