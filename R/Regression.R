#' NNS Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x a vector, matrix or data frame of variables of numeric or factor data types.
#' @param y a numeric or factor vector with compatible dimensions to \code{x}.
#' @param factor.2.dummy logical; \code{TRUE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.
#' @param order integer; Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{(order = ...)} integer settings with \code{(noise.reduction = "off")}.  \code{(order = "max")} will force a limit condition perfect fit.
#' @param stn numeric [0, 1]; Signal to noise parameter, sets the threshold of \code{(NNS.dep)} which reduces \code{("order")} when \code{(order = NULL)}.  Defaults to 0.95 to ensure high dependence for higher \code{("order")} and endpoint determination.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all", "equal", \code{numeric vector}, NULL) method for determining synthetic X* coefficients.  Selection of a method automatically engages the dimension reduction regression.  The default is \code{NULL} for full multivariate regression.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "cor")} uses standard linear correlation for weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.  \code{(dim.red.method = "equal")} uses unit weights.  Alternatively, user can specify a numeric vector of coefficients.
#' @param tau options("ts", NULL); \code{NULL}(default) To be used in conjunction with \code{(dim.red.method = "NNS.caus")} or \code{(dim.red.method = "all")}.  If the regression is using time-series data, set \code{(tau = "ts")} for more accurate causal analysis.
#' @param type \code{NULL} (default).  To perform a classification, set to \code{(type = "CLASS")}.  Like a logistic regression, it is not necessary for target variable of two classes e.g. [0, 1].
#' @param point.est a numeric or factor vector with compatible dimensions to \code{x}.  Returns the fitted value \code{y.hat} for any value of \code{x}.
#' @param location Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.
#' @param return.values logical; \code{TRUE} (default), set to \code{FALSE} in order to only display a regression plot and call values as needed.
#' @param plot logical; \code{TRUE} (default) To plot regression.
#' @param plot.regions logical; \code{FALSE} (default).  Generates 3d regions associated with each regression point for multivariate regressions.  Note, adds significant time to routine.
#' @param residual.plot logical; \code{TRUE} (default) To plot \code{y.hat} and \code{Y}.
#' @param confidence.interval numeric [0, 1]; \code{NULL} (default) Plots the associated confidence interval with the estimate and reports the standard error for each individual segment.  Also applies the same level for the prediction intervals.
#' @param threshold  numeric [0, 1]; \code{(threshold = 0)} (default) Sets the threshold for dimension reduction of independent variables when \code{(dim.red.method)} is not \code{NULL}.
#' @param n.best integer; \code{NULL} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at \code{sqrt(# of regressors)}.  \code{(n.best = "all")} will select and weight all generated regression points.  Analogous to \code{k} in a
#' \code{k Nearest Neighbors} algorithm.  Different values of \code{n.best} are tested using cross-validation in \link{NNS.stack}.
#' @param noise.reduction the method of determining regression points options: ("mean", "median", "mode", "off"); In low signal:noise situations,\code{(noise.reduction = "mean")}  uses means for \link{NNS.dep} restricted partitions, \code{(noise.reduction = "median")} uses medians instead of means for \link{NNS.dep} restricted partitions, while \code{(noise.reduction = "mode")}  uses modes instead of means for \link{NNS.dep} restricted partitions.  \code{(noise.reduction = "off")} uses an overall central tendency measure for partitions.
#' @param dist options:("L1", "L2", "FACTOR") the method of distance calculation; Selects the distance calculation used. \code{dist = "L2"} (default) selects the Euclidean distance and \code{(dist = "L1")} selects the Manhattan distance; \code{(dist = "FACTOR")} uses a frequency.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param multivariate.call Internal argument for multivariate regressions.
#' @param point.only Internal argument for abbreviated output.
#' @return UNIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:
#' \itemize{
#'  \item{\code{"R2"}} provides the goodness of fit;
#'
#'  \item{\code{"SE"}} returns the overall standard error of the estimate between \code{y} and \code{y.hat};
#'
#'  \item{\code{"Prediction.Accuracy"}} returns the correct rounded \code{"Point.est"} used in classifications versus the categorical \code{y};
#'
#'  \item{\code{"derivative"}} for the coefficient of the \code{x} and its applicable range;
#'
#'  \item{\code{"Point.est"}} for the predicted value generated;
#'  
#'  \item{\code{"pred.int"}} lower and upper prediction intervals for the \code{"Point.est"} returned using the \code{"confidence.interval"} provided;
#'  
#'  \item{\code{"regression.points"}} provides the points used in the regression equation for the given order of partitions;
#'
#'  \item{\code{"Fitted.xy"}} returns a \link{data.table} of \code{x}, \code{y}, \code{y.hat}, \code{resid}, \code{NNS.ID}, \code{gradient};
#' }
#'
#'
#' MULTIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:
#' \itemize{
#'  \item{\code{"R2"}} provides the goodness of fit;
#'
#'  \item{\code{"equation"}} returns the numerator of the synthetic X* dimension reduction equation as a \link{data.table} consisting of regressor and its coefficient.  Denominator is simply the length of all coefficients > 0, returned in last row of \code{equation} data.table.
#'
#'  \item{\code{"x.star"}} returns the synthetic X* as a vector;
#'
#'  \item{\code{"rhs.partitions"}} returns the partition points for each regressor \code{x};
#'
#'  \item{\code{"RPM"}} provides the Regression Point Matrix, the points for each \code{x} used in the regression equation for the given order of partitions;
#'
#'  \item{\code{"Point.est"}} returns the predicted value generated;
#'  
#'  \item{\code{"pred.int"}} lower and upper prediction intervals for the \code{"Point.est"} returned using the \code{"confidence.interval"} provided;
#'
#'  \item{\code{"Fitted.xy"}} returns a \link{data.table} of \code{x},\code{y}, \code{y.hat}, \code{gradient}, and \code{NNS.ID}.
#' }
#'
#' @note
#' \itemize{
#'  \item Please ensure \code{point.est} is of compatible dimensions to \code{x}, error message will ensue if not compatible.
#'
#'  \item Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].  The response variable base category should be 1 for classification problems.
#'
#'  \item For low signal:noise instances, increasing the dimension may yield better results using \code{NNS.stack(cbind(x,x), y, method = 1, ...)}.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
#' \url{https://link.springer.com/article/10.1007/s10614-017-9713-5}
#'
#' Vinod, H. and Viole, F. (2018) "Clustering and Curve Fitting by Line Segments"
#' \url{https://www.preprints.org/manuscript/201801.0090/v1}
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.reg(x, y)
#'
#' ## Manual {order} selection
#' NNS.reg(x, y, order = 2)
#'
#' ## Maximum {order} selection
#' NNS.reg(x, y, order = "max")
#'
#' ## x-only paritioning (Univariate only)
#' NNS.reg(x, y, type = "XONLY")
#'
#' ## For Multiple Regression:
#' x <- cbind(rnorm(100), rnorm(100), rnorm(100)) ; y <- rnorm(100)
#' NNS.reg(x, y, point.est = c(.25, .5, .75))
#'
#' ## For Multiple Regression based on Synthetic X* (Dimension Reduction):
#' x <- cbind(rnorm(100), rnorm(100), rnorm(100)) ; y <- rnorm(100)
#' NNS.reg(x, y, point.est = c(.25, .5, .75), dim.red.method = "cor", ncores = 1)
#'
#' ## IRIS dataset examples:
#' # Dimension Reduction:
#' NNS.reg(iris[,1:4], iris[,5], dim.red.method = "cor", order = 5, ncores = 1)
#'
#' # Dimension Reduction using causal weights:
#' NNS.reg(iris[,1:4], iris[,5], dim.red.method = "NNS.caus", order = 5, ncores = 1)
#'
#' # Multiple Regression:
#' NNS.reg(iris[,1:4], iris[,5], order = 2, noise.reduction = "off")
#'
#' # Classification:
#' NNS.reg(iris[,1:4], iris[,5], point.est = iris[1:10, 1:4], type = "CLASS")$Point.est
#'
#' ## To call fitted values:
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.reg(x, y)$Fitted
#'
#' ## To call partial derivative (univariate regression only):
#' NNS.reg(x, y)$derivative
#' }
#' @export


NNS.reg = function (x, y,
                    factor.2.dummy = TRUE, order = NULL,
                    stn = .95,
                    dim.red.method = NULL, tau = NULL,
                    type = NULL,
                    point.est = NULL,
                    location = "top",
                    return.values = TRUE,
                    plot = TRUE, plot.regions = FALSE, residual.plot = TRUE,
                    confidence.interval = NULL,
                    threshold = 0,
                    n.best = NULL,
                    noise.reduction = "off",
                    dist = "L2",
                    ncores = NULL,
                    point.only = FALSE,
                    multivariate.call = FALSE){
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")
  
  if(plot.regions && !is.null(order) && order == "max") stop('Please reduce the "order" or set "plot.regions = FALSE".')
  
  dist <- tolower(dist)
  
  if(any(class(x)%in%c("tbl","data.table")) && dim(x)[2]==1) x <- as.vector(unlist(x))
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
 
  n <- length(y)
  original.x <- x
  
  if(n < 2000) ncores <- 1
  if(!is.null(dim(x))){
    if(ncol(x) < 5) ncores <- 1
  }
  
  if(!is.null(dim.red.method)){
    if(is.null(dim(x)) || dim(x)[1]==1){
      dim.red.method <- NULL
    }
  }
  
  synthetic.x.equation <- NULL
  x.star <- NULL
  
  if(!is.null(type)){
    type <- tolower(type)
    noise.reduction <- "mode_class"
  }
  
  if(is.discrete(y) && length(unique(y)) < sqrt(length(y))){
    type <- "class"
    noise.reduction <- "mode_class"
  }
  
  if(any(class(y)==c("tbl", "data.table"))) y <- as.vector(unlist(y))
  
  if(!plot) residual.plot <- FALSE
  
  # Variable names
  original.names <- colnames(x)
  original.columns <- ncol(x)
  
  if(!is.null(original.columns) & is.null(colnames(x))) x <- data.frame(x)
 
  y.label <- deparse(substitute(y))
  if(is.null(y.label)) y.label <- "y"
  
  
  if(factor.2.dummy){
    if(is.list(x) & !is.data.frame(x)) x <- do.call(cbind, x)
    
    
    if(!is.null(point.est)){
      if(!is.null(dim(x)) && dim(x)[2]>1){
        if(is.null(dim(point.est))) point.est <- data.frame(t(point.est)) else point.est <- data.frame(point.est)
        new_x <- data.table::rbindlist(list(data.frame(x), point.est), use.names = FALSE)
      } else {
        new_x <- unlist(list(x, point.est))
      }
    }
 
    
    if(!is.null(dim(x)) && dim(x)[2] > 1) x <- apply(x, 2, function(z) factor_2_dummy_FR(z)) else x <- factor_2_dummy_FR(x)
       
    
    x <- data.matrix(x)
    
    if(!is.null(point.est)){
      point.est.y <- numeric()
      
      if(!is.null(dim(x)) && dim(x)[2]>1){
        new_x <- apply(new_x, 2, function(z) factor_2_dummy_FR(z))
      } else {
        new_x <- factor_2_dummy_FR(new_x)
      }
      
      if(is.null(dim(point.est))) l_point.est <- length(point.est) else l_point.est <- dim(point.est)[1]
      
      
      point.est <- tail(new_x, l_point.est)
      
      if(is.null(dim(point.est)) || dim(point.est)[2]==1) point.est <- as.vector(unlist(point.est))
      
    } else { # is.null(point.est)
      point.est.y <- NULL
    }
  } #if(factor.2.dummy)
  
  # Variable names
  original.names <- colnames(x)
  original.columns <- ncol(x)
  
  y <- as.numeric(y)
  original.y <- y
  
  
  if(!factor.2.dummy){
    if(is.null(ncol(x))){
      x <- as.double(x)
      if(!is.null(point.est)){
        point.est <- as.double(point.est)
        point.est.y <- numeric()
      } else {
        point.est.y <- NULL
      }
    } else {
      x <- data.matrix(x)
      if(!is.null(point.est)){
        if(is.null(ncol(point.est))){
          point.est <- as.double(point.est)
          point.est.y <- numeric()
        } else {
          point.est <- data.matrix(point.est)
          point.est.y <- numeric()
        }
      } else {
        point.est.y <- NULL
      }
    }
  } # !factor to dummy
  original.variable <- x
  
  np <- nrow(point.est)
  
  if(!is.null(type) && type == "class" ){
    if(is.null(n.best)) n.best <- 1
  }
  
  
  if(!is.null(ncol(original.variable))){
    if(ncol(original.variable) == 1){
      x <- original.variable
    } else {
      if(is.null(dim.red.method)){
        
        colnames(x) <- make.unique(colnames(x), sep = "_")
        
        return(NNS.M.reg(x, y, factor.2.dummy = factor.2.dummy, point.est = point.est, plot = plot,
                         residual.plot = residual.plot, order = order, n.best = n.best, type = type,
                         location = location, noise.reduction = noise.reduction,
                         dist = dist, stn = stn, return.values = return.values, plot.regions = plot.regions,
                         point.only = point.only, ncores = ncores, confidence.interval = confidence.interval))
        
      } else { # Multivariate dim.red == FALSE
        if(is.null(original.names)){
          colnames.list <- lapply(1 : ncol(x), function(i) paste0("x", i))
        } else {
          colnames.list <- original.names
        }
        
        x <- apply(data.matrix(x), 2, as.numeric)
        y <- as.numeric(y)
        
        if(!is.null(dim.red.method) & !is.null(dim(x))){
          if(!is.numeric(dim.red.method)) dim.red.method <- tolower(dim.red.method)
          x.star.matrix <- matrix(nrow = length(y))
          
          if(!is.numeric(dim.red.method) && dim.red.method!="cor" && dim.red.method!="equal"){
            if(!is.null(type)) fact <- TRUE else fact <- FALSE
            
            x.star.dep <-  sapply(1:dim(x)[2], function(i) NNS.dep(x[,i], y, print.map = FALSE, asym = TRUE)$Dependence)
            
            x.star.dep[is.na(x.star.dep)] <- 0
          }

          x.star.cor <-  sapply(1:dim(x)[2], function(i) cor(x[,i], y, method = "spearman"))

          x.star.cor[is.na(x.star.cor)] <- 0
          
          if(!is.numeric(dim.red.method) && dim.red.method == "nns.dep"){
            x.star.coef <- x.star.dep
            x.star.coef[is.na(x.star.coef)] <- 0
          }
          
          if(!is.numeric(dim.red.method) && dim.red.method == "cor"){
            x.star.coef <- x.star.cor
            x.star.coef[is.na(x.star.coef)] <- 0
          }
          
          if(!is.numeric(dim.red.method) && dim.red.method == "nns.caus"){
            if(is.null(tau)){
              tau <- "cs"
            }
            x.star.coef <- numeric()

            cause <- sapply(1:dim(x)[2], function(i) Uni.caus(y, x[,i], tau = tau, plot = FALSE))

            cause[is.na(cause)] <- 0
            
            x.star.coef <- cause
          }
          
          if(!is.numeric(dim.red.method) && dim.red.method == "all"){
            if(is.null(tau)) tau <- "cs"
            
            x.star.coef.1 <- numeric()
            
            x.star.coef.1 <- sapply(1:dim(x)[2], function(i) Uni.caus(y, x[,i], tau = tau, plot = FALSE))
            
            
            x.star.coef.3 <- x.star.cor
            x.star.coef.3[is.na(x.star.coef.3)] <- 0
            x.star.coef.2 <- x.star.dep
            x.star.coef.2[is.na(x.star.coef.2)] <- 0
            x.star.coef.4 <- rep(1, ncol(x))
            x.star.coef <- apply(cbind(x.star.coef.1, x.star.coef.2, x.star.coef.3, x.star.coef.4), 1, function(x) mode(x)) 
            x.star.coef[is.na(x.star.coef)] <- 0
          }

          if(!is.numeric(dim.red.method) && dim.red.method == "equal")  x.star.coef <- rep(1, ncol(x))
          
          if(is.numeric(dim.red.method)) x.star.coef <- as.numeric(dim.red.method)
          
          preserved.coef <- x.star.coef
          x.star.coef[abs(x.star.coef) < threshold] <- 0
          
          norm.x <- apply(original.variable, 2, function(b) (b - min(b)) / (max(b) - min(b)))
          

          
          x.star.matrix <- Rfast::eachrow(norm.x, x.star.coef, "*") #t( t(norm.x) * x.star.coef)
          x.star.matrix[is.na(x.star.matrix)] <- 0
     
          #In case all IVs have 0 correlation to DV
          if(all(x.star.matrix == 0)){
            x.star.matrix <- x
            x.star.coef[x.star.coef == 0] <- preserved.coef
          }
          
          xn <- sum( abs( x.star.coef) > 0)
      
          if(is.numeric(dim.red.method)) DENOMINATOR <- sum(dim.red.method) else DENOMINATOR <- sum( abs( x.star.coef) > 0)
        
          synthetic.x.equation.coef <- data.table::data.table(Variable = colnames.list, Coefficient = x.star.coef)
        
          synthetic.x.equation <- data.table::rbindlist( list( synthetic.x.equation.coef, list("DENOMINATOR", DENOMINATOR)))
          
          
          if(!is.null(point.est)){
            new.point.est <- numeric()
            points.norm <- rbind(point.est, x)
            
            if(dist!="FACTOR"){
              points.norm <- apply(points.norm, 2, function(b) (b - min(b)) / ifelse((max(b) - min(b))==0, 1, (max(b) - min(b))))
            }
            if(is.null(np) || np==1){
              new.point.est <- sum(points.norm[1,] * x.star.coef) / xn
              
            } else {
              point.est2 <- points.norm[1:np,]
              new.point.est <- apply(point.est2, 1, function(i) as.numeric(as.vector(i)[!is.na(i)|!is.nan(i)] %*% x.star.coef[!is.na(i)|!is.nan(i)])
                                     / xn)
            }
            
            point.est <- new.point.est
            
          }
          
          x <- Rfast::rowsums(x.star.matrix / sum( abs( x.star.coef) > 0), parallel = FALSE)
          x.star <- data.table::data.table(x)
          
        }
       } # Multivariate Not NULL type
      
    } # Univariate
    
  } # Multivariate
  

  x.label <- names(x)
  if(is.null(x.label)) x.label <- "x"
   
  dependence <- tryCatch(NNS.dep(x, y, print.map = FALSE, asym = TRUE)$Dependence, error = function(e) .1)
  dependence <- dependence^2
  dependence[is.na(dependence)] <- 0
  
  rounded_dep <- ifelse((dependence*10)%%1 < .5, floor(dependence*10), ceiling(dependence*10))
  
  dep.reduced.order <- max(1, ifelse(is.null(order), rounded_dep, order))
  
  if(!is.null(order)) dep.reduced.order <- order

  if(dependence == 1 || dep.reduced.order == "max"){
    if(is.null(order)) dep.reduced.order <- "max"
    part.map <- NNS.part(x, y, order = dep.reduced.order, obs.req = 0)
  } else {
    if(is.null(type)){
      noise.reduction2 <- ifelse(noise.reduction=="mean", "off", noise.reduction)
    } else {
      if(type == "class") noise.reduction2 <- "mode_class" else noise.reduction2 <- noise.reduction
    }
    
    if(dep.reduced.order == "max"){
      part.map <- NNS.part(x, y, order = dep.reduced.order, obs.req = 0)
    } else {
      part.map <- NNS.part(x, y, noise.reduction = noise.reduction2, order = dep.reduced.order, type = "XONLY", obs.req = 0)
      
      part.map1 <- NNS.part(c(x, part.map$regression.points$x), c(y, part.map$regression.points$y), noise.reduction = noise.reduction2, order = dep.reduced.order, type = "XONLY", obs.req = 0)
      part.map2 <- NNS.part(c(x, part.map$regression.points$x, part.map1$regression.points$x), c(y, part.map$regression.points$y, part.map1$regression.points$y), noise.reduction = noise.reduction2, order = dep.reduced.order, type = "XONLY", obs.req = 0)
      part.map3 <- NNS.part(c(x, part.map$regression.points$x, part.map1$regression.points$x, part.map2$regression.points$x), c(y, part.map$regression.points$y, part.map1$regression.points$y, part.map2$regression.points$y), noise.reduction = noise.reduction2, order =  dep.reduced.order, type =  "XONLY", obs.req = 0)
      part.map4 <- NNS.part(c(x, part.map$regression.points$x, part.map1$regression.points$x, part.map2$regression.points$x, part.map3$regression.points$x), c(y, part.map$regression.points$y, part.map1$regression.points$y, part.map2$regression.points$y, part.map3$regression.points$y), noise.reduction = noise.reduction2, order =  dep.reduced.order, type =  "XONLY", obs.req = 0)

      if(Reduce(identical, lapply(list(part.map$regression.points$x, part.map1$regression.points$x, part.map2$regression.points$x, part.map3$regression.points$x, part.map4$regression.points$x), length)) &&
         Reduce(identical, lapply(list(part.map$regression.points$y, part.map1$regression.points$y, part.map2$regression.points$y, part.map3$regression.points$y, part.map4$regression.points$y), length))){
        
        part.map$regression.points$x <- apply(cbind(c(part.map$regression.points$x, part.map1$regression.points$x, part.map2$regression.points$x, part.map3$regression.points$x, part.map4$regression.points$x)),1, function(x) gravity(x))
        part.map$regression.points$y <- apply(cbind(c(part.map$regression.points$y, part.map1$regression.points$y, part.map2$regression.points$y, part.map3$regression.points$y, part.map4$regression.points$y)),1, function(x) gravity(x))
      } else {
        part.map <- NNS.part(c(part.map$regression.points$x, part.map1$regression.points$x, part.map2$regression.points$x, part.map3$regression.points$x, part.map4$regression.points$x),
                             c(part.map$regression.points$y, part.map1$regression.points$y, part.map2$regression.points$y, part.map3$regression.points$y, part.map4$regression.points$y),
                             noise.reduction = noise.reduction2, order = dep.reduced.order, type = "XONLY", obs.req = 0)
      }
      
      if(length(part.map$regression.points$x) > length(y)){
        i <- 0
        while(length(part.map$regression.points$x) > length(y)){
          part.map <- NNS.part(c(part.map$regression.points$x, part.map1$regression.points$x, part.map2$regression.points$x, part.map3$regression.points$x, part.map4$regression.points$x),
                               c(part.map$regression.points$y, part.map1$regression.points$y, part.map2$regression.points$y, part.map3$regression.points$y, part.map4$regression.points$y),
                               noise.reduction = noise.reduction2, order =  max(c(1, (dep.reduced.order - i))), type = "XONLY", obs.req = 0)
          i <- i + 1
        }
      }
      
      if(length(part.map$regression.points$x) == 0){
        part.map <- NNS.part(x, y, type =  "XONLY", noise.reduction = noise.reduction2, order = min( nchar(part.map$dt$quadrant)), obs.req = 0)
      }
    }
  }
  
  if(length(part.map$dt$y) > length(y)){
    part.map$dt$x <- pmax(min(x), pmin(part.map$dt$x, max(x)))
    part.map$dt[, y := gravity(y), by = "x"]
    data.table::setkey(part.map$dt, x)
    part.map$dt <- unique(part.map$dt, by = "x")
  }

  Regression.Coefficients <- data.frame(matrix(ncol = 3))
  colnames(Regression.Coefficients) <- c('Coefficient', 'X Lower Range', 'X Upper Range')
  
  regression.points <- part.map$regression.points[,.(x,y)]

  regression.points$x <- pmin(max(x), pmax(regression.points$x, min(x)))
 
  data.table::setkey(regression.points,x)
  regression.points <- regression.points[, y := gravity(y), by = "x"]
  regression.points <- unique(regression.points)
  
  
  if(type!="class" || is.null(type)){
    central_rows <- c(floor(median(1:nrow(regression.points))), ceiling(median(1:nrow(regression.points))))
    central_x <- regression.points[central_rows,]$x
    ifelse(length(unique(central_rows))>1, central_y <- gravity(y[x>=central_x[1] & x<=central_x[2]]), central_y <- regression.points[central_rows[1],]$y)
    central_x <- gravity(central_x)
    med.rps <- t(c(central_x, central_y))
  } else {
    med.rps <- t(c(NA, NA))
  }
 
  regression.points <- data.table::rbindlist(list(regression.points,data.table::data.table(do.call(rbind, list(med.rps)))), use.names = FALSE)

  regression.points <- regression.points[complete.cases(regression.points),]
  regression.points <- regression.points[ , .(x,y)]
  data.table::setkey(regression.points, x, y)
  
  ### Consolidate possible duplicated points
  regression.points <- regression.points[, y := gravity(y), by = "x"]
  regression.points <- unique(regression.points)
  

  if(dependence < 1){
    min.range <- min(regression.points$x)
    max.range <- max(regression.points$x)
    
    mid.min.range <- mean(c(min(x),min(regression.points$x)))
    mid.max.range <- mean(c(max(x),max(regression.points$x)))
    
    y.min <-  na.omit(y[x <= min.range])
    l_y.min <- length(y.min)
    l_y.min_unique <- length(unique(y.min))
    
    y.mid.min <- na.omit(y[x <= mid.min.range])
    l_y.mid.min <- length(y.mid.min)
    l_y.mid.min_unique <- length(unique(y.mid.min))
    
    x.mid.min <- na.omit(x[x <= mid.min.range])
    l_x.mid.min <- length(x.mid.min)
    l_x.mid.min_unique <- length(unique(x.mid.min))
    
    y.max <- na.omit(y[x >= max.range])
    l_y.max <- length(y.max)
    l_y.max_unique <- length(unique(y.max))
    
    y.mid.max <- na.omit(y[x >= mid.max.range])
    l_y.mid.max <- length(y.mid.max)
    l_y.mid.max_unique <- length(unique(y.mid.max))
    
    x.mid.max <- na.omit(x[x >= mid.max.range])
    l_x.mid.max <- length(x.mid.max)
    l_x.mid.max_unique <- length(unique(x.mid.max))
    
    
    ### Endpoints
    if(l_x.mid.min_unique > 1 && l_y.min > 5){
      if(dependence < stn){
        if(!is.null(type)){
          if(type=="class") x0 <- mode_class(y.min) else x0 <- unique(gravity(y[x == min(x)]))
        } else {
          if(l_y.min>1 && l_y.mid.min>1){
            x0 <- sum(lm((y[which(x <= min.range)]) ~  (x[which(x <= min.range)]))$fitted.values[which.min(x[which(x <= min.range)])]*l_y.min,
                      lm((y[which(x <= mid.min.range)]) ~  (x[which(x <= mid.min.range)]))$fitted.values[which.min(x[which(x <= mid.min.range)])]*l_y.mid.min) /
              sum(l_y.min, l_y.mid.min)
          } else {
            x0 <- y.min
          }
        }
      } else {
        if(!is.null(type)){
          if(type=="class") x0 <- mode_class(y.min) else x0 <- unique(y[x == min(x)])
        } else {
          x0 <- unique(y[x == min(x)])
        }
      }
    } else {
      if(!is.null(type)){
        if(type=="class") x0 <- mode_class(y.min) else x0 <- unique(gravity(y[x == min(x)]))
      } else {
        x0 <- unique(gravity(y[x == min(x)]))
      }
    }
    
    
    if(l_x.mid.max_unique > 1 && l_y.max > 5){
      if(dependence < stn){
        if(!is.null(type)){
          if(type=="class") x.max <- mode_class(y.max) else x.max <- unique(gravity(y[x == max(x)]))
        } else {
          if(l_y.max > 1 && l_y.mid.max > 1){
            x.max <- sum(lm(y[which(x >= max.range)] ~ x[which(x >= max.range)])$fitted.values[which.max(x[which(x >= max.range)])]*l_y.max,
                         lm(y[which(x >= mid.max.range)] ~ x[which(x >= mid.max.range)])$fitted.values[which.max(x[which(x >= mid.max.range)])]*l_y.mid.max) /
              sum(l_y.max, l_y.mid.max)
          } else{
            x.max <- y.max
          }
        }
      } else {
        if(!is.null(type)){
          if(type=="class") x.max <- mode_class(y.max) else x.max <- unique(gravity(y[x == max(x)]))
        } else {
          x.max <- unique(y[x == max(x)])
        }
      }
    } else {
      if(!is.null(type)){
        if(type=="class") x.max <- mode_class(y.max) else x.max <- unique(gravity(y[x == max(x)]))
      } else{
        x.max <- unique(gravity(y[x == max(x)]))
      }
    }
    
    ### Endpoints
    max.rps <- t(c(max(x), mean(x.max)))
    min.rps <- t(c(min(x), mean(x0)))
  } else {
    ### Endpoints
    max.rps <- t(c(max(x), y[x == max(x)][1]))
    min.rps <- t(c(min(x), y[x == min(x)][1]))
  }
  

  
  regression.points <- data.table::rbindlist(list(regression.points,data.table::data.table(do.call(rbind, list(min.rps, max.rps, med.rps )))), use.names = FALSE)
  
  regression.points <- regression.points[complete.cases(regression.points),]
  regression.points <- regression.points[ , .(x,y)]
  data.table::setkey(regression.points, x, y)
  
  ### Consolidate possible duplicated points
  regression.points <- regression.points[, y := gravity(y), by = "x"]
  regression.points <- unique(regression.points)
  

  if(dim(regression.points)[1] > 1){
    rise <- regression.points[ , 'rise' := y - data.table::shift(y)]
    run <- regression.points[ , 'run' := x - data.table::shift(x)]
  } else {
    rise <- max(y) - min(y)
    rise <- regression.points[ , 'rise' := rise]
    run <- max(x) - min(x)
    if(run==0) run <- 1
    run <- regression.points[ , 'run' := run]
    regression.points <- data.table::rbindlist(list(regression.points, regression.points, regression.points), use.names = FALSE)
  }
  
  Regression.Coefficients <- regression.points[ , .(rise,run)]
  
  Regression.Coefficients <- Regression.Coefficients[complete.cases(Regression.Coefficients), ]
  
  upper.x <- regression.points[(2 : .N), x]
  
  if(length(unique(upper.x)) > 1){
    Regression.Coefficients <- Regression.Coefficients[ , `:=` ('Coefficient'=(rise / run),'X.Lower.Range' = regression.points[-.N, x], 'X.Upper.Range' = upper.x)]
  } else {
    Regression.Coefficients <- Regression.Coefficients[ , `:=` ('Coefficient'= 0,'X.Lower.Range' = unique(upper.x), 'X.Upper.Range' = unique(upper.x))]
  }
  
  Regression.Coefficients <- Regression.Coefficients[ , .(Coefficient,X.Lower.Range, X.Upper.Range)]
  
  
  Regression.Coefficients <- unique(Regression.Coefficients)
  Regression.Coefficients[Regression.Coefficients == Inf] <- 1
  Regression.Coefficients[is.na(Regression.Coefficients)] <- 0
  
  ### Fitted Values
  p <- length(unlist(regression.points[ , 1]))
 
  
  if(is.na(Regression.Coefficients[1, Coefficient])){
    Regression.Coefficients[1, Coefficient := Regression.Coefficients[2, Coefficient] ]
  }
  if(is.na(Regression.Coefficients[.N, Coefficient])){
    Regression.Coefficients[.N, Coefficient := Regression.Coefficients[.N-1, Coefficient] ]
  }
  
  coef.interval <- findInterval(x, Regression.Coefficients[ , (X.Lower.Range)], left.open = FALSE)
  reg.interval <- findInterval(x, regression.points[, x], left.open = FALSE)

  
  if(is.fcl(order) || ifelse(is.null(order), FALSE, ifelse(order >= length(y), TRUE, FALSE))){
    estimate <- y
  } else {
    estimate <- ((x - regression.points[reg.interval, x]) * Regression.Coefficients[coef.interval, Coefficient]) + regression.points[reg.interval, y]
  }
  
  if(!is.null(point.est)){
    coef.point.interval <- findInterval(point.est, Regression.Coefficients[ , (X.Lower.Range)], left.open = FALSE, rightmost.closed = TRUE)
    reg.point.interval <- findInterval(point.est, regression.points[ , x], left.open = FALSE, rightmost.closed = TRUE)
    coef.point.interval[coef.point.interval == 0] <- 1
    reg.point.interval[reg.point.interval == 0] <- 1
    point.est.y <- as.vector(((point.est - regression.points[reg.point.interval, x]) * Regression.Coefficients[coef.point.interval, Coefficient]) + regression.points[reg.point.interval, y])
    
    if(any(point.est > max(x) | point.est < min(x) ) & length(na.omit(point.est)) > 0){
      upper.slope <- mean(tail(Regression.Coefficients[, unique(Coefficient)], 2))
      point.est.y[point.est>max(x)] <- ((point.est[point.est>max(x)] - max(x)) * upper.slope + mode(y[which.max(x)]))
      
      lower.slope <- mean(head(Regression.Coefficients[, unique(Coefficient)], 2))
      point.est.y[point.est<min(x)] <- ((point.est[point.est<min(x)] - min(x)) * lower.slope + mode(y[which.min(x)]))
    }
    
    if(!is.null(type)){
      if(type=="class") point.est.y <- pmax(min(y), pmin(max(y), ifelse(point.est.y%%1 < .5, floor(point.est.y), ceiling(point.est.y))))
    }
  }
  
  colnames(estimate) <- NULL
  if(!is.null(type)){
    if(type=="class") estimate <- pmin(max(y), pmax(min(y), ifelse(estimate%%1 < .5, floor(estimate), ceiling(estimate))))
  }

  fitted <- data.table::data.table(x = x,
                                   y = original.y,
                                   y.hat = estimate,
                                   NNS.ID = part.map$dt$quadrant)
  
  colnames(fitted) <- gsub("y.hat.V1", "y.hat", colnames(fitted))
  
  fitted$y.hat[is.na(fitted$y.hat)] <- gravity(na.omit(fitted$y.hat))
  
  Values <- cbind(x, Fitted = fitted[ , y.hat], Actual = original.y, Difference = fitted[ , y.hat] - original.y,  Accuracy = abs(round(fitted[ , y.hat]) - original.y))
  
  SE <- sqrt( sum(fitted[ , ( (y.hat - y)^2) ]) / (length(y) - 1 ))
  
  gradient <- Regression.Coefficients$Coefficient[findInterval(fitted$x, Regression.Coefficients$X.Lower.Range)]

  fitted <- cbind(fitted, gradient)
  fitted$residuals <- original.y - fitted$y.hat
  
  if(dependence < stn && mean(c(length(unique(diff(x))), length(unique(x)))) > .33*length(x)){
    bias <- fitted
    
    data.table::setkey(bias, NNS.ID)
    
    bias <- bias[, gravity(residuals)*-1, by = gradient]
    
    fitted <- fitted[bias, on=.(gradient), y.hat := y.hat + V1]
    
    bias_r <- c(bias[, bias_r := lapply(.SD, data.table::frollmean, n = 2, fill = 0, align = 'right'), .SDcols = 2]$bias_r, 0)
    bias_l <- c(0, bias[, bias_l := lapply(.SD, data.table::frollmean, n = 2, fill = 0, align = 'left'), .SDcols = 2]$bias_l)
    bias_c <- bias[, bias_c := lapply(.SD, data.table::frollmean, n = 3, fill = 0, align = 'center'), .SDcols = 2]$bias_c
    
    bias <- suppressWarnings((bias_r + bias_l + bias_c)/3)
    bias[is.na(bias)] <- 0
    
    if(!is.null(type)){
      if(type=="class") suppressWarnings(regression.points[, y := ifelse((y + bias)%%1 < 0.5, floor(y + bias), ceiling(y + bias))]) else suppressWarnings(regression.points[, y := y + bias])
    } else {
      suppressWarnings(regression.points[, y := y + bias])
    }
  }
  
  regression.points$x <- pmin(regression.points$x, max(x))
  regression.points$x <- pmax(regression.points$x, min(x))
  
  regression.points$y <- pmin(regression.points$y, max(y))
  regression.points$y <- pmax(regression.points$y, min(y))
  
  if(!is.numeric(order) && !is.null(order)){
    regression.points <- part.map$dt[, .(x,y)]
    data.table::setkey(regression.points, x)
  }
  
  ### Regression Equation
  if(multivariate.call)  return(regression.points)
  
  
  rise <- regression.points[ , 'rise' := y - data.table::shift(y)]
  run <- regression.points[ , 'run' := x - data.table::shift(x)]
  
  
  Regression.Coefficients <- regression.points[ , .(rise,run)]
  
  Regression.Coefficients <- Regression.Coefficients[complete.cases(Regression.Coefficients), ]
  
  upper.x <- regression.points[(2 : .N), x]
  
  Regression.Coefficients <- Regression.Coefficients[ , `:=` ('Coefficient'=(rise / run),'X.Lower.Range' = regression.points[-.N, x], 'X.Upper.Range' = upper.x)]
  
  Regression.Coefficients <- Regression.Coefficients[ , .(Coefficient,X.Lower.Range, X.Upper.Range)]
  
  
  Regression.Coefficients <- unique(Regression.Coefficients)
  Regression.Coefficients$Coefficient[Regression.Coefficients$Coefficient==Inf] <- 1
  Regression.Coefficients$Coefficient[is.na(Regression.Coefficients$Coefficient)] <- 0
  
  
  ### Fitted Values
  p <- length(unlist(regression.points[ , 1]))
  
  if(is.na(Regression.Coefficients[1, Coefficient])){
    Regression.Coefficients[1, Coefficient := Regression.Coefficients[2, Coefficient] ]
  }
  if(is.na(Regression.Coefficients[.N, Coefficient])){
    Regression.Coefficients[.N, Coefficient := Regression.Coefficients[.N-1, Coefficient] ]
  }
  
  coef.interval <- findInterval(x, Regression.Coefficients[ , (X.Lower.Range)], left.open = FALSE)
  reg.interval <- findInterval(x, regression.points[, x], left.open = FALSE)
  
  if(!is.null(order) && is.character(order)){
    estimate <- y
  } else{
    estimate <- ((x - regression.points[reg.interval, x]) * Regression.Coefficients[coef.interval, Coefficient]) + regression.points[reg.interval, y]
  }
  
  if(!is.null(point.est)){
    coef.point.interval <- findInterval(point.est, Regression.Coefficients[ , (X.Lower.Range)], left.open = FALSE, rightmost.closed = TRUE)
    reg.point.interval <- findInterval(point.est, regression.points[ , x], left.open = FALSE, rightmost.closed = TRUE)
    coef.point.interval[coef.point.interval == 0] <- 1
    reg.point.interval[reg.point.interval == 0] <- 1
    point.est.y <- as.vector(((point.est - regression.points[reg.point.interval, x]) * Regression.Coefficients[coef.point.interval, Coefficient]) + regression.points[reg.point.interval, y])
    
    if(any(point.est > max(x) | point.est < min(x) ) & length(na.omit(point.est)) > 0){
      upper.slope <- mean(tail(Regression.Coefficients[, unique(Coefficient)], 2))
      point.est.y[point.est>max(x)] <- (point.est[point.est>max(x)] - max(x)) * upper.slope +  regression.points[.N, y]
      
      lower.slope <- mean(head(Regression.Coefficients[, unique(Coefficient)], 2))
      point.est.y[point.est<min(x)] <- (point.est[point.est<min(x)] - min(x)) * lower.slope +  regression.points[1, y]
    }
    
    if(!is.null(type)){
      if(type=="class") point.est.y <- pmin(max(y), pmax(min(y), ifelse(point.est.y%%1 < .5, floor(point.est.y), ceiling(point.est.y))))
    }
  }
  
  colnames(estimate) <- NULL
  if(!is.null(type)){
    if(type=="class") estimate <- pmin(max(y), pmax(min(y), ifelse(estimate%%1 < 0.5, floor(estimate), ceiling(estimate))))
  }
  
  fitted <- data.table::data.table(x = x,
                                   y = original.y,
                                   y.hat = estimate,
                                   NNS.ID = part.map$dt$quadrant)
  
  colnames(fitted) <- gsub(".V1", "", colnames(fitted))
  
  fitted$y.hat[is.na(fitted$y.hat)] <- mode(na.omit(fitted$y.hat))
  
  Values <- cbind(x, Fitted = fitted[ , y.hat], Actual = fitted[ , y], Difference = fitted[ , y.hat] - fitted[ , y],  Accuracy = abs(round(fitted[ , y.hat]) - fitted[ , y]))
  
  SE <- sqrt( sum(fitted[ , ( (y.hat - y)^2) ]) / (length(y) - 1 ))
  
  gradient <- Regression.Coefficients$Coefficient[findInterval(fitted$x, Regression.Coefficients$X.Lower.Range)]
  
  fitted <- cbind(fitted, gradient)
  fitted$residuals <-  original.y - fitted$y.hat
  
  if(!is.null(type)){
    if(type=="class") Prediction.Accuracy <- (length(y) - sum( abs( round(fitted$y.hat) - (y)) > 0)) / length(y) else Prediction.Accuracy <- NULL
  } else {
    Prediction.Accuracy <- NULL
  }
  
  y.mean <- mean(y)
  R2 <- (sum((fitted$y - y.mean)*(fitted$y.hat - y.mean))^2)/(sum((fitted$y - y.mean)^2)*sum((fitted$y.hat - y.mean)^2))

  
  ###Standard errors estimation
  fitted[, `:=` ( 'standard.errors' = sqrt( sum((y.hat - y) ^ 2) / ( max(1,(.N - 1))) ) ), by = gradient]
  
  
  ###Confidence and prediction intervals
  pred.int = NULL
  if(is.numeric(confidence.interval)){
    fitted[, `:=` ( 'conf.int.pos' = abs(UPM.VaR((1-confidence.interval)/2, degree = 1, residuals)) + y.hat) , by = gradient]
    fitted[, `:=` ( 'conf.int.neg' = y.hat - abs(LPM.VaR((1-confidence.interval)/2, degree = 1, residuals))) , by = gradient]
  
    if(!is.null(point.est)){
      lower.pred.int <- point.est.y - abs(LPM.VaR((1-confidence.interval)/2, degree = 1, fitted$residuals))
      upper.pred.int <- abs(UPM.VaR((1-confidence.interval)/2, degree = 1, fitted$residuals)) + point.est.y
    
      pred.int <- data.table::data.table(lower.pred.int, upper.pred.int)
      if(!is.null(type)&&type=="class") pred.int <- data.table::data.table(apply(pred.int, 2, function(x) ifelse(x%%1 <0.5, floor(x), ceiling(x))))
    }
  }
  
  ###Plotting and regression equation
  if(plot){
    if(!is.null(type) && type=="class") r2.leg <- paste("Accuracy: ", format(Prediction.Accuracy, digits = 4)) else  r2.leg <- bquote(bold(R ^ 2 == .(format(R2, digits = 4))))
    xmin <- min(c(point.est, x))
    xmax <- max(c(point.est, x))
    ymin <- min(c(point.est.y, y, fitted$y.hat, regression.points$y))
    ymax <- max(c(point.est.y, y, fitted$y.hat, regression.points$y))
    
    if(is.null(order)){
      plot.order <- max(1, part.map$order)
    } else {
      plot.order <- max(1, order)
    }
    
    if(is.numeric(confidence.interval)){
      plot(x, y, xlim = c(xmin, xmax), pch = 1, lwd = 2,
           ylim = c(min(c(fitted$conf.int.neg, ymin)), max(c(fitted$conf.int.pos,ymax))),
           col ='steelblue', main = paste(paste0("NNS Order = ", plot.order), sep = "\n"),
           xlab = if(!is.null(original.columns)){
             if(original.columns > 1){
               "Synthetic X*"
             } else { x.label }
           } else {
             x.label
           },
           ylab = y.label, mgp = c(2.5, 0.5, 0),
           cex.lab = 1.5, cex.main = 2)

      idx <- order(fitted$x)
      polygon(c(x[idx], x[rev(idx)]), c(na.omit(fitted$conf.int.pos[idx]), (na.omit(fitted$conf.int.neg[rev(idx)]))), col = "pink", border = NA)
      points(x[idx], y[idx], pch = 1, lwd = 2, col = "steelblue")
    } else {
      plot(x, y, pch = 1, lwd = 2, xlim = c(xmin, xmax), ylim = c(ymin, ymax),col = 'steelblue', main = paste(paste0("NNS Order = ", plot.order), sep = "\n"),
           xlab = if(!is.null(original.columns)){
             if(original.columns > 1){
               "Synthetic X*"
             } else { x.label }
           } else {
             x.label
           },
           ylab = y.label, mgp = c(2.5, 0.5, 0),
           cex.lab = 1.5, cex.main = 2)
    } # !confidence.intervals
    
    ### Plot Regression points and fitted values and legend
    points(na.omit(regression.points[ , .(x,y)]), col = 'red', pch = 15)
    lines(na.omit(regression.points[ , .(x,y)]), col = 'red', lwd = 2, lty = 2)
  
    
    if(!is.null(point.est)){
      points(point.est, point.est.y, col='green', pch = 18, cex = 1.5)
      legend(location, bty = "n", y.intersp = 0.75, legend = r2.leg)
    } else {
      legend(location, bty = "n", y.intersp = 0.75, legend = r2.leg)
    }
    
    if(!is.null(point.est)){
      if(any(point.est > max(x))){
        segments(point.est[point.est > max(x)], point.est.y[point.est > max(x)], regression.points[.N, x], regression.points[.N, y], col = "green", lty = 2)
      }
      
      if(any(point.est < min(x))){
        segments(point.est[point.est < min(x)], point.est.y[point.est < min(x)], regression.points[1, x], regression.points[1, y], col = "green", lty = 2)
      }
    }
  }# plot TRUE bracket
  
  options(warn = oldw)
  
  
  ### Return Values
  if(return.values){
    return(list("R2" = R2,
                "SE" = SE,
                "Prediction.Accuracy" = Prediction.Accuracy,
                "equation" = synthetic.x.equation,
                "x.star" = x.star,
                "derivative" = Regression.Coefficients[],
                "Point.est" = point.est.y,
                "pred.int" = pred.int,
                "regression.points" = regression.points[, .(x,y)],
                "Fitted.xy" = fitted))
  } else {
    invisible(list("R2" = R2,
                   "SE" = SE,
                   "Prediction.Accuracy" = Prediction.Accuracy,
                   "equation" = synthetic.x.equation,
                   "x.star" = x.star,
                   "derivative" = Regression.Coefficients[],
                   "Point.est" = point.est.y,
                   "pred.int" = pred.int,
                   "regression.points" = regression.points[ ,.(x,y)],
                   "Fitted.xy" = fitted))
  }
  
}
