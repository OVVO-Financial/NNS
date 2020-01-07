#' NNS Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x a vector, matrix or data frame of variables of numeric or factor data types.
#' @param y a numeric or factor vector with compatible dimsensions to \code{x}.
#' @param factor.2.dummy logical; \code{TRUE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.
#' @param order integer; Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{(order = ...)} integer settings with \code{(noise.reduction = "off")}.  \code{(order = "max")} will force a limit condition perfect fit.
#' @param stn numeric [0, 1]; Signal to noise parameter, sets the threshold of \code{(NNS.dep)} which reduces \code{("order")} when \code{(order = NULL)}.  Defaults to 0.96 to ensure high dependence for higher \code{("order")} and endpoint determination.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all", NULL) method for determining synthetic X* coefficients.  Selection of a method automatically engages the dimension reduction regression.  The default is \code{NULL} for full multivariate regression.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "cor")} uses standard linear correlation for weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param tau options("ts", NULL); \code{NULL}(default) To be used in conjuction with \code{(dim.red.method = "NNS.caus")} or \code{(dim.red.method = "all")}.  If the regression is using time-series data, set \code{(tau = "ts")} for more accurate causal analysis.
#' @param type \code{NULL} (default).  To perform a classification, set to \code{(type = "CLASS")}.  Like a logistic regression, it is not necessary for target variable of two classes e.g. [0, 1].
#' @param point.est a numeric or factor vector with compatible dimsensions to \code{x}.  Returns the fitted value \code{y.hat} for any value of \code{x}.
#' @param location Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.
#' @param return.values logical; \code{TRUE} (default), set to \code{FALSE} in order to only display a regression plot and call values as needed.
#' @param plot  logical; \code{TRUE} (default) To plot regression.
#' @param plot.regions logical; \code{FALSE} (default).  Generates 3d regions associated with each regression point for multivariate regressions.  Note, adds significant time to routine.
#' @param residual.plot logical; \code{TRUE} (default) To plot \code{y.hat} and \code{Y}.
#' @param std.errors logical; \code{FALSE} (default) To provide standard errors of each linear segment in the \code{"Fitted.xy"} output.
#' @param confidence.interval numeric [0, 1]; \code{NULL} (default) Plots the associated confidence interval with the estimate and reports the standard error for each individual segment.
#' @param threshold  numeric [0, 1]; \code{(threshold = 0)} (default) Sets the threshold for dimension reduction of independent variables when \code{(dim.red.method)} is not \code{NULL}.
#' @param n.best integer; \code{NULL} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at \code{sqrt(# of regressors)}.  \code{(n.best = "all")} will select and weight all generated regression points.  Analogous to \code{k} in a
#' \code{k Nearest Neighbors} algorithm.  Different values of \code{n.best} are tested using cross-validation in \link{NNS.stack}.
#' @param noise.reduction the method of determing regression points options: ("mean", "median", "mode", "off"); In low signal:noise situations,\code{(noise.reduction = "mean")}  uses means for \link{NNS.dep} restricted partitions, \code{(noise.reduction = "median")}  uses medians instead of means for \link{NNS.dep} restricted partitions, while \code{(noise.reduction = "mode")}  uses modes instead of means for \link{NNS.dep} restricted partitions.  \code{(noise.reduction = "off")}  allows for maximum possible fit with a specific \code{order}.
#' @param dist options:("L1", "L2", "DTW", "FACTOR") the method of distance calculation; Selects the distance calculation used. \code{dist = "L2"} (default) selects the Euclidean distance and \code{(dist = "L1")} seclects the Manhattan distance; \code{(dist = "DTW")} selects the dynamic time warping distance; \code{(dist = "FACTOR")} uses a frequency.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param multivariate.call Internal parameter for multivariate regressions.
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
#'  \item{\code{"Fitted.xy"}} returns a \link{data.table} of \code{x},\code{y}, \code{y.hat}, \code{gradient}, and \code{NNS.ID}.
#' }
#'
#' @note Please ensure \code{point.est} is of compatible dimensions to \code{x}, error message will ensue if not compatible.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
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
#' NNS.reg(x, y, point.est = c(.25, .5, .75), dim.red.method = "cor")
#'
#' ## IRIS dataset examples:
#' # Dimension Reduction:
#' NNS.reg(iris[,1:4], iris[,5], dim.red.method = "cor", order = 5)
#'
#' # Dimension Reduction using causal weights:
#' NNS.reg(iris[,1:4], iris[,5], dim.red.method = "NNS.caus", order = 5)
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
#' NNS.reg(x, y)$derivative}
#'
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
                    std.errors = FALSE, confidence.interval = NULL,
                    threshold = 0,
                    n.best = NULL,
                    noise.reduction = "mean",
                    dist = "L2", ncores = NULL,
                    multivariate.call = FALSE){

  oldw <- getOption("warn")
  options(warn = -1)

  if(plot.regions && !is.null(order) && order == 'max'){
    stop('Please reduce the "order" or set "plot.regions = FALSE".')
  }

  if(!is.null(confidence.interval) && std.errors == FALSE){
    std.errors <- TRUE
  }

  if(!is.null(type)){
    type <- tolower(type)
    noise.reduction <- "mode"
  }

  if(class(y) == "factor"){
    type <- "class"
    noise.reduction <- "mode"
  }

  if(!plot){
    residual.plot <- FALSE
  }

  # Variable names
  original.names <- colnames(x)
  original.columns <- ncol(x)

  y.label <- deparse(substitute(y))
  if(is.null(y.label)){
    y.label <- "y"
  }


  if(is.null(original.columns)){
      x.label <- deparse(substitute(x))
      if(is.null(x.label)){
          x.label <- "x"
      }
  }

  if(factor.2.dummy && !multivariate.call){
    if(is.list(x) & !is.data.frame(x)){
      x <- do.call(cbind, x)
    }

    if(!is.null(dim(x))){
      x <- do.call(cbind, lapply(data.frame(x), factor_2_dummy_FR))
      if(is.null(colnames(x))) {colnames(x) <- colnames(x, do.NULL = FALSE)}
      colnames(x) <- make.unique(colnames(x),sep = "_")

    } else {
      x <- factor_2_dummy_FR(x)
    }

    x <- data.matrix(x)

    if(!is.null(point.est)){
        point.est.y <- numeric()
        if(is.list(point.est) & !is.data.frame(point.est)){
            point.est <- do.call(cbind, point.est)
        }

        if(!is.null(dim(x))){
            if(!is.null(dim(point.est)) && dim(point.est)[1]>=1) {
                point.est <- do.call(cbind, lapply(data.frame(point.est), factor_2_dummy_FR))
            } else {
                point.est <- t(point.est)
                point.est <- do.call(cbind, lapply(data.frame(point.est), factor_2_dummy_FR))
            }

            if(is.null(colnames(point.est)) && !is.null(dim(point.est))){
                colnames(point.est) <- colnames(x)
            }

            point.est <- as.matrix(point.est)
            l <- dim(point.est)[2]
            colnames(point.est) <- colnames(x)
        } else { # !is.null(dim(x))...implying univariate regression

            point.est <- factor_2_dummy_FR(data.frame(point.est, row.names = FALSE))
            l <- dim(t(point.est))[2]

            if(is.null(colnames(point.est))) {
              colnames(point.est) <- colnames(x)}
        }

      ### Add 0's to data for missing regressors
        if(length(colnames(x)[!(colnames(x)%in%colnames(point.est))]) > 0 ){
            Missing <- colnames(x)[!(colnames(x)%in%colnames(point.est))]
            point.est <- data.frame(point.est)
            point.est[, Missing] <- 0
            point.est <- point.est[, colnames(x)]
        }

        point.est <- data.matrix(point.est)

    } else { # is.null(point.est)
        point.est.y <- NULL
    }

  } #if(factor.2.dummy && !multivariate.call)

  # Variable names
  original.names <- colnames(x)
  original.columns <- ncol(x)

  y.label <- deparse(substitute(y))
  if(is.null(y.label)){
    y.label <- "y"
  }


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

  if(noise.reduction == 'off'){
    stn <- 0
  }

  if(!is.null(type) && type == "class" && is.null(n.best)){
    n.best <- 1
  }

  if(!is.null(ncol(original.variable))){
    if(ncol(original.variable) == 1){
      x <- original.variable
    } else {
      if(is.null(dim.red.method)){
        if(!is.null(original.columns)){
          if(is.null(n.best)){
            n.best <- floor(sqrt(original.columns))
          }
        } else {
          if(is.null(n.best)){
            n.best <- 2
          }
        }


        mult_stn <- 1 - NNS.dep.hd(cbind(x, y))$Dependence^(1/exp(1))
        if(is.na(mult_stn)) mult_stn <- stn

        return(NNS.M.reg(x, y, factor.2.dummy = factor.2.dummy, point.est = point.est, plot = plot, residual.plot = residual.plot, order = order, n.best = n.best, type = type, location = location, noise.reduction = noise.reduction, dist = dist, stn = mult_stn, return.values = return.values, plot.regions = plot.regions, ncores = ncores))
      } else { # Multivariate dim.red == FALSE

        if(is.null(original.names)){
          colnames.list <- list()
          for(i in 1 : ncol(x)){
            colnames.list[i] <- paste0("X", i)
          }
        } else {
          colnames.list <- original.names
        }

        x <- apply(data.matrix(x), 2, as.numeric)
        y <- as.numeric(y)

        if(!is.null(dim.red.method)){
          dim.red.method <- tolower(dim.red.method)
          x.star.matrix <- matrix(nrow = length(y))

          if(dim.red.method!="cor"){
            x.star.dep <- NNS.dep(cbind(x, y),print.map = FALSE)
            x.star.dep[is.na(x.star.dep)] <- 0
            x.star.cor <- cor(cbind(x, y))
            x.star.cor[is.na(x.star.cor)] <- 0
          } else {
            x.star.cor <- cor(cbind(x, y))
            x.star.cor[is.na(x.star.cor)] <- 0
          }

          if(dim.red.method == "nns.dep"){
            x.star.coef <- x.star.dep$Dependence[- (ncol(x) + 1), (ncol(x) + 1)]
            x.star.coef[is.na(x.star.coef)] <- 0
          }

          if(dim.red.method == "cor"){
            x.star.coef <- x.star.cor[- (ncol(x) + 1), (ncol(x) + 1)]
            x.star.coef[is.na(x.star.coef)] <- 0
          }

          if(dim.red.method == "nns.caus"){
            if(is.null(tau)){
              tau <- "cs"
            }
            x.star.coef <- numeric()
            for(i in 1:ncol(x)){
              cause <- NNS.caus(x[ , i], y, tau = tau, plot = FALSE)
              cause[is.na(cause)] <- 0
              x.star.coef[i] <- cause[1] - cause[2]
            }
          }

          if(dim.red.method == "all"){
            if(is.null(tau)){
              tau <- "cs"
            }

            x.star.coef.1 <- numeric()

            for(i in 1 : ncol(x)){
              cause <- NNS.caus(x[ , i], y, tau = tau, plot = FALSE)
              cause[is.na(cause)] <- 0
              x.star.coef.1[i] <- cause[1] - cause[2]
            }

            x.star.coef.3 <- x.star.cor[- (ncol(x) + 1), (ncol(x) + 1)]
            x.star.coef.3[is.na(x.star.coef.3)] <- 0
            x.star.coef.2 <- x.star.dep$Correlation[- (ncol(x) + 1), (ncol(x) + 1)]
            x.star.coef.2[is.na(x.star.coef.2)] <- 0
            x.star.coef[is.na(x.star.coef)] <- 0
            x.star.coef <- rowMeans(cbind(x.star.coef.1, x.star.coef.2, x.star.coef.3))
          }

          signs <- sign(x.star.coef)
          x.star.coef[abs(x.star.coef) < threshold] <- 0

          norm.x <- apply(original.variable, 2, function(b) (b - min(b)) / (max(b) - min(b)))

          x.star.matrix <- t( t(norm.x) * x.star.coef)
          x.star.matrix[is.na(x.star.matrix)] <- 0

          #In case all IVs have 0 correlation to DV
          if(all(x.star.matrix == 0)){
            x.star.matrix <- x
            x.star.coef[x.star.coef == 0] <- signs
          }

          DENOMINATOR <- sum( abs( x.star.coef) > 0)

          synthetic.x.equation.coef <- data.table(Variable = colnames.list, Coefficient = x.star.coef)

          synthetic.x.equation <- rbindlist( list( synthetic.x.equation.coef, list("DENOMINATOR", DENOMINATOR)))


          if(!is.null(point.est)){
            new.point.est <- numeric()
            points.norm <- rbind(point.est, original.variable)

            if(dist!="FACTOR"){
                points.norm <- apply(points.norm, 2, function(b) (b - min(b)) / (max(b) - min(b)))
            }
            if(is.null(np)|np==1){
              new.point.est <- sum(points.norm[1,] * x.star.coef) / sum( abs( x.star.coef) > 0)

            } else {
              point.est2 <- points.norm[1:np,]
              new.point.est <- sapply(1 : np, function(i) sum(point.est2[i, ] * x.star.coef) / sum( abs( x.star.coef) > 0))
            }

            point.est <- new.point.est
          }

          x <- rowSums(x.star.matrix / sum( abs( x.star.coef) > 0))

          x.star <- data.table(x.star = x)
        }

      } # Multivariate Not NULL type

    } # Univariate

  } # Multivariate

  dependence <- NNS.dep(x, y, print.map = FALSE)$Dependence ^ (1/exp(1))
  dependence[is.na(dependence)] <- 0

  if(is.null(original.columns) | is.null(dim.red.method)){
    synthetic.x.equation <- NULL
    x.star <- NULL
  }

  if(is.null(order)){
    dep.reduced.order <- floor(floor(log(length(y))) * dependence)
  } else {
    dep.reduced.order <- order
  }

  if(dependence > stn ){
    if(is.null(type)){
      part.map <- NNS.part(x, y, noise.reduction = 'median', order = dep.reduced.order, obs.req = 0, min.obs.stop = FALSE)
      if(length(part.map$regression.points$x) == 0){
        part.map <- NNS.part(x, y, noise.reduction ='median', order = min( nchar( part.map$dt$quadrant)), obs.req = 0, min.obs.stop = FALSE)
      }
    } else {
      part.map <- NNS.part(x, y, type = "XONLY",
                           noise.reduction = 'median', order = dep.reduced.order, obs.req = 0, min.obs.stop = FALSE)
      if(length(part.map$regression.points$x) == 0){
        part.map <- NNS.part(x, y,noise.reduction = 'median',type = "XONLY", order = min(nchar(part.map$dt$quadrant)), obs.req = 0, min.obs.stop = FALSE)
      }
    } # type
  } else {
    if(is.null(type)){
        noise.reduction <- "median"
        type2 <- "XONLY"
    } else {
        if(type == "class"){
          type2 = NULL
          noise.reduction <- "mode"
        } else {
            noise.reduction <- "median"
            type2 <- "XONLY"
        }
    }

    noise.reduction2 <- ifelse(noise.reduction=="mean", "median", noise.reduction)

    part.map <- NNS.part(x, y, noise.reduction = noise.reduction2,
                         order = dep.reduced.order, type = type2, min.obs.stop = FALSE)
    if(length(part.map$regression.points$x) == 0){
      part.map <- NNS.part(x, y, type =  type2, noise.reduction = noise.reduction2, order = min( nchar(part.map$dt$quadrant)), obs.req = 0, min.obs.stop = FALSE)
    }
  } # Dependence < stn


  Regression.Coefficients <- data.frame(matrix(ncol = 3))
  colnames(Regression.Coefficients) <- c('Coefficient', 'X Lower Range', 'X Upper Range')

  regression.points <- part.map$regression.points[,.(x,y)]

  min.range <- min(regression.points$x)
  max.range <- max(regression.points$x)

  mid.min.range <- mean(c(min(x),min(regression.points$x)))
  mid.max.range <- mean(c(max(x),max(regression.points$x)))

  y.mid.min <-  na.omit(y[x <= min.range])
  l_y.mid.min <- length(y.mid.min)

  y.min <- na.omit(y[x <= mid.min.range])
  l_y.min <- length(y.min)

  if(l_y.min <= 1){
    a <- median(y.min)
    b <- a
    f <- a
  } else {
    a <- median(y.min)
    if(!is.null(type)){
      if(type == "class"){
          b <- mode_class(y.min)
      } else {
          b <- mode(y.min)
      }
    } else{
      b <- mode(y.min)
    }
    f <- mean(c(max(y.min), min(y.min)))
  }

  if(l_y.mid.min <= 1){
    a1 <- median(y.mid.min)
    b1 <- a1
  } else {
    a1 <- median(y.mid.min)
    if(!is.null(type)){
      if(type == "class"){
        b1 <- mode_class(y.mid.min)
      } else {
        b1 <- mode(y.mid.min)
      }
    } else{
      b1 <- mode(y.mid.min)
    }
  }


  y.mid.max <- na.omit(y[x >= max.range])
  l_y.mid.max <- length(y.mid.max)

  y.max <- na.omit(y[x >= mid.max.range])
  l_y.max <- length(y.max)

  if(l_y.max <= 1){
    d <- median(y.max)
    e <- d
    g <- e
  } else {
    d <- median(y.max)
    if(!is.null(type)){
      if(type == "class"){
        e <- mode_class(y.max)
      } else {
        e <- mode(y.max)
      }
    } else{
      e <- mode(y.max)
    }
    g <- mean(c(max(y.max), min(y.max)))
  }

  if(l_y.mid.max <= 1){
    d1 <- median(y.mid.max)
    e1 <- d1
  } else {
    d1 <- median(y.mid.max)
    if(!is.null(type)){
      if(type == "class"){
        e1 <- mode_class(y.mid.max)
      } else {
        e1 <- mode(y.mid.max)
      }
    } else{
      e1 <- mode(y.mid.max)
    }
  }

  Dynamic.average.min <- mean(c(a, b, f))
  Dynamic.average.max <- mean(c(d, e, g))

  Dynamic.average.mid.min <- mean(c(a1, b1, Dynamic.average.min))
  Dynamic.average.mid.max <- mean(c(d1, e1, Dynamic.average.max))

  ### Endpoints
  if(length(x[x < mid.min.range]) > 1){
    if(dependence < stn){
      x0 <- Dynamic.average.min
    } else {
      x0 <- unique(y[x == min(x)])
    }
  } else {
    x0 <- unique(y[x == min(x)])
  }

  if(length(x[x > mid.max.range]) > 1){
    if(dependence < stn){
      x.max <- Dynamic.average.max
    } else {
      x.max <- unique(y[x == max(x)])
    }
  } else {
    x.max <- unique(y[x == max(x)])
  }

  ### Mid Endpoints
  x.mid.min <- Dynamic.average.mid.min

  x.mid.max <- Dynamic.average.mid.max

  mid.max.rps <- data.table(do.call(rbind,list(c(mid.max.range, mean(x.mid.max)),
                                               c(max(x),mean(x.max)))))

  mid.min.rps <- data.table(do.call(rbind,list(c(min(x), mean(x0)),
                                               c(mid.min.range, mean(x.mid.min)))))

  regression.points <- rbindlist(list(regression.points, mid.max.rps ), use.names = FALSE)

  regression.points <- rbindlist(list(regression.points, mid.min.rps ), use.names = FALSE)

  regression.points <- regression.points[complete.cases(regression.points),]
  setkey(regression.points, x, y)

### Consolidate possible duplicated points
  if(any(duplicated(regression.points$x))){
      if(noise.reduction == "mean" | noise.reduction == "off"){
          regression.points <- regression.points[, lapply(.SD, mean), .SDcols = 2, by = .(x)]
      }

      if(noise.reduction == "median"){
          regression.points <- regression.points[, lapply(.SD, median), .SDcols = 2, by = .(x)]
      }

      if(noise.reduction == "mode"){
          regression.points <- regression.points[, lapply(.SD, mode), .SDcols = 2, by = .(x)]
      }
  }

  ### Regression Equation

  if(multivariate.call){
    return(regression.points)
  }

  rise <- regression.points[ , 'rise' := y - shift(y)]
  run <- regression.points[ , 'run' := x - shift(x)]


  Regression.Coefficients <- regression.points[ , .(rise,run)]

  Regression.Coefficients <- Regression.Coefficients[complete.cases(Regression.Coefficients), ]

  upper.x <- regression.points[(2 : .N), x]

  Regression.Coefficients <- Regression.Coefficients[ , `:=` ('Coefficient'=(rise / run),'X.Lower.Range' = regression.points[-.N, x], 'X.Upper.Range' = upper.x)]

  Regression.Coefficients <- Regression.Coefficients[ , .(Coefficient,X.Lower.Range, X.Upper.Range)]


  Regression.Coefficients <- unique(Regression.Coefficients)
  Regression.Coefficients[Regression.Coefficients == Inf] <- 1
  regression.points <- regression.points[ , .(x,y)]
  setkey(regression.points)
  regression.points <- unique(regression.points)

  ### Fitted Values
  p <- length(regression.points[ , x])

  if(is.na(Regression.Coefficients[1, Coefficient])){
    Regression.Coefficients[1, Coefficient := Regression.Coefficients[2, Coefficient] ]
  }
  if(is.na(Regression.Coefficients[.N, Coefficient])){
    Regression.Coefficients[.N, Coefficient := Regression.Coefficients[.N-1, Coefficient] ]
  }

  coef.interval <- findInterval(x, Regression.Coefficients[ , (X.Lower.Range)], left.open = FALSE)
  reg.interval <- findInterval(x, regression.points[, x], left.open = FALSE)

  estimate <- ((x - regression.points[reg.interval, x]) * Regression.Coefficients[coef.interval, Coefficient]) + regression.points[reg.interval, y]

  if(!is.null(point.est)){
    coef.point.interval <- findInterval(point.est, Regression.Coefficients[ , (X.Lower.Range)], left.open = FALSE)
    reg.point.interval <- findInterval(point.est, regression.points[ , x], left.open = FALSE)
    coef.point.interval[coef.point.interval == 0] <- 1
    reg.point.interval[reg.point.interval == 0] <- 1
    point.est.y <- as.vector(((point.est - regression.points[reg.point.interval, x]) * Regression.Coefficients[coef.point.interval, Coefficient]) + regression.points[reg.point.interval, y])
    if(!is.null(type)){
        point.est.y <- round(point.est.y)
    }
  }

  colnames(estimate) <- NULL
  if(!is.null(type)){
    estimate <- round(estimate)
  }

  fitted <- data.table(x = part.map$dt$x,
                       y = part.map$dt$y,
                       y.hat = estimate,
                       NNS.ID = part.map$dt$quadrant)

  colnames(fitted) <- gsub("y.hat.V1", "y.hat", colnames(fitted))

  Values <- cbind(x, Fitted = fitted[ , y.hat], Actual = fitted[ , y], Difference = fitted[ , y.hat] - fitted[ , y],  Accuracy = abs(round(fitted[ , y.hat]) - fitted[ , y]))

  deg.fr <- length(y) - 2

  SE <- sqrt( sum(fitted[ , ( (y.hat - y)^2) ]) / deg.fr )

  y.fitted <- fitted[ , y.hat]

  gradient <- Regression.Coefficients$Coefficient[findInterval(fitted$x, Regression.Coefficients$X.Lower.Range)]

  fitted <- cbind(fitted, gradient)

  if(!is.null(type)){
    Prediction.Accuracy <- (length(y) - sum( abs( round(y.fitted) - (y)) > 0)) / length(y)
  } else {
    Prediction.Accuracy <- NULL
  }

  R2 <- (sum((fitted[ , y.hat] - mean(y)) * (y - mean(y))) ^ 2) / (sum((y - mean(y)) ^ 2) * sum((fitted[ , y.hat] - mean(y)) ^ 2))

  fitted$residuals <- fitted$y.hat - fitted$y

  ###Standard errors estimatation
  if(std.errors){
    l <- length(regression.points[ , x])
    for(i in 1 : l){
      if(i < l){
        obs <- fitted[ x >= regression.points[ , x][i] & x < regression.points[ , x][i + 1], which = TRUE]

        se.denominator <- length(obs - 2)

        if(se.denominator > 0){
          fitted[obs, `:=`
                 ( 'standard.errors' = sqrt( sum(((y.hat - y) ^ 2)) / se.denominator ) )]
        } else {
          fitted[obs, `:=` ( 'standard.errors' = 0 )]
        }
      }

      if(i == l){
        obs <- fitted[x >= regression.points[ , x][i - 1], which = TRUE]
        se.denominator <- length(obs - 2)
        fitted[obs, `:=`
               ( 'standard.errors' = sqrt(sum( ((y.hat - y) ^ 2))/ se.denominator ) )]
      }
    }
  } # std.errors

  ###Plotting and regression equation
  if(plot){
    r2.leg <- bquote(bold(R ^ 2 == .(format(R2, digits = 4))))
    xmin <- min(c(point.est, x))
    xmax <- max(c(point.est, x))
    ymin <- min(c(point.est.y, y))
    ymax <- max(c(point.est.y, y))

    if(is.null(order)){
      plot.order <- dep.reduced.order
    } else {
      plot.order <- order
    }

    if(is.numeric(confidence.interval)){
      pval <- 1 - confidence.interval
      se.max <- max(fitted[ , y.hat + qnorm(1 - (pval / 2)) * standard.errors])
      se.min <- min(fitted[, y.hat - qnorm(1 - (pval / 2)) * standard.errors])

      plot(x, y, xlim = c(xmin, xmax),
           ylim = c(min(c(se.min, ymin)), max(c(se.max,ymax))),
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

      points(na.omit(fitted[ , .(x,y.hat + qnorm(1 - (pval / 2)) * standard.errors)]), col = 'pink', pch = 19)
      points(na.omit(fitted[ , .(x,y.hat - qnorm(1 - (pval / 2)) * standard.errors)]), col = 'pink', pch = 19)

    } else {

      plot(x, y, xlim = c(xmin, xmax), ylim = c(ymin, ymax),col = 'steelblue', main = paste(paste0("NNS Order = ", plot.order), sep = "\n"),
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
      if(sum(point.est > max(x)) > 0){
        segments(point.est[point.est > max(x)], point.est.y[point.est > max(x)], regression.points[.N, x], regression.points[.N, y], col = "green", lty = 2)
      }

      if(sum(point.est < min(x)) > 0){
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
                "regression.points" = regression.points[],
                "Fitted.xy" = fitted))
  } else {
    invisible(list("R2" = R2,
                   "SE" = SE,
                   "Prediction.Accuracy" = Prediction.Accuracy,
                   "equation" = synthetic.x.equation,
                   "x.star" = x.star,
                   "derivative" = Regression.Coefficients[],
                   "Point.est" = point.est.y,
                   "regression.points" = regression.points[],
                   "Fitted.xy" = fitted))
  }

}
