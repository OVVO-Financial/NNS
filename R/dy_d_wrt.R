#' Partial Derivative dy/d_[wrt]
#'
#' Returns the numerical partial derivative of \code{y} with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{NNS.reg} estimates as \code{f(x + h)} and \code{f(x - h)} values.
#'
#' @param x a numeric matrix or data frame.
#' @param y a numeric vector with compatible dimensions to \code{x}.
#' @param wrt integer; Selects the regressor to differentiate with respect to (vectorized).
#' @param eval.points numeric or options: ("mean", median", "last", "all"); Regressor points to be evaluated.  \code{(eval.points = "median")} (default) to find the average partial derivative at the median of the variable with respect to.  Set to \code{(eval.points = "last")} to find the average partial derivative at the last observation of the variable with respect to (relevant for time-series data).  Set to \code{(eval.points="mean")} to find the average partial derivative at the mean of the variable with respect to. Set to \code{(eval.points = "all")} to find the average partial derivative at every observation of the variable with respect to.
#' @param mixed logical; \code{FALSE} (default) If mixed derivative is to be evaluated, set \code{(mixed = TRUE)}.
#' @param folds integer; 5 (default) Sets the number of \code{folds} in the \link{NNS.stack} procedure for optimal \code{n.best} parameter.
#' @param plot logical; \code{FALSE} (default) Set to \code{(plot = TRUE)} to view plot.
#' Default setting is \code{(noise.reduction = "mean")}.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param messages logical; \code{TRUE} (default) Prints status messages of cross-validation on \code{n.best} parameter for \link{NNS.reg}.
#' @return Returns column-wise matrix of wrt regressors:
#' \itemize{
#' \item{\code{dy.d_(...)["First", ]} the 1st derivative
#' \item{\code{dy.d_(...)["Second", ]} the 2nd derivative
#' \item{\code{dy.d_(...)["Mixed", ]} the mixed derivative (for two independent variables only).
#' }
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
#' @examples
#' \dontrun{
#' set.seed(123) ; x_1 <- runif(100) ; x_2 <- runif(100) ; y <- x_1 ^ 2 * x_2 ^ 2
#' B <- cbind(x_1, x_2)
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' only supply 1 value in [eval.points]:
#' dy.d_(B, y, wrt = 1, eval.points = c(.5))
#'
#' dy.d_(B, y, wrt = 1, eval.points = mean(B[, 1]))
#'
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' for every oberservation of 1st regressor:
#' apd <- dy.d_(B, y, wrt = 1, eval.points = "all")
#' plot(B[,1], apd["First",][[1]])
#'
#' ## 95% Confidence Interval to test if 0 is within
#' ### Lower CI
#' LPM.VaR(.025, 0, apd["First",][[1]])
#'
#' ### Upper CI
#' UPM.VaR(.025, 0, apd["First",][[1]])
#'
#' ## To find derivatives of y wrt 1st regressor and specified 2nd regressor points for both regressors
#' dy.d_(B, y, wrt = c(1, 2), eval.points = c(.5, .5))
#'
#'
#' ## Known function analysis: [y = a ^ 2 * b ^ 2]
#' x_1 <- seq(0, 1, .1) ; x_2 <- seq(0, 1, .1)
#' B <- expand.grid(x_1, x_2) ; y <- B[ , 1] ^ 2 * B[ , 2] ^ 2
#' dy.d_(B, y, wrt = 1, eval.points = c(.5, .5))
#'
#'
#' }
#' @export


dy.d_<- function(x, y, wrt,
                 eval.points = "median",
                 folds = 5,
                 mixed = FALSE,
                 plot = FALSE,
                 ncores = NULL,
                 messages = TRUE){




  n <- dim(x)[1]
  l <- dim(x)[2]


  if(is.null(l)) stop("Please ensure (x) is a matrix or data.frame type object.")

  h_step <-  abs(diff(range(x[, wrt])))/exp(1)

  order <- NULL

  if(l != 2) mixed <- FALSE

  if(is.character(eval.points)){
    if(eval.points == "median"){
      eval.points <- median(x[ , wrt])
    } else {
      if(eval.points == "last"){
        eval.points <- as.numeric(tail(x[ , wrt], 1))
      } else {
        if(eval.points == "mean"){
          eval.points <- mean(x[ , wrt])
        } else {
          if(eval.points == "all"){
            eval.points <- x[ , wrt, drop = FALSE]
          }
        }
      }
    }
  }

  original.eval.points.min <- eval.points
  original.eval.points.max <- eval.points
  original.eval.points <- eval.points

  if(messages){
    message("Currently generating NNS.reg finite difference estimates...","\r",appendLF=TRUE)
  }

  if(any(is.null(dim(eval.points)) || dim(eval.points)[2]==1)){

    if(length(eval.points)==dim(x)[2]){
        original.eval.points.min[wrt] <- original.eval.points.min[wrt] - h_step
        original.eval.points.max[wrt] <- h_step + original.eval.points.max[wrt]
    } else {
        original.eval.points.min <- original.eval.points.min - h_step
        original.eval.points.max <- h_step + original.eval.points.max
    }

    if(!is.null(dim(eval.points)) && dim(eval.points)[2]==1){
      index <- apply(sapply(quantile(unlist(eval.points), seq(.01,1,.05)), function(z) abs(z - unlist(eval.points))), 2, which.min)
      sampsize <- length(index)

      deriv.points <- x[index,]


      deriv.points <- do.call(rbind, replicate(3*sampsize*n, deriv.points, simplify=FALSE))

      deriv.points[, wrt] <- rep(c(rbind(rep(original.eval.points.min),
                                             rep(eval.points),
                                             rep(original.eval.points.max))
                                     ), each = sampsize, length.out = dim(deriv.points)[1] )


      deriv.points <- data.table::data.table(deriv.points)

      distance_wrt <- 2*h_step # (original.eval.points.max - original.eval.points.min)[1]


      position <- rep(rep(c("l", "m", "u"), each = sampsize), length.out = dim(deriv.points)[1])
      id <- rep(1:n, each = 3*sampsize, length.out = dim(deriv.points)[1])


      if(messages){
          message(paste("Currently evaluating the ", dim(deriv.points)[1], " required points"  ),"\r",appendLF=TRUE)
      }
    }


    if(length(unlist(eval.points)) == 1){
      set.seed(317)
        index <- sample.int(n = n, size = 100, replace = FALSE)
        deriv.points <- x[index, ]
        deriv.points <- do.call(rbind, replicate(3, deriv.points, simplify = FALSE))
        deriv.points[, wrt] <- c(rep(original.eval.points.min, 100),
                                 rep(eval.points, 100),
                                 rep(original.eval.points.max, 100))

        distance_wrt <- 2 * h_step

    }

    if((!is.null(dim(original.eval.points)[2]) && dim(original.eval.points)[2] > 1)  || (length(original.eval.points) > 1 && is.null(dim(original.eval.points)))){
        deriv.points <- matrix(c(original.eval.points.min, original.eval.points, original.eval.points.max), ncol = dim(x)[2], byrow = TRUE)
        distance_wrt <- 2 * h_step
    }

    estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "cor", plot = FALSE, threshold = 0)$Point.est


    if(length(unlist(eval.points)) == 1){
        lower <- mean(estimates[1:100])
        two.f.x <- 2 * mean(estimates[101:200])
        upper <- mean(estimates[201:300])
    }

    if(!is.null(dim(eval.points)) && dim(eval.points)[2] == 1){
        estimates <- data.table::data.table(cbind(estimates = estimates,
                                                  position = position,
                                                  id = id))

        lower_msd <- estimates[position=="l", sapply(.SD, function(x) list(mean=mean(as.numeric(x)), sd=sd(as.numeric(x)))), .SDcols = "estimates", by = id]
        lower <- lower_msd$V1
        lower_sd <- lower_msd$V2

        fx_msd <- estimates[position=="m", sapply(.SD, function(x) list(mean=mean(as.numeric(x)), sd=sd(as.numeric(x)))), .SDcols = "estimates", by = id]
        two.f.x <- 2* fx_msd$V1
        two.f.x_sd <- fx_msd$V2

        upper_msd <- estimates[position=="u", sapply(.SD, function(x) list(mean=mean(as.numeric(x)), sd=sd(as.numeric(x)))), .SDcols = "estimates", by = id]
        upper <- upper_msd$V1
        upper_msd <- upper_msd$V2

    }

    if((!is.null(dim(original.eval.points)[2]) && dim(original.eval.points)[2] > 1)  || (length(original.eval.points) > 1 && is.null(dim(original.eval.points)))){
        lower <- estimates[1]
        two.f.x <- 2 * estimates[2]
        upper <- estimates[3]

    }

    rise <- upper - lower

  } else {
    n <- dim(eval.points)[1]
    original.eval.points <- eval.points
    original.eval.points.min[ , wrt] <- original.eval.points.min[ , wrt] - h_step
    original.eval.points.max[ , wrt] <- h_step + original.eval.points.max[ , wrt]

    original.eval.points <- rbind(original.eval.points.min,
                                  original.eval.points,
                                  original.eval.points.max)


    estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "cor", plot = FALSE, threshold = 0)$Point.est


    lower <- head(estimates,n)
    two.f.x <- 2 * estimates[(n+1):(2*n)]
    upper <- tail(estimates,n)

    rise <- upper - lower

    distance_wrt <- 2 * h_step
  }


  if(mixed){
    if(is.null(dim(eval.points))){
      if(length(eval.points)!=2){
        stop("Mixed Derivatives are only for 2 IV")
      }
    } else {
      if(ncol(eval.points) != 2){
        stop("Mixed Derivatives are only for 2 IV")
      }
    }

    if(!is.null(dim(eval.points))){
      h_step_1 <- abs(diff(range(x[, 1])))/exp(1)
      h_step_2 <- abs(diff(range(x[, 2])))/exp(1)
      mixed.deriv.points <- matrix(c(h_step_1 + eval.points[,1], h_step_2 + eval.points[,2],
                                     eval.points[,1] - h_step_1, h_step_2 + eval.points[,2],
                                     h_step_1 + eval.points[,1], eval.points[,2] - h_step_2,
                                     eval.points[,1] - h_step_1, eval.points[,2] - h_step_2), ncol = 2, byrow = TRUE)

      mixed.distances <- 2 * (h_step_1) * 2 * (h_step_2)

    } else {
      mixed.deriv.points <- matrix(c(h_step + eval.points,
                                     eval.points[1] - h_step, h_step + eval.points[2],
                                     h_step + eval.points[1], eval.points[2] - h_step,
                                     eval.points - h_step), ncol = 2, byrow = TRUE)

      mixed.distances <- (2 * h_step) * (2 * h_step)
    }


    mixed.estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "cor", plot = FALSE, threshold = 0)$Point.est


    if(messages){
      message("Done :-)","\r",appendLF=TRUE)
    }

    z <- matrix(mixed.estimates, ncol=4, byrow=TRUE)
    z <- z[,1] + z[,4] - z[,2] - z[,3]
    mixed <- (z / mixed.distances)

    results <- list("First" = as.numeric(unlist(rise / distance_wrt)),
                    "Second" = as.numeric(unlist((upper - two.f.x + lower) / ((distance_wrt) ^ 2))),
                    "Mixed" = mixed)

    return(results)
  } else {

    results <- list("First" = as.numeric(unlist(rise / distance_wrt)),
                    "Second" = as.numeric(unlist((upper - two.f.x + lower) / ((distance_wrt) ^ 2) )))

    return(results)

  }


}

dy.d_ <- Vectorize(dy.d_, vectorize.args = "wrt")
