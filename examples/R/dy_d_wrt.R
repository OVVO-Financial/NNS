#' Partial Derivative dy/d_[wrt]
#'
#' Returns the numerical partial derivate of \code{y} with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{NNS.reg} estimates as \code{f(x + h)} and \code{f(x - h)} values.
#'
#' @param x a numeric matrix or data frame.
#' @param y a numeric vector with compatible dimsensions to \code{x}.
#' @param wrt integer; Selects the regressor to differentiate with respect to.
#' @param eval.points numeric or options: ("mean", median", "last"); Regressor points to be evaluated.  \code{(eval.points = "median")} (default) to find partial derivatives at the median of every variable.  Set to \code{(eval.points = "last")} to find partial derivatives at the last value of every variable.  Set to \code{(eval.points="mean")} to find partial derivatives at the mean value of every variable. Set to \code{(eval.points = "all")} to find partial derivatives at every observation.
#' @param mixed logical; \code{FALSE} (default) If mixed derivative is to be evaluated, set \code{(mixed = TRUE)}.
#' @param folds integer; 5 (default) Sets the number of \code{folds} in the \link{NNS.stack} procedure for optimal \code{n.best} parameter.
#' @param plot logical; \code{FALSE} (default) Set to \code{(plot = TRUE)} to view plot.
#' Default setting is \code{(noise.reduction = "mean")}.
#' @param messages logical; \code{TRUE} (default) Prints status messages of cross-validation on \code{n.best} parameter for \link{NNS.reg}.
#' @return Returns:
#' \itemize{
#' \item{\code{dy.d_(...)$"First Derivative"}} the 1st derivative
#' \item{\code{dy.d_(...)$"Second Derivative"}} the 2nd derivative
#' \item{\code{dy.d_(...)$"Mixed Derivative"}} the mixed derivative (for two independent variables only).
#' }
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' \dontrun{
#' set.seed(123) ; x_1 <- runif(100) ; x_2 <- runif(100) ; y <- x_1 ^ 2 * x_2 ^ 2
#' B <- cbind(x_1, x_2)
#' ## To find derivatives of y wrt 1st regressor
#' dy.d_(B, y, wrt = 1, eval.points = c(.5, .5))
#'
#' ## Known function analysis: [y = a ^ 2 * b ^ 2]
#' x_1 <- seq(0, 1, .1) ; x_2 <- seq(0, 1, .1)
#' B <- expand.grid(x_1, x_2) ; y <- B[ , 1] ^ 2 * B[ , 2] ^ 2
#' dy.d_(B, y, wrt = 1, eval.points = c(.5, .5))}
#' @export


dy.d_<- function(x, y, wrt,
                 eval.points = "median",
                 folds = 5,
                 mixed = FALSE,
                 plot = FALSE,
                 messages = TRUE){

  order = NULL

  h = max(0.01, 1 - NNS.dep.hd(cbind(x,y))$Dependence^(1/exp(1)))

  if(messages){
    message("Currently determining [n.best] clusters...","\r",appendLF=TRUE)
  }

  n.best <- NNS.stack(x, y, folds = folds,
                      status = messages, method = 1,
                      order = order)$NNS.reg.n.best

  if(is.character(eval.points)){
    if(eval.points == "median"){
      eval.points = apply(x, 2, median)
    } else {
      if(eval.points == "last"){
        eval.points = as.numeric(x[length(x[ , 1]), ])
      } else {
        if(eval.points == "mean"){
          eval.points = apply(x, 2, mean)
        } else {
          if(eval.points == "all"){
            eval.points=x
          }
        }
      }
    }
  }


  original.eval.points.min <- eval.points
  original.eval.points.max <- eval.points

  if(messages){
    message("Currently generating NNS.reg finite difference estimates...","\r",appendLF=TRUE)
  }

  if(is.null(dim(eval.points))){
    original.eval.points.min[wrt] <- (1 - h) * original.eval.points.min[wrt]
    original.eval.points.max[wrt] <- (1 + h) * original.eval.points.max[wrt]

    deriv.points <- matrix(c(original.eval.points.min, eval.points, original.eval.points.max), ncol = length(eval.points), byrow = TRUE)

    estimates <- NNS.reg(x, y, order = order, point.est = deriv.points,
                         n.best = n.best,
                         residual.plot = plot, plot = plot)$Point.est

    lower <- estimates[1]
    two.f.x <- 2 * estimates[2]
    upper <- estimates[3]

    rise <- upper - lower

    distance_wrt <-  original.eval.points.max[wrt] - original.eval.points.min[wrt]
  } else {
    n <- dim(eval.points)[1]
    original.eval.points <- eval.points
    original.eval.points.min[ , wrt] <- (1 - h) * original.eval.points.min[ , wrt]
    original.eval.points.max[ , wrt] <- (1 + h) * original.eval.points.max[ , wrt]

    original.eval.points <- rbind(original.eval.points.min,
                                  original.eval.points,
                                  original.eval.points.max)

    estimates <- NNS.reg(x, y, order = order, point.est = original.eval.points,
                         n.best = n.best,
                         residual.plot = plot, plot = plot)$Point.est


    lower <- head(estimates,n)
    two.f.x <- 2 * estimates[(n+1):(2*n)]
    upper <- tail(estimates,n)

    rise <- upper - lower

    distance_wrt <- original.eval.points.max[ , wrt] - original.eval.points.min[ , wrt]
  }


  if(mixed){
    if(is.null(dim(eval.points))){
      if(length(eval.points)!=2){
        return("Mixed Derivatives are only for 2 IV")
      }
    } else {
      if(ncol(eval.points) != 2){
        return("Mixed Derivatives are only for 2 IV")
      }
    }

    if(!is.null(dim(eval.points))){

      mixed.deriv.points <- matrix(c((1 + h) * eval.points[,1],(1 + h) * eval.points[,2],
                                     (1 - h) * eval.points[,1], (1 + h) * eval.points[,2],
                                     (1 + h) * eval.points[,1], (1 - h) * eval.points[,2],
                                     (1 - h) * eval.points[,1], (1 - h) * eval.points[,2]), ncol = 2, byrow = TRUE)
      mixed.distances <- 2 * (h * abs(eval.points[,1])) * 2 * (h * abs(eval.points[,2]))

    } else {
      mixed.deriv.points <- matrix(c((1 + h) * eval.points,
                                     (1 - h) * eval.points[1], (1 + h) * eval.points[2],
                                     (1 + h) * eval.points[1], (1 - h) * eval.points[2],
                                     (1 - h) * eval.points), ncol = 2, byrow = TRUE)

      mixed.distances <- (2 * (h * abs(eval.points[1]))) * (2 * (h * abs(eval.points[2])))
    }

    mixed.estimates <- NNS.reg(x, y, order = order, point.est = mixed.deriv.points,
                               n.best = n.best,
                               residual.plot = plot, plot = plot)$Point.est

    if(messages){
      message("Done :-)","\r",appendLF=TRUE)
    }

    z <- matrix(mixed.estimates, ncol=4, byrow=TRUE)
    z <- z[,1]+z[,4]-z[,2]-z[,3]
    mixed <- (z/mixed.distances)


    return(list("First Derivative" = rise / distance_wrt,
                "Second Derivative" = (upper - two.f.x + lower) / ((distance_wrt) ^ 2),
                "Mixed Derivative" = mixed))
  } else {
    return(list("First Derivative" = rise / distance_wrt,
                "Second Derivative" = (upper - two.f.x + lower) / ((distance_wrt) ^ 2)))
  }


}
