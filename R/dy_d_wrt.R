#' Partial Derivative dy/d_[wrt]
#'
#' Returns the numerical partial derivate of \code{y} with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{NNS.reg} estimates as \code{f(x + h)} and \code{f(x - h)} values.
#'
#' @param x a numeric matrix or data frame.
#' @param y a numeric vector with compatible dimsensions to \code{x}.
#' @param wrt integer; Selects the regressor to differentiate with respect to.
#' @param order integer; \link{NNS.reg} \code{"order"}, defaults to NULL.
#' @param stn numeric [0, 1]; Signal to noise parameter, sets the threshold of \link{NNS.dep} which reduces \code{"order"} when \code{(order = NULL)}.  Defaults to 0.99 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param eval.points numeric or options: ("mean", median", "last"); Regressor points to be evaluated.  \code{(eval.points = "median")} (default) to find partial derivatives at the median of every variable.  Set to \code{(eval.points = "last")} to find partial derivatives at the last value of every variable.  Set to \code{(eval.points="mean")} to find partial derivatives at the mean value of every variable. Set to \code{(eval.points = "all")} to find partial derivatives at every observation.
#' @param h numeric [0,...]; Percentage step used for finite step method.  Defaults to \code{h = .05} representing a 5 percent step from the value of the regressor.
#' @param n.best integer; Sets the number of closest regression points to use in estimating finite difference points in \link{NNS.reg}.  \code{NULL} (default) Uses \code{ceiling(sqrt(ncol(x)))}.
#' @param mixed logical; \code{FALSE} (default) If mixed derivative is to be evaluated, set \code{(mixed = TRUE)}.  Only for single valued \code{eval.points}.
#' @param plot logical; \code{FALSE} (default) Set to \code{(plot = TRUE)} to view plot.
#' @param noise.reduction the method of determing regression points options: ("mean", "median", "mode", "off"); In low signal to noise situations, \code{(noise.reduction = "median")} uses medians instead of means for partitions, while \code{(noise.reduction = "mode")} uses modes instead of means for partitions.  \code{(noise.reduction = "off")}  allows for maximum possible fit in \link{NNS.reg}.
#' Default setting is \code{(noise.reduction = "mean")}.
#' @return Returns:
#' \itemize{
#' \item{\code{dy.d_(...)$"First Derivative"}} the 1st derivative
#' \item{\code{dy.d_(...)$"Second Derivative"}} the 2nd derivative
#' \item{\code{dy.d_(...)$"Mixed Derivative"}} the mixed derivative (for two independent variables only).
#' }
#' Retuns a vector of partial derivatives when \code{(eval.points = "all")}.
#' @note For known function testing and analysis, regressors should be transformed via \link{expand.grid} to fill the dimensions with \code{(order = "max")}.  Example provided below.
#' @keywords multivaiate partial derivative
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123) ; x_1 <- runif(100) ; x_2 <- runif(100) ; y <- x_1 ^ 2 * x_2 ^ 2
#' B = cbind(x_1, x_2)
#' ## To find derivatives of y wrt 1st regressor
#' dy.d_(B, y, wrt = 1, eval.points = c(.5, .5))
#'
#' ## Known function analysis: [y = a ^ 2 * b ^ 2]
#' x_1 <- seq(0, 1, .1) ; x_2 <- seq(0, 1, .1)
#' B = expand.grid(x_1, x_2) ; y <- B[ , 1] ^ 2 * B[ , 2] ^ 2
#' dy.d_(B, y, wrt = 1, eval.points = c(.5, .5), order = "max")
#' @export


dy.d_<- function(x, y, wrt, eval.points = "median", order = NULL, stn = 0.99, h = .05, n.best = NULL, mixed = FALSE, plot = FALSE, noise.reduction = 'mean'){
  if(is.null(n.best)){
    n.best = ceiling(sqrt(ncol(x)))
  } else {
      n.best = n.best
  }

  if(is.character(eval.points)){
    if(eval.points == "median"){
      eval.points = numeric()
      eval.points = apply(x, 2, median)
    } else {
    if(eval.points == "last"){
      eval.points = numeric()
      eval.points = as.numeric(x[length(x[ , 1]), ])
    } else {
    if(eval.points == "mean"){
      eval.points = numeric()
      eval.points = apply(x, 2, mean)
    } else {
    if(eval.points == "all"){
      eval.points=x
    }
    }
    }
    }
  }


  original.eval.points.min = eval.points
  original.eval.points.max = eval.points

  if(is.null(dim(eval.points))){
  original.eval.points.min[wrt] = (1 - h) * original.eval.points.min[wrt]
  original.eval.points.max[wrt] = (1 + h) * original.eval.points.max[wrt]

  deriv.points = matrix(c(original.eval.points.min, eval.points, original.eval.points.max), ncol = length(eval.points), byrow = TRUE)

  estimates = NNS.reg(x, y, order = order, point.est = deriv.points, n.best = n.best, stn = stn, plot = plot, noise.reduction = noise.reduction)$Point.est

  lower = estimates[1]
  two.f.x = 2 * estimates[2]
  upper = estimates[3]

  rise = upper - lower

  distance.1 = sqrt(sum(sweep(t(c(original.eval.points.max)), 2, t(c(eval.points))) ^ 2))
  distance.2 = sqrt(sum(sweep(t(c(original.eval.points.min)), 2, t(c(eval.points))) ^ 2))
  run = distance.1 + distance.2
  } else {

    original.eval.points = eval.points
    original.eval.points.min[ , wrt] = (1 - h) * original.eval.points.min[ , wrt]
    original.eval.points.max[ , wrt] = (1 + h) * original.eval.points.max[ , wrt]

    estimates = NNS.reg(x, y, order = order, point.est = original.eval.points, n.best = n.best, stn = stn, plot = plot, noise.reduction = noise.reduction)$Point.est
    estimates.min = NNS.reg(x, y, order = order, point.est = original.eval.points.min, n.best = n.best, stn = stn, plot = plot, noise.reduction = noise.reduction)$Point.est
    estimates.max = NNS.reg(x, y, order = order, point.est = original.eval.points.max, n.best = n.best, stn = stn, plot = plot, noise.reduction = noise.reduction)$Point.est

    lower = estimates.min
    two.f.x = 2 * estimates
    upper = estimates.max

    rise = upper - lower

    distance.1 = rowSums(abs(original.eval.points.max - eval.points))
    distance.2 = rowSums(abs(original.eval.points.min - eval.points))
    run = distance.1 + distance.2
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

    mixed.deriv.points = matrix(c((1 + h) * eval.points,
                                (1 - h) * eval.points[1], (1 + h) * eval.points[2],
                                (1 + h) * eval.points[1], (1 - h) * eval.points[2],
                                (1 - h) * eval.points), ncol = 2, byrow = TRUE)


  mixed.estimates = NNS.reg(x, y, order = order, point.est = mixed.deriv.points, n.best = n.best, stn = stn, plot=plot, noise.reduction = noise.reduction)$Point.est

  mixed.first = mixed.estimates[1]

  mixed.second = mixed.estimates[2]

  mixed.third = mixed.estimates[3]

  mixed.fourth = mixed.estimates[4]



  return(list("First Derivative" = rise / run,
              "Second Derivative" = (upper - two.f.x + lower) / ((.5 * run) ^ 2),
              "Mixed Derivative" = (mixed.first - mixed.second - mixed.third + mixed.fourth) / (4 * ((.5 * run) ^ 2))))
  } else {
    return(list("First Derivative" = rise / run,
                "Second Derivative" = (upper - two.f.x + lower) / ((.5 * run) ^ 2)))
  }


}
