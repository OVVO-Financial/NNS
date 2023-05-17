#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point numeric or ("overall"); \code{x} point to be evaluated.  Defaults to \code{(eval.point = median(x))}.  Set to \code{(eval.point = "overall")} to find an overall partial derivative estimate (1st derivative only).
#' @return Returns a list of both 1st and 2nd derivative:
#' \itemize{
#' \item{\code{dy.dx(...)$First}} the 1st derivative.
#' \item{\code{dy.dx(...)$Second}} the 2nd derivative.
#' }
#'
#' @note If a vector of derivatives is required, ensure \code{(deriv.method = "FD")}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
#' \url{https://link.springer.com/article/10.1007/s10614-017-9713-5}
#'
#' @examples
#' \dontrun{
#' x <- seq(0, 2 * pi, pi / 100) ; y <- sin(x)
#' dy.dx(x, y, eval.point = 1.75)
#'
#' # Vector of derivatives
#' dy.dx(x, y, eval.point = c(1.75, 2.5))
#' }
#' @export

dy.dx <- function(x, y, eval.point = median(x)){
  
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
  if(any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))
  
  if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")
  
  order <- NULL
  
  dep <- NNS.dep(x, y, asym = TRUE)$Dependence
  
  h <- mean(abs(diff(LPM.VaR(seq(0, 1, min(.05, max(.01, 1-dep))), 1, x))))
  
  if(dep < .85) h <- 2*h
  
  if(!is.null(ncol(x)) && is.null(colnames(x))){
    x <- data.frame(x)
    x <- unlist(x)
  }
  
  if(is.character(eval.point)){
    return("First" = mean(NNS.reg(x, y, order = order, plot = FALSE, ncores = 1)$Fitted.xy$gradient))
  } else {
    original.eval.point.min <- eval.point
    original.eval.point.max <- eval.point
    
    h_step <- LPM.ratio(1, unlist(eval.point), x)
    h_step <- LPM.VaR(h_step + h, 1, x) - LPM.VaR(h_step - h, 1, x)
    
    eval.point.min <- original.eval.point.min - h_step
    eval.point.max <- h_step + original.eval.point.max
    
    deriv.points <- cbind(eval.point.min, eval.point, eval.point.max)
    
    n <- dim(deriv.points)[1]
    
    run <- eval.point.max - eval.point.min
    
    if(any(run == 0)) {
      z <- which(run == 0)
      eval.point.max[z] <- (abs((max(x) - min(x)) * h)) + eval.point[z]
      eval.point.max[z] <- eval.point[z] - (abs((max(x) - min(x)) * h))
      run[z] <- eval.point.max[z] - eval.point.min[z]
    }
    
    reg.output <- NNS.reg(x, y, plot = FALSE, point.est = as.vector(deriv.points), type = NULL, point.only = TRUE, ncores = 1)
    
    estimates.min <- reg.output$Point.est[1:n]
    estimates.max <- reg.output$Point.est[(2*n+1):(3*n)]
    estimates <- reg.output$Point.est[(n+1):(2*n)]
    
    
    rise <- estimates.max - estimates.min
    
    first.deriv <-  rise / run
    
    
    ## Second derivative form:
    # [f(x+h) - 2(f(x)) + f(x-h)] / h^2
    f.x__h <- estimates.min
    
    two_f.x <- 2 * estimates
    
    f.x_h <- estimates.max
    
    second.deriv <- (f.x_h - two_f.x + f.x__h) / (h_step ^ 2)
    
    return(list("First" = first.deriv,
                "Second" = second.deriv))
    
  }
}
