#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point numeric or ("overall"); \code{x} point to be evaluated, must be provided.  Defaults to \code{(eval.point = NULL)}.  Set to \code{(eval.point = "overall")} to find an overall partial derivative estimate (1st derivative only).
#' @param messages logical; \code{TRUE} (default) Prints status messages.
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

dy.dx <- function(x, y, eval.point = NULL, messages = TRUE){

  if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
  if(any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))
  
  if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")
  
  order <- NULL
  
  if(!is.null(ncol(x)) && is.null(colnames(x))){
    x <- data.frame(x)
    x <- unlist(x)
  }
  
  if(is.character(eval.point)){
    return("First" = mean(NNS.reg(x, y, order = order, plot = FALSE, ncores = 1)$Fitted.xy$gradient))
  } else {

    original.eval.point.min <- eval.point
    original.eval.point.max <- eval.point
    
    eval.point.idx <- which(eval.point==eval.point)

    h_s <- c(1:5, seq(10, 20, 5), 30)[1:min(length(x),9)]  
    
    results <- vector(mode = "list", length(h_s))
    first.deriv <- vector(mode = "list", length(h_s))
    second.deriv <- vector(mode = "list", length(h_s))
    deriv.points <- vector(mode = "list", length(h_s))
    grads <- vector(mode = "numeric", length(h_s))
  
    for(h in h_s){
      index <- which(h == h_s)
      h_step <- gravity(abs(diff(x))) * h_s[index]
      
      eval.point.min <- original.eval.point.min - h_step
      eval.point.max <- h_step + original.eval.point.max
      
      deriv.points[[index]] <- cbind(eval.point.min, eval.point, eval.point.max)
    }

    deriv.points <- do.call(rbind.data.frame, deriv.points)
    
      n <- nrow(deriv.points)

      run <- gravity(deriv.points[,3]) - gravity(deriv.points[,1])
     
      if(any(run == 0)) {
        z <- which(run == 0)
        eval.point.max[z] <- ((abs((max(x) - min(x)) ))/length(x)) * index + eval.point[z]
        eval.point.max[z] <- eval.point[z] - ((abs((max(x) - min(x)) ))/length(x)) * index
        run[z] <- eval.point.max[z] - eval.point.min[z]
      }
      
      reg.output <- NNS.reg(x, y, plot = FALSE, point.est = unlist(deriv.points), point.only = TRUE, ncores = 1)

      estimates.min <- gravity(reg.output$Point.est[1:n])
      estimates.max <- gravity(reg.output$Point.est[(2*n+1):(3*n)])
      estimates <- gravity(reg.output$Point.est[(n+1):(2*n)])
      
      rise <- estimates.max - estimates.min

      first.deriv <-  rise / run
      
      
      ## Second derivative form:
      # [f(x+h) - 2(f(x)) + f(x-h)] / h^2
      f.x__h <- estimates.min
      
      two_f.x <- 2 * estimates
      
      f.x_h <- estimates.max
      
      second.deriv <- (f.x_h - two_f.x + f.x__h) / (h_step ^ 2)
  }

  bandwidths <- list("First" = first.deriv, "Second" = second.deriv)

  return(bandwidths)
  
}

dy.dx <- Vectorize(dy.dx, vectorize.args = "eval.point")
