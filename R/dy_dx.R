#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point numeric or ("overall"); \code{x} point to be evaluated, must be provided.  Defaults to \code{(eval.point = NULL)}.  Set to \code{(eval.point = "overall")} to find an overall partial derivative estimate (1st derivative only).
#' @return Returns a \code{data.table} of eval.point along with both 1st and 2nd derivative.
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
#' # First derivative
#' dy.dx(x, y, eval.point = 1.75)[ , first.derivative]
#' 
#' # Second derivative
#' dy.dx(x, y, eval.point = 1.75)[ , second.derivative]
#' 
#' # Vector of derivatives
#' dy.dx(x, y, eval.point = c(1.75, 2.5))
#' }
#' @export

dy.dx <- function(x, y, eval.point = NULL){

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
    deriv.points <- data.table::data.table(deriv.points, key = "eval.point")
    
      n <- nrow(deriv.points)

      run_1 <- deriv.points[,3] - deriv.points[,2]
      run_2 <- deriv.points[,2] - deriv.points[,1]
     
      if(any(run_1 == 0)||any(run_2 == 0)) {
        z_1 <- which(run_1 == 0); z_2 <- which(run_2 == 0)
        eval.point.max[z_1] <- ((abs((max(x) - min(x)) ))/length(x)) * index + eval.point[z_1]; eval.point.max[z_2] <- ((abs((max(x) - min(x)) ))/length(x)) * index + eval.point[z_2]
        eval.point.max[z_1] <- eval.point[z_1] - ((abs((max(x) - min(x)) ))/length(x)) * index; eval.point.max[z_2] <- eval.point[z_2] - ((abs((max(x) - min(x)) ))/length(x)) * index
        run_1[z_1] <- eval.point.max[z_1] - eval.point[z_1]; run_2[z_2] <- eval.point[z_2] - eval.point.min[z_2]
      }
    
      reg.output <- NNS.reg(x, y, plot = FALSE, point.est = unlist(deriv.points), point.only = TRUE, ncores = 1)
      
      combined.matrices <- cbind(deriv.points, matrix(unlist(reg.output$Point.est), ncol = 3, byrow = F))
      colnames(combined.matrices) <- c(colnames(deriv.points), "estimates.min", "estimates", "estimates.max")
      
      combined.matrices[, `:=` (
        run_1 = eval.point.max - eval.point,
        run_2 = eval.point - eval.point.min,
        rise_1 = estimates.max - estimates,
        rise_2 = estimates - estimates.min
      )]
      
      
      
      naive_first_gradient <- function(x){
        idx <- max(1, findInterval(x, reg.output$derivative$X.Lower.Range))
        return(reg.output$derivative[idx, "Coefficient"])
      }
      
      naive_second_gradient <- function(x){
        idx <- max(1, findInterval(x, reg.output$derivative$X.Lower.Range)-1)
        d <- reg.output$derivative[-1, ] - reg.output$derivative[-nrow(reg.output$derivative), ]
        return(d[idx, "Coefficient"]/d[idx, "X.Lower.Range"])
      }
      
      combined.matrices$naive.first.grad <- naive_first_gradient(combined.matrices$eval.point)
      combined.matrices$naive.second.grad <- naive_second_gradient(combined.matrices$eval.point)
      
      combined.matrices[, `:=` (
        first.deriv = (rise_1 + rise_2) / (run_1 + run_2),
        second.deriv = (rise_1 / run_1 - rise_2 / run_2) / mean(c(run_1, run_2))
      )]
      
      first.deriv <- tryCatch(combined.matrices[ , mean(c(first.deriv, naive.first.grad)), by = eval.point],
                              error = function(e) combined.matrices[ , mean(first.deriv), by = eval.point])
      second.deriv <- tryCatch(combined.matrices[ , mean(c(second.deriv, naive.second.grad)), by = eval.point], 
                               error = function(e) combined.matrices[ , mean(second.deriv), by = eval.point])
  }

  colnames(first.deriv) <- c("eval.point", "first.derivative")
  colnames(second.deriv) <- c("eval.point", "second.derivative")
  
  return(merge(first.deriv, second.deriv, by = "eval.point"))
}


