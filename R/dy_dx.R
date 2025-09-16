#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point numeric or ("overall"); \code{x} point to be evaluated, must be provided.  Defaults to \code{(eval.point = NULL)}.  Set to \code{(eval.point = "overall")} to find an overall partial derivative estimate (1st derivative only).
#' @return Returns a \code{data.table} of eval.point along with both 1st and 2nd derivative.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"  \doi{10.1007/s10614-017-9713-5}
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
  
  if(anyNA(cbind(x,y))) stop("You have some missing values, please address.")
  
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

    n <- length(x)
    root_n <- floor(sqrt(n))
    h_s <- round(exp(seq(log(2), log(root_n), length.out = 5)))

    results <- vector(mode = "list", length(h_s))
    first.deriv <- vector(mode = "list", length(h_s))
    second.deriv <- vector(mode = "list", length(h_s))
    deriv.points <- vector(mode = "list", length(h_s))
    grads <- vector(mode = "numeric", length(h_s))
  
    for(h in h_s){
      index <- which(h == h_s)

      h_step <- gravity(abs(diff(x))) * h_s[index]

      eval.point.min <- pmax(min(x), original.eval.point.min - h_step)
      eval.point.max <- pmin(max(x), h_step + original.eval.point.max)
      
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
    
      reg.output <- NNS.reg(x, y, plot = FALSE, point.est = unlist(deriv.points), point.only = TRUE, ncores = 1, smooth = TRUE)
     
      combined.matrices <- cbind(deriv.points, matrix(unlist(reg.output$Point.est), ncol = 3, byrow = F))
      colnames(combined.matrices) <- c(colnames(deriv.points), "estimates.min", "estimates", "estimates.max")
   
      combined.matrices[, `:=` (
        run_1 = eval.point.max - eval.point,
        run_2 = eval.point - eval.point.min,
        rise_1 = estimates.max - estimates,
        rise_2 = estimates - estimates.min
      )]

      
      combined.matrices[, `:=` (
        first.deriv = (rise_1 + rise_2) / (run_1 + run_2),
        second.deriv = (rise_1 / run_1 - rise_2 / run_2) / ((run_1 + run_2)/2)
      )]
      
      first.deriv  <- combined.matrices[, .(first.derivative  = mean(first.deriv)), by = eval.point]
      second.deriv <- combined.matrices[, .(second.derivative = mean(second.deriv)), by = eval.point]

  }

  colnames(first.deriv) <- c("eval.point", "first.derivative")
  colnames(second.deriv) <- c("eval.point", "second.derivative")
  
  return(merge(first.deriv, second.deriv, by = "eval.point"))
}


