#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point numeric or ("overall"); \code{x} point to be evaluated, must be provided.  Defaults to \code{(eval.point = NULL)}.  Set to \code{(eval.point = "overall")} to find an overall partial derivative estimate (1st derivative only).
#' @return Returns a \code{data.frame} of eval.point along with both 1st and 2nd derivative.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"  \doi{10.1007/s10614-017-9713-5}
#'
#' @examples
#' \dontrun{
#' x <- seq(0, 2 * pi, pi / 100) ; y <- sin(x)
#' dy.dx(x, y, eval.point = 1.75)
#' 
#' # First derivative
#' dy.dx(x, y, eval.point = 1.75)$first.derivative
#'
#' # Second derivative
#' dy.dx(x, y, eval.point = 1.75)$second.derivative
#' 
#' # Vector of derivatives
#' dy.dx(x, y, eval.point = c(1.75, 2.5))
#' }
#' @export

dy.dx <- function(x, y, eval.point = NULL){
  
  if(is.null(eval.point) || length(eval.point) == 0)
    stop("eval.point must be provided: numeric point(s) or 'overall'")

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
    # Typical spacing between neighbouring x values, independent of row order.
    x.spacing <- gravity(diff(sort(unique(x))))

    results <- vector(mode = "list", length(h_s))
    first.deriv <- vector(mode = "list", length(h_s))
    second.deriv <- vector(mode = "list", length(h_s))
    deriv.points <- vector(mode = "list", length(h_s))
    grads <- vector(mode = "numeric", length(h_s))
  
    for(index in seq_along(h_s)){
      h_step <- x.spacing * h_s[index]

      eval.point.min <- pmax(min(x), original.eval.point.min - h_step)
      eval.point.max <- pmin(max(x), h_step + original.eval.point.max)
      
      deriv.points[[index]] <- cbind(eval.point.min, eval.point, eval.point.max)
    }

    deriv.points <- do.call(rbind.data.frame, deriv.points)

    deriv.points <- deriv.points[order(deriv.points[, "eval.point"], method = "radix"), , drop = FALSE]
    rownames(deriv.points) <- NULL

      # At the boundary one run is 0 (the window is clipped to [min(x), max(x)]);
      # the first derivative then uses the one-sided slope below.

      reg.output <- NNS.reg(x, y, plot = FALSE, point.est = unlist(deriv.points), point.only = TRUE, ncores = 1, smooth = TRUE)

      combined.matrices <- cbind(deriv.points, matrix(unlist(reg.output$Point.est), ncol = 3, byrow = F))
      colnames(combined.matrices) <- c(colnames(deriv.points), "estimates.min", "estimates", "estimates.max")

      # Per-row run/rise computed from the combined columns (matches prior data.table `:=` scoping)
      eval.point.col <- combined.matrices[, "eval.point"]
      run_1c  <- combined.matrices[, "eval.point.max"] - eval.point.col
      run_2c  <- eval.point.col - combined.matrices[, "eval.point.min"]
      rise_1c <- combined.matrices[, "estimates.max"] - combined.matrices[, "estimates"]
      rise_2c <- combined.matrices[, "estimates"]     - combined.matrices[, "estimates.min"]

      first.deriv.vec  <- (rise_1c + rise_2c) / (run_1c + run_2c)
      second.deriv.vec <- (rise_1c / run_1c - rise_2c / run_2c) / ((run_1c + run_2c) / 2)

      # Grouped means by eval.point (ascending), NA-preserving like data.table mean()
      unique.eval  <- sort(unique(eval.point.col))
      first.deriv  <- data.frame(eval.point = unique.eval,
                                 first.derivative  = vapply(unique.eval,
                                   function(p) mean(first.deriv.vec[eval.point.col == p]), numeric(1)))
      second.deriv <- data.frame(eval.point = unique.eval,
                                 second.derivative = vapply(unique.eval,
                                   function(p) mean(second.deriv.vec[eval.point.col == p]), numeric(1)))

  }

  colnames(first.deriv) <- c("eval.point", "first.derivative")
  colnames(second.deriv) <- c("eval.point", "second.derivative")
  
  return(.NNS.out(merge(first.deriv, second.deriv, by = "eval.point")))
}


