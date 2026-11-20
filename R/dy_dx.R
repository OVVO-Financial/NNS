#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point Numeric points in the observed range of \code{x}, or \code{"overall"} for the average first derivative.
#' @return For numeric points, a \code{data.frame} with the first and second derivatives in the supplied order, including repeated points. For \code{"overall"}, a single first derivative.
#'
#' @details
#' \code{eval.point = "overall"} returns the mean of the \code{gradient} column of \code{NNS.reg(x, y)$Fitted.xy}. Each observation receives the slope of the unsmoothed regression segment that contains it, so the result is an observation-weighted average of segment slopes, not an average of smoothed derivatives. Very short segments between nearly coincident regression points can have extreme slopes; \code{mean(dy.dx(x, y, eval.point = x)$first.derivative)} averages the smoothed pointwise estimates instead.
#'
#' Numeric evaluation points use central differences of smoothed \link{NNS.reg} estimates over five steps from 2 to \code{floor(sqrt(n))} typical spacings of \code{x}, clipped to the observed range. \code{second.derivative} is a second difference of those estimates, and second differences amplify noise: for \code{y = x^2} plus noise with an R-squared near 0.99 and n = 400, estimates of the true value 2 at three interior points had a root mean squared error of about 5, and about one in five had the wrong sign. Treat it as a descriptive diagnostic on noisy data.
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
  
  if(!is.null(ncol(x)) && is.null(colnames(x))){
    x <- data.frame(x)
    x <- unlist(x)
  }
  
  if(identical(unname(eval.point), "overall")){
    return(mean(NNS.reg(x, y, type = "REG", order = NULL, plot = FALSE, ncores = 1)$Fitted.xy$gradient))
  }
  if(!is.numeric(eval.point))
    stop("eval.point must be numeric point(s) or 'overall'")
  if(any(!is.finite(eval.point)) || any(eval.point < min(x) | eval.point > max(x)))
    stop("Derivative evaluation points must lie inside the observed coordinate support.")

  points <- unique(as.numeric(eval.point))
  root_n <- floor(sqrt(length(x)))
  h_s <- round(exp(seq(log(2), log(root_n), length.out = 5)))
  x.spacing <- gravity(diff(sort(unique(x))))
  windows <- lapply(h_s * x.spacing, function(h) {
    cbind(pmax(min(x), points - h), points, pmin(max(x), points + h))
  })
  windows <- do.call(rbind, windows)
  fit <- NNS.reg(x, y, type = "REG", plot = FALSE, point.est = as.vector(windows),
                 point.only = TRUE, ncores = 1, smooth = TRUE)
  estimates <- matrix(as.numeric(fit$Point.est), ncol = 3L)

  left <- windows[, 2L] - windows[, 1L]
  right <- windows[, 3L] - windows[, 2L]
  rise.left <- estimates[, 2L] - estimates[, 1L]
  rise.right <- estimates[, 3L] - estimates[, 2L]
  first <- (rise.right + rise.left) / (right + left)
  second <- (rise.right / right - rise.left / left) / ((right + left) / 2)

  idx <- match(eval.point, points)
  .NNS.out(data.frame(eval.point = points[idx],
                      first.derivative = rowMeans(matrix(first, nrow = length(points)))[idx],
                      second.derivative = rowMeans(matrix(second, nrow = length(points)))[idx]))
}
