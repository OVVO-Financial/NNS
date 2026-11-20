#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivative of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point Finite numeric points, or \code{"overall"} for the average first derivative.
#' @return For numeric points, a \code{data.frame} with the first and second derivatives in the supplied order, including repeated points. For \code{"overall"}, a single first derivative.
#'
#' @details
#' \code{eval.point = "overall"} returns the mean of the \code{gradient} column of \code{NNS.reg(x, y)$Fitted.xy}. Each observation receives the slope of the unsmoothed regression segment that contains it, so the result is an observation-weighted average of segment slopes, not an average of smoothed derivatives. Very short segments between nearly coincident regression points can have extreme slopes; \code{mean(dy.dx(x, y, eval.point = x)$first.derivative)} averages the smoothed pointwise estimates instead.
#'
#' Numeric evaluation points use central differences of smoothed \link{NNS.reg} estimates over five steps from 2 to \code{floor(sqrt(n))} typical spacings of \code{x}. The left point is bounded below by \code{min(x)} and the right point above by \code{max(x)}. At coordinates outside that range, the differences use the extrapolated fit. \code{second.derivative} is a second difference of those estimates, and second differences amplify noise: for \code{y = x^2} plus noise with an R-squared near 0.99 and n = 400, estimates of the true value 2 at three interior points had a root mean squared error of about 5, and about one in five had the wrong sign. Treat it as a descriptive diagnostic on noisy data.
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
  
  points <- .nns_derivative_point(eval.point)
  x <- .nns_derivative_x(x)
  y <- .nns_table_values(y)
  if(anyNA(cbind(x,y))) stop("You have some missing values, please address.")

  if(identical(points, "overall")){
    return(mean(NNS.reg(x, y, type = "REG", order = NULL, plot = FALSE, ncores = 1)$Fitted.xy$gradient))
  }
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

  unique.points <- unique(points)
  groups <- match(windows[, 2L], unique.points)
  result <- data.frame(eval.point = unique.points,
                       first.derivative = vapply(split(first, groups), mean, numeric(1), USE.NAMES = FALSE),
                       second.derivative = vapply(split(second, groups), mean, numeric(1), USE.NAMES = FALSE))
  result <- result[match(points, unique.points), , drop = FALSE]
  rownames(result) <- NULL
  .NNS.out(result)
}
