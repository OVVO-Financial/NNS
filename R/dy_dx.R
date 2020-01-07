#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivate of \code{y} wrt \code{x} for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param eval.point numeric; \code{x} point to be evaluated.  Defaults to \code{(eval.point = median(x))}.  Set to \code{(eval.point = "overall")} to find an overall partial derivative estimate.
#' @param deriv.order numeric options: (1, 2); 1 (default) for first derivative.  For second derivative estimate of \code{f(x)}, set \code{(deriv.order = 2)}.
#' @param h numeric [0, ...]; Percentage step used for finite step method.  Defaults to \code{h = .05} representing a 5 percent step from the value of the independent variable.
#' @param deriv.method method of derivative estimation, options: ("NNS", "FS"); Determines the partial derivative from the coefficient of the \link{NNS.reg} output when \code{(deriv.method = "NNS")} or generates a partial derivative using the finite step method \code{(deriv.method = "FS")} (Defualt).
#' @return Returns the value of the partial derivative estimate for the given order.
#' @note If a vector of derivatives is required, ensure \code{(deriv.method = "FS")}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
#' \url{https://link.springer.com/article/10.1007/s10614-017-9713-5}
#'
#' @examples
#' \dontrun{
#' x <- seq(0, 2 * pi, pi / 100) ; y <-sin(x)
#' dy.dx(x, y, eval.point = 1.75)
#'
#' # Vector of derivatives
#' dy.dx(x, y, eval.point = c(1.75, 2.5), deriv.method = "FS")}
#' @export

dy.dx <- function(x, y, eval.point = median(x), deriv.order = 1, h = .05, deriv.method = "FS"){

  if(length(eval.point) > 1 & deriv.method == "NNS"){
    deriv.method <- "FS"
  }

  if(class(eval.point) == "character"){
    ranges <- NNS.reg(x, y, order = order, plot = FALSE)$derivative
    ranges[ , interval := seq(1 : length(ranges$Coefficient))]

    range.weights <- numeric()
    range.weights <- data.table(x, 'interval' = findInterval(x, ranges[ , X.Lower.Range]))
    ranges <- ranges[interval %in% range.weights$interval, ]

    range.weights <- range.weights[ , .N, by = 'interval']

    range.weights <- range.weights$N / sum(range.weights$N)

    return(sum(ranges[,Coefficient]*range.weights))

  } else {

    original.eval.point.min <- eval.point
    original.eval.point.max <- eval.point

    h_step <- abs(h * diff(range(x)))

    eval.point.min <- original.eval.point.min - h_step
    eval.point.max <- h_step + original.eval.point.max

    run <- eval.point.max - eval.point.min


    if(any(run == 0)) {
      z <- which(run == 0)
      eval.point.max[z] <- (abs((max(x) - min(x)) * h)) + eval.point[z]
      eval.point.max[z] <- eval.point[z] - (abs((max(x) - min(x)) * h))
      run[z] <- eval.point.max[z] - eval.point.min[z]
    }

    if(deriv.order == 1){
      if(deriv.method == "FS"){
        estimates.min <- NNS.reg(x, y, plot = FALSE, order = order, point.est = eval.point.min)$Point.est
        estimates.max <- NNS.reg(x, y, plot = FALSE, order = order, point.est = eval.point.max)$Point.est

        rise <- estimates.max - estimates.min

        return(rise / run)
      } else {

        reg.output <- NNS.reg(x, y, plot = FALSE, return.values = TRUE, order = order)

        output <- reg.output$derivative
        if(length(output[ , Coefficient]) == 1){
          return(output[ , Coefficient])
        }

        if((output[ , X.Upper.Range][which(eval.point < output[ , X.Upper.Range]) - 1][1]) < eval.point){
          return(output[ , Coefficient][which(eval.point < output[ , X.Upper.Range])][1])
        } else {
          return(mean(c(output[ , Coefficient][which(eval.point < output[ , X.Upper.Range])][1], output[ , X.Lower.Range][which(eval.point < output[ , X.Upper.Range]) - 1][1])))
        }
      }
    } else {
      ## Second derivative form:
      # f(x+h) - 2(f(x)) +f(x-h) / h^2

      h_step <- abs(h * diff(range(x)))
      deriv.points <- cbind(h_step + eval.point, eval.point, eval.point - h_step)

      second.deriv.estimates.1 <- NNS.reg(x, y, plot = FALSE, return.values = TRUE, point.est = deriv.points[ , 1])$Point.est
      second.deriv.estimates.2 <- NNS.reg(x, y, plot = FALSE, return.values = TRUE, point.est = deriv.points[ , 2])$Point.est
      second.deriv.estimates.3 <- NNS.reg(x, y, plot = FALSE, return.values = TRUE, point.est = deriv.points[ , 3])$Point.est


      f.x_h <- second.deriv.estimates.1

      two_f.x <- 2 * second.deriv.estimates.2

      f.x__h <- second.deriv.estimates.3

      run <- ((1 + h) * eval.point) - eval.point
      return((f.x_h - two_f.x + f.x__h) / (run ^ 2))

    }

  }

}
