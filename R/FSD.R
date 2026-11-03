#' NNS FSD Test
#'
#' Compares empirical first-order stochastic dominance or degree-one area-ratio (Omega/expectile) dominance using partial moments.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type \code{"discrete"} (default) compares empirical CDFs for ordinary first-order stochastic dominance. The legacy option \code{"continuous"} compares degree-one ratios \eqn{L_1/(L_1+U_1)}, equivalently Omega/expectile order; it does not carry the ordinary FSD guarantee for every increasing utility.
#' @param plot logical; \code{TRUE} (default) plots the FSD test.
#' @return With \code{type = "discrete"}, returns \code{"X FSD Y"}, \code{"Y FSD X"}, or \code{"NO FSD EXISTS"}. With \code{type = "continuous"}, returns \code{"X OMEGA Y"}, \code{"Y OMEGA X"}, or \code{"NO OMEGA DOMINANCE EXISTS"}; the numerical area-order criterion is unchanged.
#' @note The curves are compared on the merged sample points and, for TSD and continuous area order, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous area order) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \doi{10.2139/ssrn.3002675}.
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.FSD(x, y)
#' }
#' @export



NNS.FSD <- function(x, y, type = "discrete", plot = TRUE){
  type <- tolower(type)

  if(!any(type%in%c("discrete", "continuous"))) warning("type needs to be either 'discrete' or 'continuous'")
  
  to_numeric_vector <- function(v, arg_name){
    if(any(class(v)%in%c("tbl","data.table")) || is.data.frame(v) || is.matrix(v) || any(class(v) %in% c("xts", "zoo"))){
      if(!is.null(dim(v)) && ncol(v) > 1){
        stop(sprintf("%s must be a single-column object or numeric vector.", arg_name))
      }
      v <- as.vector(unlist(v, use.names = FALSE))
    }
    
    as.numeric(v)
  }
  
  x <- to_numeric_vector(x, "x")
  y <- to_numeric_vector(y, "y")


  if(anyNA(x) || anyNA(y)) stop("You have some missing values, please address.")

  Combined_sort <- sort(c(x, y), decreasing = FALSE)

  ## Indicator function ***for all values of x and y*** as the continuous CDF target
  if(type == "discrete"){
    degree <- 0
  } else {
    degree <- 1
  }

  LPM_x_sort <- LPM.ratio(degree, Combined_sort, x)
  LPM_y_sort <- LPM.ratio(degree, Combined_sort, y)



  order.name <- if (degree == 0) "FSD" else "OMEGA"
  if(plot){
    plot(Combined_sort, LPM_x_sort, type = "l", lwd = 3,col = "red",
         main = if (degree == 0) "FSD" else "Omega / expectile dominance",
         ylab = if (degree == 0) "Empirical cumulative probability" else "Lower deviation-area ratio",
         ylim = c(0, 1))
    lines(Combined_sort, LPM_y_sort, type = "l", lwd = 3,col = "steelblue")
    legend("topleft", c("X", "Y"), lwd = 10, col = c("red", "steelblue"))
  }
  
  ## Verification of ***0 instances*** of CDFx > CDFy, and conversely of CDFy > CDFx
  # Verdict from the pairwise test (NNS.FSD.uni both ways), which also checks
  # the continuous curves between sample points; the grid curves are for the plot.
  type.c <- if (degree == 0) "discrete" else "continuous"
  x.dom <- .Call(`_NNS_NNS_FSD_uni_cpp`, x, y, type.c) == 1L
  y.dom <- .Call(`_NNS_NNS_FSD_uni_cpp`, y, x, type.c) == 1L
  if (x.dom) return(paste("X", order.name, "Y"))
  if (y.dom) return(paste("Y", order.name, "X"))
  if (degree == 0) "NO FSD EXISTS" else "NO OMEGA DOMINANCE EXISTS"

}
