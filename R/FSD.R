#' NNS FSD Test
#'
#' Bi-directional test of first degree stochastic dominance using lower partial moments.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @param plot logical; \code{TRUE} (default) plots the FSD test.
#' @return Returns one of the following FSD results: \code{"X FSD Y"}, \code{"Y FSD X"}, or \code{"NO FSD EXISTS"}.
#' @note The curves are compared on the merged sample points and, for TSD and continuous FSD, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous FSD) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
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



  if(plot){
    plot(Combined_sort, LPM_x_sort, type = "l", lwd = 3,col = "red", main = "FSD", ylab = "Probability of Cumulative Distribution", ylim = c(0, 1))
    lines(Combined_sort, LPM_y_sort, type = "l", lwd = 3,col = "steelblue")
    legend("topleft", c("X", "Y"), lwd = 10, col = c("red", "steelblue"))
  }
  
  ## Verification of ***0 instances*** of CDFx > CDFy, and conversely of CDFy > CDFx
  # Verdict from the pairwise test (NNS.FSD.uni both ways), which also checks
  # the continuous curves between sample points; the grid curves are for the plot.
  type.c <- if (degree == 0) "discrete" else "continuous"
  x.dom <- .Call(`_NNS_NNS_FSD_uni_cpp`, x, y, type.c) == 1L
  y.dom <- .Call(`_NNS_NNS_FSD_uni_cpp`, y, x, type.c) == 1L
  ifelse (x.dom,
          "X FSD Y",
          ifelse (y.dom,
                  "Y FSD X",
                  "NO FSD EXISTS"))

}
