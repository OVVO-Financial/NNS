#' NNS SSD Test
#'
#' Bi-directional test of second degree stochastic dominance using lower partial moments.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param plot logical; \code{TRUE} (default) plots the SSD test.
#' @return Returns one of the following SSD results: \code{"X SSD Y"}, \code{"Y SSD X"}, or \code{"NO SSD EXISTS"}.
#' @note The curves are compared on the merged sample points and, for TSD and continuous FSD, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous FSD) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}.
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.SSD(x, y)
#' }
#' @export


NNS.SSD <- function(x, y, plot = TRUE){

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

  LPM_x_sort <- LPM(1, Combined_sort,x)
  LPM_y_sort <- LPM(1, Combined_sort,y)


  if(plot){
    plot(Combined_sort, LPM_x_sort, type = "l", lwd = 3,col = "red", main = "SSD", ylab = "Area of Cumulative Distribution",
          ylim = c(min(c(LPM_y_sort, LPM_x_sort)), max(c(LPM_y_sort, LPM_x_sort))))

    lines(Combined_sort, LPM_y_sort, type = "l", lwd = 3,col = "steelblue")
    legend("topleft", c("X", "Y"), lwd = 10, col = c("red", "steelblue"))
  }

  # Verdict from the pairwise test (NNS.SSD.uni both ways); the grid curves are for the plot.
  x.dom <- .Call(`_NNS_NNS_SSD_uni_cpp`, x, y) == 1L
  y.dom <- .Call(`_NNS_NNS_SSD_uni_cpp`, y, x) == 1L
  ifelse(x.dom,
         "X SSD Y",
          ifelse (y.dom,
                  "Y SSD X",
                   "NO SSD EXISTS"))

}

