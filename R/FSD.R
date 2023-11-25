#' NNS FSD Test
#'
#' Bi-directional test of first degree stochastic dominance using lower partial moments.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @param plot logical; \code{TRUE} (default) plots the FSD test.
#' @return Returns one of the following FSD results: \code{"X FSD Y"}, \code{"Y FSD X"}, or \code{"NO FSD EXISTS"}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  DOI: \doi{10.4236/jmf.2016.61012}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \url{https://www.ssrn.com/abstract=3002675}.
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


  if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
  if(any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))

  if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")

  Combined_sort <- sort(c(x, y), decreasing = FALSE)

  ## Indicator function ***for all values of x and y*** as the continuous CDF target
  if(type == "discrete"){
    degree <- 0
  } else {
    degree <- 1
  }

  LPM_x_sort <- LPM.ratio(degree, Combined_sort, x)
  LPM_y_sort <- LPM.ratio(degree, Combined_sort, y)


  x.fsd.y <- any(LPM_x_sort > LPM_y_sort)

  y.fsd.x <- any(LPM_y_sort > LPM_x_sort)


  if(plot){
    plot(Combined_sort, LPM_x_sort, type = "l", lwd = 3,col = "red", main = "FSD", ylab = "Probability of Cumulative Distribution", ylim = c(0, 1))
    lines(Combined_sort, LPM_y_sort, type = "l", lwd = 3,col = "steelblue")
    legend("topleft", c("X", "Y"), lwd = 10, col = c("red", "steelblue"))
  }
  
  ## Verification of ***0 instances*** of CDFx > CDFy, and conversely of CDFy > CDFx
  ifelse (!x.fsd.y && min(x) >= min(y) && !identical(LPM_x_sort, LPM_y_sort),
          "X FSD Y",
          ifelse (!y.fsd.x && min(y) >= min(x) && !identical(LPM_x_sort, LPM_y_sort),
                  "Y FSD X",
                  "NO FSD EXISTS"))

}
