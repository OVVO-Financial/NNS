#' NNS FSD Test uni-directional
#'
#' Uni-directional test of first degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @return Returns (1) if \code{"X FSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance."  \doi{10.2139/ssrn.3002675}
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.FSD.uni(x, y)
#' }
#' @export

NNS.FSD.uni <- function(x, y, type = "discrete"){
    if(any(class(x)%in%c("tbl","data.table"))) { 
      x <- as.vector(unlist(x))
    }
    if(any(class(y)%in%c("tbl","data.table"))) {
      y <- as.vector(unlist(y))
    }
    if(sum(is.na(cbind(x,y))) > 0) {
      stop("You have some missing values, please address.")
    }
    type <- tolower(type)
    if(!any(type %in% c("discrete", "continuous"))) {
      warning("type needs to be either discrete or continuous")
    }
    .Call(`_NNS_NNS_FSD_uni_cpp`, x, y, as.character(type))
}

#' NNS SSD Test uni-directional
#'
#' Uni-directional test of second degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X SSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}.
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.SSD.uni(x, y)
#' }
#' @export

NNS.SSD.uni <- function(x, y){
    if(any(class(x) %in% c("tbl","data.table"))){
	  x <- as.vector(unlist(x))
	}
    if(any(class(y) %in% c("tbl","data.table"))){
	  y <- as.vector(unlist(y))
	}
    if(sum(is.na(cbind(x,y))) > 0){
	  stop("You have some missing values, please address.")
	}
  .Call(`_NNS_NNS_SSD_uni_cpp`, x, y)
}


#' NNS TSD Test uni-directional
#'
#' Uni-directional test of third degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X TSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}.
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.TSD.uni(x, y)
#' }
#' @export

NNS.TSD.uni <- function(x, y){
    if(any(class(x)%in%c("tbl","data.table"))){
	  x <- as.vector(unlist(x))
	}
    if(any(class(y)%in%c("tbl","data.table"))){
	  y <- as.vector(unlist(y))
	}
    if(sum(is.na(cbind(x,y))) > 0){
	  stop("You have some missing values, please address.")
	}
  .Call(`_NNS_NNS_TSD_uni_cpp`, x, y)
}
