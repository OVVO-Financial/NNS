#' NNS FSD Test uni-directional
#'
#' Uni-directional test of first degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @return Returns (1) if \code{"X FSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  DOI: \doi{10.4236/jmf.2016.61012}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \url{https://www.ssrn.com/abstract=3002675}.
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
    if(!(min(x) >= min(y))){
      return(0)
    }
    Combined_sort <- sort(c(x, y), decreasing = FALSE)
    if(type == "discrete"){
      degree <- 0
    } else {
      degree <- 1
    }
    L.x <- LPM(degree, Combined_sort, x)
    LPM_x_sort <- L.x / (UPM(degree, Combined_sort, x) + L.x)
    L.y <- LPM(degree, Combined_sort, y)
    LPM_y_sort <- L.y / (UPM(degree, Combined_sort, y) + L.y)
    if (identical(LPM_x_sort, LPM_y_sort))
      return (0)
    
    x.fsd.y <- any(LPM_x_sort > LPM_y_sort)
    if(!x.fsd.y){
      return(1)
    }
    return(0)
}

#' NNS SSD Test uni-directional
#'
#' Uni-directional test of second degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X SSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  DOI: \doi{10.4236/jmf.2016.61012}.
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
    if(!(min(x) >= min(y)) | mean(y) > mean(x)){
        return(0)
    }
	Combined_sort <- sort(c(x, y), decreasing = FALSE)
	LPM_x_sort <- LPM(1, Combined_sort, x)
	LPM_y_sort <- LPM(1, Combined_sort, y)
    if (identical(LPM_x_sort, LPM_y_sort))
      return (0)

	x.ssd.y <- any(LPM_x_sort > LPM_y_sort)
	if(!x.ssd.y){
		return(1)
	}
	return(0)
}


#' NNS TSD Test uni-directional
#'
#' Uni-directional test of third degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X TSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  DOI: \doi{10.4236/jmf.2016.61012}.
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
    if(!(min(x) >= min(y)) | mean(y) > mean(x)) {
        return(0)
    }
	Combined_sort <- sort(c(x, y), decreasing = FALSE)
	LPM_x_sort <- LPM(2, Combined_sort, x)
	LPM_y_sort <- LPM(2, Combined_sort, y)
    if (identical(LPM_x_sort, LPM_y_sort))
      return (0)
	x.tsd.y <- any(LPM_x_sort > LPM_y_sort)
	if(!x.tsd.y){
	  return(1)
	}
	return(0)
}
