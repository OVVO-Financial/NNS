#' NNS FSD Test uni-directional
#'
#' Uni-directional test of first degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type \code{"discrete"} (default) uses empirical FSD. At degree 1, \code{"continuous"} retains the degree-one area-ratio (Omega/expectile) ordering, which need not imply FSD or preference for every increasing utility. SSD and TSD are unaffected by \code{type}.
#' @return Returns (1) if X dominates Y in the selected order (ordinary FSD or Omega/expectile order), else (0).
#' @note The curves are compared on the merged sample points and, for TSD and continuous area order, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous area order) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
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
  
  if(anyNA(x) || anyNA(y)) {
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
#' @note The curves are compared on the merged sample points and, for TSD and continuous area order, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous area order) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
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
  
  if(anyNA(x) || anyNA(y)) {
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
#' @note The curves are compared on the merged sample points and, for TSD and continuous area order, also at their closed-form extrema between those points. Differences smaller than \code{64 * .Machine$double.eps * S^d} (\code{S} = largest absolute value in \code{x} and \code{y}, \code{d} = 1 for SSD and the mean condition, 2 for TSD; a ratio analogue for continuous area order) are treated as ties, so that rounding noise does not decide dominance. Discrete FSD compares exact empirical CDF values.
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
  
  if(anyNA(x) || anyNA(y)) {
    stop("You have some missing values, please address.")
  }
  
  .Call(`_NNS_NNS_TSD_uni_cpp`, x, y)
}
