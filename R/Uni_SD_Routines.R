#' NNS FSD Test uni-directional
#'
#' Uni-directional test of first degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X FSD Y"}, else (0).
#' @keywords stochastic dominance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.FSD.uni(x,y)
#' @export

NNS.FSD.uni <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = LPM(0,Combined_sort,x)
  LPM_y_sort = LPM(0,Combined_sort,y)

  if(min(y)>min(x)) {return(0)} else {

    x.fsd.y=sum((LPM_y_sort-LPM_x_sort)>=0)

  ifelse(x.fsd.y==length(Combined) & min(x)>=min(y) & !identical(LPM_x_sort,LPM_y_sort),return(1),return(0))

  }
}

#' NNS SSD Test uni-directional
#'
#' Uni-directional test of second degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X SSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.SSD.uni(x,y)
#' @export


NNS.SSD.uni <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = LPM(1,Combined_sort,x)
  LPM_y_sort = LPM(1,Combined_sort,y)

  if(min(y)>min(x) | mean(y)>mean(x)) {return(0)} else {

    x.ssd.y=sum((LPM_y_sort-LPM_x_sort)>=0)

    ifelse(x.ssd.y==length(Combined) & min(x)>=min(y) & !identical(LPM_x_sort,LPM_y_sort),return(1),return(0))

  }
}


#' NNS TSD Test uni-directional
#'
#' Uni-directional test of third degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X TSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.TSD.uni(x,y)
#' @export


NNS.TSD.uni <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = LPM(2,Combined_sort,x)
  LPM_y_sort = LPM(2,Combined_sort,y)

  if(min(y)>min(x) | mean(y)>mean(x)) {return(0)} else {

    x.tsd.y=sum((LPM_y_sort-LPM_x_sort)>=0)

    ifelse(x.tsd.y==length(Combined) & min(x)>=min(y) & !identical(LPM_x_sort,LPM_y_sort),return(1),return(0))

  }
}
