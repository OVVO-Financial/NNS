#' NNS FSD Test
#'
#' Bi-directional test of first degree stochastic dominance using lower partial moments.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns one of the following FSD results: \code{"X FSD Y"}, \code{"Y FSD X"}, or \code{"NO FSD EXISTS"}.
#' @keywords stochastic dominance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.FSD(x,y)
#' @export



NNS.FSD <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

 ## Indicator function ***for all values of x and y*** as the CDF target

  x.fsd.y=sum((LPM(0,Combined_sort,y)-LPM(0,Combined_sort,x))>=0)

  y.fsd.x=sum((LPM(0,Combined_sort,x)-LPM(0,Combined_sort,y))>=0)

  LPM_x_sort=LPM(0,Combined_sort,x)
  LPM_y_sort=LPM(0,Combined_sort,y)

    plot(Combined_sort,LPM_x_sort, type = "l", lwd =3,col = "red", main = "FSD", ylab = "Probability of Cumulative Distribution")
    lines(Combined_sort,LPM_y_sort, type = "l", lwd =3,col = "blue")
    legend("topleft", c("X","Y"), lwd=10,
           col=c("red","blue"))

     ## Verification of ***0 instances*** of CDFx > CDFy, and conversely of CDFy > CDFx
    ifelse (x.fsd.y==length(Combined) & min(x)>=min(y),"X FSD Y",
           ifelse (y.fsd.x==length(Combined) & min(y)>=min(x),"Y FSD X","NO FSD EXISTS"))
}

