#' NNS SSD Test
#'
#' Bi-directional test of second degree stochastic dominance using lower partial moments.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns one of the following SSD results: \code{"X SSD Y"}, \code{"Y SSD X"}, or \code{"NO SSD EXISTS"}.
#' @keywords stochastic dominance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.SSD(x,y)
#' @export


NNS.SSD <- function(x,y){


  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  x.ssd.y=sum((LPM(1,Combined_sort,y)-LPM(1,Combined_sort,x))>=0)

  y.ssd.x=sum((LPM(1,Combined_sort,x)-LPM(1,Combined_sort,y))>=0)

  LPM_x_sort=LPM(1,Combined_sort,x)
  LPM_y_sort=LPM(1,Combined_sort,y)


    plot(Combined_sort,LPM_x_sort, type = "l", lwd =3,col = "red", main = "SSD", ylab = "Area of Cumulative Distribution",
         ylim = c(min(c(LPM_y_sort,LPM_x_sort)),max(c(LPM_y_sort,LPM_x_sort))))
    lines(Combined_sort,LPM_y_sort, type = "l", lwd =3,col = "blue")
    legend("topleft", c("X","Y"), lwd=10,
                 col=c("red","blue"))


    ifelse (x.ssd.y==length(Combined) & min(x)>=min(y),"X SSD Y",
            ifelse (y.ssd.x==length(Combined) & min(y)>=min(x),"Y SSD X","NO SSD EXISTS"))
}

