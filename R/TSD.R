#' NNS TSD Test
#'
#' Bi-directional test of third degree stochastic dominance using lower partial moments.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param plot logical; \code{TRUE} (default) plots the TSD test.
#' @return Returns one of the following TSD results: \code{"X TSD Y"}, \code{"Y TSD X"}, or \code{"NO TSD EXISTS"}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{https://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.TSD(x, y)
#' @export

NNS.TSD <- function(x, y, plot = TRUE){

    if(any(class(x)=="tbl")) x <- as.vector(unlist(x))
    if(any(class(y)=="tbl")) y <- as.vector(unlist(y))

    Combined_sort <- sort(c(x, y), decreasing = FALSE)

    LPM_x_sort <- LPM(1, Combined_sort,x)
    LPM_y_sort <- LPM(1, Combined_sort,y)

    x.tsd.y <- any(LPM_x_sort > LPM_y_sort)

    y.tsd.x <- any(LPM_y_sort > LPM_x_sort)



    plot(LPM_x_sort, type = "l", lwd = 3, col = "red", main = "TSD", ylab = "Area of Cumulative Distribution",
       ylim = c(min(c(LPM_y_sort, LPM_x_sort)), max(c(LPM_y_sort, LPM_x_sort))))

    lines(LPM_y_sort, type = "l", lwd =3,col = "blue")
    legend("topleft", c("X","Y"), lwd = 10, col=c("red","blue"))

    ifelse (!x.tsd.y && min(x) >= min(y) && mean(x) >= mean(y) && !identical(LPM_x_sort, LPM_y_sort),
            "X TSD Y",
            ifelse (!y.tsd.x && min(y) >= min(x) && mean(y) >= mean(x) && !identical(LPM_x_sort, LPM_y_sort),
                    "Y TSD X",
                    "NO TSD EXISTS"))

}

