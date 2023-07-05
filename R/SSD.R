#' NNS SSD Test
#'
#' Bi-directional test of second degree stochastic dominance using lower partial moments.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param plot logical; \code{TRUE} (default) plots the SSD test.
#' @return Returns one of the following SSD results: \code{"X SSD Y"}, \code{"Y SSD X"}, or \code{"NO SSD EXISTS"}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. DOI: \doi{10.4236/jmf.2016.61012}.
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.SSD(x, y)
#' }
#' @export


NNS.SSD <- function(x, y, plot = TRUE){

    if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
    if(any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))

    if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")

    Combined_sort <- sort(c(x, y), decreasing = FALSE)

    LPM_x_sort <- LPM(1, Combined_sort,x)
    LPM_y_sort <- LPM(1, Combined_sort,y)

    x.ssd.y <- any(LPM_x_sort > LPM_y_sort)

    y.ssd.x <- any(LPM_y_sort > LPM_x_sort)


    plot(Combined_sort, LPM_x_sort, type = "l", lwd = 3,col = "red", main = "SSD", ylab = "Area of Cumulative Distribution",
         ylim = c(min(c(LPM_y_sort, LPM_x_sort)), max(c(LPM_y_sort, LPM_x_sort))))

    lines(Combined_sort, LPM_y_sort, type = "l", lwd = 3,col = "steelblue")
    legend("topleft", c("X", "Y"), lwd = 10, col = c("red", "steelblue"))


    ifelse(!x.ssd.y && min(x) >= min(y) && mean(x) >= mean(y) && !identical(LPM_x_sort, LPM_y_sort),
           "X SSD Y",
            ifelse (!y.ssd.x && min(y) >= min(x) && mean(y) >= mean(x) && !identical(LPM_x_sort, LPM_y_sort),
                    "Y SSD X",
                    "NO SSD EXISTS"))

}

