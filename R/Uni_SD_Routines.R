#' NNS FSD Test uni-directional
#'
#' Uni-directional test of first degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @return Returns (1) if \code{"X FSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{https://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \url{https://www.ssrn.com/abstract=3002675}.
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.FSD.uni(x, y)
#' @export

NNS.FSD.uni <- function(x, y, type = "discrete"){

    if(any(class(x)=="tbl")) x <- as.vector(unlist(x))
    if(any(class(y)=="tbl")) y <- as.vector(unlist(y))

    type <- tolower(type)

    if(!any(type%in%c("discrete", "continuous"))){
        warning("type needs to be either discrete or continuous")
    }

    if(min(y) > min(x)){
        return(0)
    } else {
        x_sort <- sort(x, decreasing = FALSE)
        y_sort <- sort(y, decreasing = FALSE)

        Combined <- c(x_sort, y_sort)
        Combined_sort <- sort(Combined, decreasing = FALSE)

        if(type == "discrete"){
            degree <- 0
        } else {
            degree <- 1
        }

        L.x <- LPM(degree, Combined_sort, x)
        LPM_x_sort <- L.x / (UPM(degree, Combined_sort, x) + L.x)
        L.y <- LPM(degree, Combined_sort, y)
        LPM_y_sort <- L.y / (UPM(degree, Combined_sort, y) + L.y)

        x.fsd.y <- any(LPM_x_sort > LPM_y_sort)

        ifelse(!x.fsd.y & min(x) >= min(y) & !identical(LPM_x_sort, LPM_y_sort),
               return(1),
               return(0))

    }
}

#' NNS SSD Test uni-directional
#'
#' Uni-directional test of second degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X SSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{https://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.SSD.uni(x, y)
#' @export


NNS.SSD.uni <- function(x, y){

    if(any(class(x)=="tbl")) x <- as.vector(unlist(x))
    if(any(class(y)=="tbl")) y <- as.vector(unlist(y))

    if(min(y) > min(x) | mean(y) > mean(x)) {
        return(0)
    } else {
        x_sort <- sort(x, decreasing = FALSE)
        y_sort <- sort(y, decreasing = FALSE)

        Combined <- c(x_sort, y_sort)
        Combined_sort <- sort(Combined, decreasing = FALSE)

        LPM_x_sort <- LPM(1, Combined_sort, x)
        LPM_y_sort <- LPM(1, Combined_sort, y)

        x.ssd.y <- any(LPM_x_sort > LPM_y_sort)

        ifelse(!x.ssd.y & min(x) >= min(y) & !identical(LPM_x_sort, LPM_y_sort),
               return(1),
               return(0))

    }
}


#' NNS TSD Test uni-directional
#'
#' Uni-directional test of third degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @return Returns (1) if \code{"X TSD Y"}, else (0).
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{https://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.TSD.uni(x, y)
#' @export


NNS.TSD.uni <- function(x, y){

    if(any(class(x)=="tbl")) x <- as.vector(unlist(x))
    if(any(class(y)=="tbl")) y <- as.vector(unlist(y))

    if(min(y) > min(x) | mean(y) > mean(x)) {
        return(0)
    } else {
        x_sort <- sort(x, decreasing = FALSE)
        y_sort <- sort(y, decreasing = FALSE)

        Combined <- c(x_sort, y_sort)
        Combined_sort <- sort(Combined, decreasing = FALSE)

        LPM_x_sort <- LPM(2, Combined_sort, x)
        LPM_y_sort <- LPM(2, Combined_sort, y)

        x.tsd.y <- any(LPM_x_sort > LPM_y_sort)

        ifelse(!x.tsd.y & min(x) >= min(y) & !identical(LPM_x_sort, LPM_y_sort),
               return(1),
               return(0))

    }
}
