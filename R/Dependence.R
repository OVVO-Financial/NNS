#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param order integer; Controls the level of quadrant partitioning.  Defaults to \code{(order = NULL)}.  Errors can generally be rectified by setting \code{(order = 1)}.  Will not partition further if less than 4 observations exist in a quadrant.
#' @param degree integer; Defaults to NULL to allow number of observations to be \code{"degree"} determinant.
#' @param print.map logical; \code{FALSE} (default) Plots quadrant means.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.dep(x, y)
#'
#' ## Correlation / Dependence Matrix
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' B <- cbind(x, y, z)
#' NNS.dep(B)
#' @export

NNS.dep = function(x,
                   y = NULL,
                   order = NULL,
                   degree = NULL,
                   print.map = FALSE,
                   ncores = NULL){

    if(!is.null(y)){
        # No dependence if only a single value
        if(length(unique(x))==1 | length(unique(y))==1){
            if(print.map==TRUE){
                NNS.part(x, y, order=1, Voronoi = TRUE)
            }

            return(list("Correlation" = 0,
                       "Dependence" = 0))
        }

 if (is.null(ncores)) {
        num_cores <- detectCores() - 1
      } else {
        num_cores <- ncores
      }

        l <- length(x)

        if(l <= 500){return(NNS.dep.base(x,y,order=order, degree = degree,print.map = print.map))}

        seg <- as.integer(.2*l)
        segs <- list(5L)

        for (i in 1:5){
            if(i == 1){
                segs[[i]] <- 1:min(l,100)
            }
            if(i > 1 & i < 5){
                segs[[i]] <- max(1,(i*seg - 50)):min(l,(i*seg+50))
            }
            if(i == 5){
                segs[[i]] <- max(1,(l-100)):max(l,100)
            }
        }

        nns.dep <- list(5L)

cl <- makeCluster(num_cores)
registerDoParallel(cl)

  nns.dep <- foreach(i = 1:5,.packages = "NNS")%dopar%{
        NNS.dep.base(x[segs[[i]]],y[segs[[i]]],print.map = FALSE)$Dependence
  }

stopCluster(cl)
registerDoSEQ()


        nns.cor <- NNS.dep.base(x,y,order=order, degree = degree, print.map = print.map)$Correlation

        return(list("Correlation" = nns.cor,
                "Dependence" = mean(unlist(nns.dep))))

    } else {
        return(NNS.dep.matrix(x, order=order, degree = degree))
    }

}
