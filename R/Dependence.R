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
#' @note p-values and confidence intervals can be obtained from sampling random permutations of \code{y_p} and running \code{NNS.dep(x,y_p)} to compare against a null hypothesis of 0 correlation or dependence between \code{x,y}.
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
#'
#' \dontrun{
#' ## p-values for [NNS.dep]
#' x=seq(-5,5,.1);y=x^2+rnorm(length(x))
#'
#'
#' nns_cor_dep = NNS.dep(x,y,print.map = TRUE)
#' nns_cor_dep
#'
#' ## Create permutations of y
#' y_p = replicate(1000,sample.int(length(y)))
#'
#' ## Generate new correlation and dependence measures on each new permutation of y
#' nns.mc = apply(y_p,2,function(g) NNS.dep(x,y[g]))
#'
#' ## Store results
#' cors = unlist(lapply(nns.mc, "[[",1))
#' deps = unlist(lapply(nns.mc, "[[",2))
#'
#' ## View results
#' hist(cors)
#' hist(deps)
#'
#' ## Left tailed correlation p-value
#' cor_p_value = LPM(0,nns_cor_dep$Correlation,cors)
#' cor_p_value
#'
#' ## Right tailed correlation p-value
#' cor_p_value = UPM(0,nns_cor_dep$Correlation,cors)
#' cor_p_value
#'
#' ## Two sided correlation p-value
#' cor_p_value = UPM(0,abs(nns_cor_dep$Correlation),abs(cors))
#' cor_p_value
#'
#' ## Left tailed dependence p-value
#' dep_p_value = LPM(0,nns_cor_dep$Dependence,deps)
#' dep_p_value
#'
#' ## Right tailed dependence p-value
#' dep_p_value = UPM(0,nns_cor_dep$Dependence,deps)
#' dep_p_value
#'
#' ## Confidence Intervals
#' ## For 95th percentile VaR (both-tails) see [LPM.VaR] and [UPM.VaR]
#' ## Lower CI
#' LPM.VaR(.975,0,cors)
#' ## Upper CI
#' UPM.VaR(.975,0,cors)
#' }
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
            if(print.map){
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

        if(l <= 500){
            return(NNS.dep.base(x, y, order = order, degree = degree, print.map = print.map))
        }

        seg <- as.integer(.2*l)
        segs <- list(5L)

        for(i in 1:5){
            if(i == 1){
                segs[[i]] <- 1 : min(l,100)
            }
            if(i > 1 & i < 5){
                segs[[i]] <- max(1, (i*seg - 50)) : min(l,(i*seg+50))
            }
            if(i == 5){
                segs[[i]] <- max(1, (l-100)) : max(l,100)
            }
        }

        nns.dep <- list(5L)

cl <- makeCluster(num_cores)
registerDoParallel(cl)

  nns.dep <- foreach(i = 1:5,.packages = "NNS")%dopar%{
        NNS.dep.base(x[segs[[i]]], y[segs[[i]]], print.map = FALSE)
  }

stopCluster(cl)
registerDoSEQ()

        if(l > 500 & print.map){
            NNS.part(x, y, order = order, Voronoi = TRUE)
        }

        return(list("Correlation" = mean(unlist(lapply(nns.dep, `[[`, 1))),
                "Dependence" = mean(unlist(lapply(nns.dep, `[[`, 2)))))

    } else {
        return(NNS.dep.matrix(x, order=order, degree = degree))
    }

}
