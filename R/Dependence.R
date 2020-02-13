#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param order integer; Controls the level of quadrant partitioning.  Defaults to \code{(order = 3)}.  Errors can generally be rectified by setting \code{(order = 1)}.
#' @param degree integer; Defaults to NULL to allow number of observations to be \code{"degree"} determinant.
#' @param asym logical; \code{FALSE} (default) Allows for asymmetrical dependencies.
#' @param print.map logical; \code{FALSE} (default) Plots quadrant means.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#' @note p-values and confidence intervals can be obtained from sampling random permutations of \code{y_p} and running \code{NNS.dep(x,y_p)} to compare against a null hypothesis of 0 correlation or independence between \code{x,y}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.dep(x, y)
#'
#' ## Correlation / Dependence Matrix
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' B <- cbind(x, y, z)
#' NNS.dep(B)
#'
#'
#' ## p-values for [NNS.dep]
#' x <- seq(-5, 5, .1); y <- x^2 + rnorm(length(x))
#'
#'
#' nns_cor_dep <- NNS.dep(x, y, print.map = TRUE)
#' nns_cor_dep
#'
#' ## Create permutations of y
#' y_p <- replicate(1000, sample.int(length(y)))
#'
#' ## Generate new correlation and dependence measures on each new permutation of y
#' nns.mc <- apply(y_p, 2, function(g) NNS.dep(x, y[g]))
#'
#' ## Store results
#' cors <- unlist(lapply(nns.mc, "[[", 1))
#' deps <- unlist(lapply(nns.mc, "[[", 2))
#'
#' ## View results
#' hist(cors)
#' hist(deps)
#'
#' ## Left tailed correlation p-value
#' cor_p_value <- LPM(0, nns_cor_dep$Correlation, cors)
#' cor_p_value
#'
#' ## Right tailed correlation p-value
#' cor_p_value <- UPM(0, nns_cor_dep$Correlation, cors)
#' cor_p_value
#'
#' ## Confidence Intervals
#' ## For 95th percentile VaR (both-tails) see [LPM.VaR] and [UPM.VaR]
#' ## Lower CI
#' LPM.VaR(.975, 0, cors)
#' ## Upper CI
#' UPM.VaR(.975, 0, cors)
#'
#' ## Left tailed dependence p-value
#' dep_p_value <- LPM(0, nns_cor_dep$Dependence, deps)
#' dep_p_value
#'
#' ## Right tailed dependence p-value
#' dep_p_value <- UPM(0, nns_cor_dep$Dependence, deps)
#' dep_p_value
#' }
#' @export

NNS.dep = function(x,
                   y = NULL,
                   order = 3,
                   degree = NULL,
                   asym = FALSE,
                   print.map = FALSE,
                   ncores = NULL){

  oldw <- getOption("warn")
  options(warn = -1)

  if(asym){type <- "XONLY"} else {type <- NULL}

  if(!is.null(y)){
    # No dependence if only a single value
    if(length(unique(x))==1 | length(unique(y))==1){
      if(print.map){
        NNS.part(x, y, order = 1, Voronoi = TRUE, type = type)
      }

      options(warn = oldw)
      return(list("Correlation" = 0,
                  "Dependence" = 0))
    }


    if (is.null(ncores)) {
      num_cores <- detectCores() - 1
    } else {
      num_cores <- ncores
    }

    l <- length(x)

    if(l < 150){
        return(NNS.dep.base(x, y, order = order, degree = degree, print.map = print.map, asym = asym, type = "XONLY"))
    }

    seg <- as.integer(.2*l)
    segs <- list(5L)
    uniques <- list(5L)

    for(i in 1:5){
      if(i == 1){
        segs[[i]] <- 1 : min(l, min(100, l/2.5))
        uniques[[i]] <- length(unique(x[segs[[i]]]))
      }
      if(i > 1 & i < 5){
        segs[[i]] <- max(1, (i*seg - min(50, l/5))) : min(l,(i*seg + min(50, l/5)))
        uniques[[i]] <- length(unique(x[segs[[i]]]))
      }
      if(i == 5){
        segs[[i]] <- max(1, (l - min(100, l/2.5))) : l
        uniques[[i]] <- length(unique(x[segs[[i]]]))
      }
    }

    nns.dep <- list(5L)

    if(any(unlist(uniques)==1)){
      DT <- data.table(x, y)
      setkey(DT[, x := x], x)

      for(i in 1:3){
        if(i==1){
          nns.dep[[i]] <- NNS.dep.base(DT[, .SD[1], by="x"]$x, DT[, .SD[1], by = "x"]$y, print.map = FALSE, order = order, asym = asym)
        }
        if(i==2) {
          nns.dep[[i]] <- NNS.dep.base(DT[, .SD[min(1,round(.N/2))], by="x"]$x, DT[, .SD[min(1,round(.N/2))], by = "x"]$y, print.map = FALSE, order = order, asym = asym)
        }
        if(i==3) {
          nns.dep[[i]] <- NNS.dep.base(DT[, .SD[.N], by = "x"]$x, DT[, .SD[.N], by = "x"]$y, print.map = FALSE, order = order, asym = asym)
        }
      }

    } else {
        for(i in 1:5){
            nns.dep[[i]] <- NNS.dep.base(x[segs[[i]]], y[segs[[i]]], print.map = FALSE, order = order, asym = asym)
        }

    }


    if(l >= 150 & print.map){
      NNS.part(x, y, order = order, min.obs.stop = TRUE, Voronoi = TRUE, type = type)
    }
    options(warn = oldw)
    return(list("Correlation" = mean(unlist(lapply(nns.dep, `[[`, 1))),
                "Dependence" = mean(unlist(lapply(nns.dep, `[[`, 2)))))

  } else {
    return(NNS.dep.matrix(x, order = order, degree = degree, asym = asym))
  }

}
