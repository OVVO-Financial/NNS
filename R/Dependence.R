#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param order integer; Controls the level of quadrant partitioning.  Defaults to \code{(order = NULL)}.  Errors can generally be rectified by setting \code{(order = 1)}.  Will not partition further if less than 4 observations exist in a quadrant.
#' @param degree integer; Defaults to NULL to allow number of observations to be \code{"degree"} determinant.
#' @param print.map logical; \code{FALSE} (default) Plots quadrant means.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#' @keywords dependence, correlation
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
                   print.map = FALSE){


  if(is.null(degree)){degree = ifelse(length(x) < 100, 0, 1)
  } else {
    degree = degree
  }

  if(length(x) < 20){
    order = 1
    max.obs = 1
    }else{
      order = order
      max.obs = NULL
      }


  if(!missing(y)){

    if(print.map == T){
      part.map = NNS.part(x, y, order = order, max.obs.req = max.obs, Voronoi = TRUE, min.obs.stop = TRUE)
    }
    else {
      part.map = NNS.part(x, y, order = order, max.obs.req = max.obs, min.obs.stop = TRUE)
    }

    part.df = part.map$dt

    part.df[ , `:=` (mean.x = mean(x), mean.y = mean(y)), by = prior.quadrant]

    part.df[ , `:=` (sub.clpm = Co.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                    sub.cupm = Co.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                    sub.dlpm = D.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                    sub.dupm = D.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                    counts = .N
                    ), by = prior.quadrant]

### Re-run with order=1 if categorical data...
    if(part.df[ , sum(sub.clpm, sub.cupm, sub.dlpm, sub.dupm)] == 0){

      mode = function(x) {
        if(length(x) > 1){
          d <- density(x)
          d$x[which.max(d$y)]
        } else {
          x
        }
      }

      if(print.map == TRUE){
        part.map = NNS.part(x, y, order = 1, Voronoi = TRUE, noise.reduction = 'mode')
      } else {
        part.map = NNS.part(x, y, order = 1, noise.reduction = 'mode')
      }

      part.df = part.map$dt

      part.df[ , `:=` (mean.x = mode(x), mean.y = mode(y)), by = prior.quadrant]

      part.df[ , `:=` (sub.clpm = Co.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                      sub.cupm = Co.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                      sub.dlpm = D.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                      sub.dupm = D.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                      counts = .N
                      ), by = prior.quadrant]

      part.df[ , c("x", "y", "quadrant", "mean.x", "mean.y") := NULL]

      setkey(part.df,prior.quadrant)
      part.df = unique(part.df[])

      part.df[ , `:=` (nns.cor = (sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm),
                       nns.dep = abs(sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm))]

      part.df = part.df[counts == 1, counts := 0]
      part.df = part.df[(sub.clpm == 0 & sub.cupm == 0 & sub.dlpm == 0 & sub.dupm == 0), counts := 0]
      zeros = length(x) - part.df[ , sum(counts)]

      part.df = part.df[ , `:=` (weight = counts / (length(x) - zeros)), by = prior.quadrant]

  } else {#Categorical re-run

      part.df[ ,c("x", "y", "quadrant", "mean.x", "mean.y") := NULL]

      setkey(part.df,prior.quadrant)
      part.df = unique(part.df[])

      part.df[ , `:=` (nns.cor = (sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm),
                     nns.dep = abs(sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm))]


    part.df = part.df[counts == 1, counts := 0]
    part.df = part.df[(sub.clpm == 0 & sub.cupm == 0 & sub.dlpm == 0 & sub.dupm == 0), counts := 0]
    zeros = length(x) - part.df[ , sum(counts)]

    part.df = part.df[ , `:=` (weight = counts / (length(x) - zeros)), by = prior.quadrant]

}


      for (j in seq_len(ncol(part.df))){
        set(part.df, which(is.na(part.df[[j]])), j, 0)}

      nns.cor = part.df[ , sum(nns.cor = weight * nns.cor)]
      nns.dep = part.df[ , sum(nns.dep = weight * nns.dep)]

    return(list("Correlation" = nns.cor,
                "Dependence" = nns.dep))

  }#Not missing Y

  else{
    NNS.dep.matrix(x, order = order, degree = degree)
  }

}
