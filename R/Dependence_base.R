#' NNS Dependence Base
#'
#' Internal function for NNS dependence \link{NNS.dep} parallel instances.
#' @param x from \link{NNS.part}
#' @param y from \link{NNS.part}
#' @param order from \link{NNS.part}
#' @param degree from \link{NNS.part}
#' @param type from \link{NNS.part}
#' @param print.map from \link{NNS.part}
#' @param asym for asymmetrical dependecies
#'
#' @return Returns NNS dependence.
#'
#'
#' @export

NNS.dep.base <- function(x,
                   y = NULL,
                   order = NULL,
                   degree = NULL,
                   type = NULL,
                   print.map = FALSE,
                   asym = FALSE){

  oldw <- getOption("warn")
  options(warn = -1)

  if(!missing(y)) {
      if(length(x) < 20 | class(x) == "factor" | class(y) == "factor") {
          order <- 2
          y <- as.numeric(y)
          asym <- TRUE
      } else {
          order <- order
      }

      if(is.null(degree)) {
          degree <- ifelse(length(x) <= 100, 0, 1)
      }
  } else {
      if(length(x[, 1]) < 20 | any(unlist(lapply(x, is.factor)))) {
          order <- 2
          x <- data.matrix(x)
          asym <- TRUE
      }

      if(is.null(degree)) {
          degree <- ifelse(length(x[, 1]) <= 100, 0, 1)
      }
  }

  if(length(unique(x)) < sqrt(length(x)) || length(unique(y)) < sqrt(length(y))) {
      order <- 1
  }

  if(!missing(y)) {
      n = length(x)
      if(asym){type <- "XONLY"} else {type <- NULL}
      if (print.map == TRUE) {
          part.map <- NNS.part(x, y, order = order, obs.req = 1, type = type,
                               Voronoi = TRUE, min.obs.stop = TRUE)
      } else {
          part.map <- NNS.part(x, y, order = order, obs.req = 1, type = type,
                               Voronoi = FALSE, min.obs.stop = TRUE)
      }

      part.df <- part.map$dt

      if(any(length(unique(x)) < sqrt(length(x)) | length(unique(y)) < sqrt(length(y))  | is.na(sd(x)) | is.na(sd(y)) | sd(x)==0 | sd(y)==0)){
            part.df[, `:=`(mean.x = gravity(x), mean.y = gravity(y)), by = prior.quadrant]
        if (degree == 0) {
            part.df <- part.df[x != mean.x & y != mean.y, ]
        }

        part.df[, `:=`(sub.clpm = Co.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                       sub.cupm = Co.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                       sub.dlpm = D.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                       sub.dupm = D.UPM(degree, degree, x, y, mean.x[1], mean.y[1]), counts = .N),
                by = prior.quadrant]

        part.df[, `:=`(c("x", "y", "quadrant", "mean.x", "mean.y"), NULL)]

        setkey(part.df, prior.quadrant)

        part.df <- unique(part.df[])

        part.df[, `:=`(nns.cor = (sub.clpm + sub.cupm - sub.dlpm - sub.dupm)/(sub.clpm + sub.cupm + sub.dlpm + sub.dupm),
                       nns.dep = abs(sub.clpm + sub.cupm - sub.dlpm - sub.dupm)/(sub.clpm + sub.cupm + sub.dlpm + sub.dupm))]

        part.df <- part.df[counts == 1, `:=`(counts, 0)]

        part.df <- part.df[(sub.clpm == 0 & sub.cupm == 0 & sub.dlpm == 0 & sub.dupm == 0), `:=`(counts, 0)]

        zeros <- length(x) - part.df[, sum(counts)]

        part.df <- part.df[, `:=`(weight = counts/(length(x) - zeros)), by = prior.quadrant]

        for (j in seq_len(ncol(part.df))) {
          set(part.df, which(is.na(part.df[[j]])), j, 0)
        }

        nns.cor <- part.df[, sum(nns.cor = weight * nns.cor)]
        nns.dep <- part.df[, sum(nns.dep = weight * nns.dep)]

        return(list(Correlation = nns.cor, Dependence = nns.dep))

      } else {
          part.df[, `:=` (weight = .N/n), by = prior.quadrant]

          if(asym){
              disp <- part.df[,.(cor(x, abs(y), method = "pearson")), by = prior.quadrant]$V1
          } else {
              disp <- part.df[,.(cor(x, y, method = "pearson")), by = prior.quadrant]$V1
          }

          disp[is.na(disp)] <- 0


          nns.cor <- sum(disp * part.df[, mean(weight), by = prior.quadrant]$V1)
          nns.dep <- sum(abs(disp) * part.df[, mean(weight), by = prior.quadrant]$V1)

          options(warn = oldw)
          return(list(Correlation = nns.cor, Dependence = nns.dep))
      }
  } else {
        NNS.dep.matrix(x, order = order, degree = degree, asym = asym)
  }
}
