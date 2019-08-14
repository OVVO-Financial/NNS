#' NNS Dependence Base
#'
#' Internal function for NNS dependence \link{NNS.dep} parallel instances.
#' @param x from \link{NNS.part}
#' @param y from \link{NNS.part}
#' @param order from \link{NNS.part}
#' @param degree from \link{NNS.part}
#' @param print.map from \link{NNS.part}
#'
#' @return Returns NNS dependence.
#'
#'
#' @export

NNS.dep.base <- function(x,
                   y = NULL,
                   order = NULL,
                   degree = NULL,
                   print.map = FALSE){

  if(!missing(y)){
      if(length(x) < 20 | class(x)=='factor' | class(y)=='factor'){
          order <- 1
          max.obs <- 1
      } else {
          order <- order
          max.obs <- NULL
      }
      if(is.null(degree)){
          degree <- ifelse(length(x) <= 100, 0, 1)
      } else {
          degree <- degree
      }

  } else {
      if(length(x[,1]) < 20 | any(unlist(lapply(x,is.factor)))){
          order <- 1
          max.obs <- 1
      } else {
          order <- order
          max.obs <- NULL
      }
      if(is.null(degree)){
          degree <- ifelse(length(x[,1]) <= 100, 0, 1)
      } else {
          degree <- degree
      }
  }


  if(length(unique(x))<=2){order <- 1}

  if(!missing(y)){
      if(is.null(order)){
          test.map <- NNS.part(x,y,order = order, obs.req = max.obs, Voronoi = FALSE, min.obs.stop = TRUE)$order
          if(test.map==1){
              noise.reduction <- 'median'
          } else{
              noise.reduction <- 'mean'
          }

      } else {
          if(order > 1){
              test.map <- NNS.part(x,y,order = order, obs.req = max.obs, Voronoi = FALSE, min.obs.stop = TRUE)$order
              if(test.map==1){
                noise.reduction <- 'median'
              } else{
                noise.reduction <- 'mean'
              }
          }
        if(order == 1){
              test.map <- NNS.part(x,y,order = order, obs.req = max.obs, Voronoi = FALSE, min.obs.stop = TRUE)$order
          if(test.map==1){
              noise.reduction <- 'median'
          } else{
              noise.reduction <- 'median'
          }
        }
      }
      if(print.map == TRUE){
          part.map <- NNS.part(x, y, order = order, obs.req = max.obs, Voronoi = TRUE, min.obs.stop = TRUE, noise.reduction = noise.reduction)
      } else {
          part.map <- NNS.part(x, y, order = order, obs.req = max.obs, min.obs.stop = TRUE, Voronoi = FALSE, noise.reduction = noise.reduction)
      }

      part.df <- part.map$dt

      part.df[ , `:=` (mean.x = mean(x), mean.y = mean(y)), by = prior.quadrant]

      if(degree==0){
          part.df <- part.df[x!=mean.x & y!=mean.y,]
      }


      part.df[ , `:=` (sub.clpm = Co.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                     sub.cupm = Co.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                     sub.dlpm = D.LPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                     sub.dupm = D.UPM(degree, degree, x, y, mean.x[1], mean.y[1]),
                     counts = .N
                      ), by = prior.quadrant]


      part.df[ ,c("x", "y", "quadrant", "mean.x", "mean.y") := NULL]

      setkey(part.df,prior.quadrant)
      part.df <- unique(part.df[])

      part.df[ , `:=` (nns.cor = (sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm),
                     nns.dep = abs(sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm))]


      part.df <- part.df[counts == 1, counts := 0]
      part.df <- part.df[(sub.clpm == 0 & sub.cupm == 0 & sub.dlpm == 0 & sub.dupm == 0), counts := 0]

      zeros <- length(x) - part.df[ , sum(counts)]

      part.df <- part.df[ , `:=` (weight = counts / (length(x) - zeros)), by = prior.quadrant]




    for (j in seq_len(ncol(part.df))){
        set(part.df, which(is.na(part.df[[j]])), j, 0)
    }

        nns.cor <- part.df[ , sum(nns.cor = weight * nns.cor)]
        nns.dep <- part.df[ , sum(nns.dep = weight * nns.dep)]

        return(list("Correlation" = nns.cor,
                    "Dependence" = nns.dep))

  }#Not missing Y
    else {
        NNS.dep.matrix(x, order = order, degree = degree)
  }

}
