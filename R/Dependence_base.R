NNS.dep.base = function(x,
                   y = NULL,
                   order = NULL,
                   degree = NULL,
                   print.map = FALSE){

  if(!missing(y)){
    if(length(x) < 20 | class(x)=='factor' | class(y)=='factor'){
      order = 1
      max.obs = 1
    } else {
      order = order
      max.obs = NULL
    }
    if(is.null(degree)){degree = ifelse(length(x) <= 100, 0, 1)
    } else {
      degree = degree
    }

  } else {
    if(length(x[,1]) < 20 | any(unlist(lapply(x,is.factor)))){
      order = 1
      max.obs = 1
    } else {
      order = order
      max.obs = NULL
    }
    if(is.null(degree)){degree = ifelse(length(x[,1]) <= 100, 0, 1)
    } else {
      degree = degree
    }
  }


  if(length(unique(x))<=2){order = 1}

  if(!missing(y)){

    if(print.map == TRUE){
      part.map = NNS.part(x, y, order = order, max.obs.req = max.obs, Voronoi = TRUE, min.obs.stop = TRUE)
    }
    else {
      part.map = NNS.part(x, y, order = order, max.obs.req = max.obs, min.obs.stop = TRUE, Voronoi = FALSE)
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

    part.df[ ,c("x", "y", "quadrant", "mean.x", "mean.y") := NULL]

    setkey(part.df,prior.quadrant)
    part.df = unique(part.df[])

    part.df[ , `:=` (nns.cor = (sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm),
                     nns.dep = abs(sub.clpm + sub.cupm - sub.dlpm - sub.dupm) / (sub.clpm + sub.cupm + sub.dlpm + sub.dupm))]


    part.df = part.df[counts == 1, counts := 0]
    part.df = part.df[(sub.clpm == 0 & sub.cupm == 0 & sub.dlpm == 0 & sub.dupm == 0), counts := 0]
    zeros = length(x) - part.df[ , sum(counts)]

    part.df = part.df[ , `:=` (weight = counts / (length(x) - zeros)), by = prior.quadrant]




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
