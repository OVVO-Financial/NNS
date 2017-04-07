#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param order integer; Controls the level of quadrant partitioning.  Defaults to \code{(order=3)}.  Errors can generally be rectified by setting \code{(order=1)}.  Will not partition further if less than 4 observations exist in a quadrant.
#' @param degree integer; Defaults to NULL to allow number of observations to be \code{"degree"} determinant.
#' @param print.map  logical; \code{FALSE} (default) Plots quadrant means.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#' @keywords dependence, correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.dep(x,y)
#'
#' ## Correlation / Dependence Matrix
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' B<-cbind(x,y,z)
#' NNS.dep(B)
#' @export

NNS.dep = function(x,y=NULL,order = 3,
                   degree=NULL,
                   print.map=FALSE){


  if(is.null(degree)){degree=ifelse(length(x)<100,0,1)}else{degree=degree}
  if(length(x)<20){
    order=1
    }else{
      order=order
      }


  if(!missing(y)){

    if(print.map==T){
      part.map = NNS.part(x,y,order=order, Voronoi=T)
    }
    else {
      part.map = NNS.part(x,y,order=order)
    }

    part.df = part.map$dt

    part.df[, `:=` (sub.clpm=Co.LPM(degree,degree,x,y,mean(x),mean(y)),
                    sub.cupm=Co.UPM(degree,degree,x,y,mean(x),mean(y)),
                    sub.dlpm=D.LPM(degree,degree,x,y,mean(x),mean(y)),
                    sub.dupm=D.UPM(degree,degree,x,y,mean(x),mean(y)),
                    counts=.N
                    ),by=prior.quadrant]

    setkey(part.df,prior.quadrant)
    part.df=(unique(part.df[,.(prior.quadrant,sub.clpm,sub.cupm,sub.dlpm,sub.dupm,
    nns.cor=(sub.clpm+sub.cupm-sub.dlpm-sub.dupm)/(sub.clpm+sub.cupm+sub.dlpm+sub.dupm),
    nns.dep=abs(sub.clpm+sub.cupm-sub.dlpm-sub.dupm)/(sub.clpm+sub.cupm+sub.dlpm+sub.dupm), counts)]))

    zeros=part.df[,sum(counts==1)]

    part.df=part.df[,`:=` (weight=counts/(length(x)-zeros)),by=prior.quadrant]

    part.df=part.df[counts==1, weight := 0]


### Re-run with order=1 if categorical data...
    if(part.df[,sum(sub.clpm,sub.cupm,sub.dlpm,sub.dupm)]==0){
      if(print.map==T){
        part.map = NNS.part(x,y,order=1, Voronoi=T)}
        else {
        part.map = NNS.part(x,y,order=1)
        }

      part.df = part.map$dt

      part.df[, `:=` (sub.clpm=Co.LPM(degree,degree,x,y,mean(x),mean(y)),
                      sub.cupm=Co.UPM(degree,degree,x,y,mean(x),mean(y)),
                      sub.dlpm=D.LPM(degree,degree,x,y,mean(x),mean(y)),
                      sub.dupm=D.UPM(degree,degree,x,y,mean(x),mean(y)),
                      counts=.N),by=prior.quadrant]

      setkey(part.df,prior.quadrant)
      part.df=(unique(part.df[,.(prior.quadrant,sub.clpm,sub.cupm,sub.dlpm,sub.dupm,
               nns.cor=(sub.clpm+sub.cupm-sub.dlpm-sub.dupm)/(sub.clpm+sub.cupm+sub.dlpm+sub.dupm),
               nns.dep=abs(sub.clpm+sub.cupm-sub.dlpm-sub.dupm)/(sub.clpm+sub.cupm+sub.dlpm+sub.dupm), counts)]))

      zeros=part.df[,sum(counts==1)]

      part.df=part.df[,`:=` (weight=counts/(length(x)-zeros)),by=prior.quadrant]

      part.df=part.df[counts==1, weight := 0]
} #Categorical re-run


      for (j in seq_len(ncol(part.df))){
        set(part.df,which(is.na(part.df[[j]])),j,0)}

      nns.cor=part.df[,sum(nns.cor=weight*nns.cor)]
      nns.dep=part.df[,sum(nns.dep=weight*nns.dep)]

    return(list("Correlation"=nns.cor,"Dependence"= nns.dep ))

  }#Not missing Y

  else{
    NNS.dep.matrix(x,order=order,degree = degree)
  }

}
