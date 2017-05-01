#' NNS Partition Map
#'
#'  Creates partitions based on partial moment quadrant means, iteratively assigning identifications to observations based on those quadrants (unsupervised partitional and hierarchial clustering method).  Basis for correlation \link{NNS.cor}, dependence \link{NNS.dep}, regression \link{NNS.reg} routines.
#' @param x a numeric vector.
#' @param y a numeric vector with compatible dimsensions to \code{x}.
#' @param Voronoi logical; \code{FALSE} (default) Displays a Voronoi type diagram using partial moment quadrants.
#' @param type \code{NULL} (default) Controls the partitioning basis.  Set to \code{(type="XONLY")} for X-axis based partitioning.  Defaults to NULL for both X and Y-axis partitioning.
#' @param order integer; Number of partial moment quadrants to be generated.  \code{(order="max")} will institute a perfect fit.
#' @param min.obs integer; (4 default) Desired observations per cluster where quadrants will not be further partitioned if observations are not greater than the entered value.  Reduces minimum number of necessary observations in a quadrant to 1 when \code{(min.obs=1)}.
#' @param noise.reduction the method of determing regression points options: ("mean","median","mode","off"); \code{(noise.reduction="median")} uses medians instead of means for partitions, while \code{(noise.reduction="mode")} uses modes instead of means for partitions.  Defaults to \code{(noise.reduction="mean")}, while \code{(noise.reduction="off")} will partition quadrant to a single observation for a given \code{(order=...)}.
#' @return Returns both a \link{data.table} \code{("dt")} of \code{x} and \code{y} observations with their partition assignment \code{"quadrant"} in the 3rd column and their prior partition assignment \code{"prior.quadrant"} in the 4th column; and the \link{data.table} of regression points \code{("regression.points")} for that given \code{(order=...)}.  Also returns the \code{"order"} of the final partition given \code{"min.obs"} constraints.
#' @keywords partitioning, cluster
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.part(x,y)
#'
#' ## Data.table of observations and partitions
#' NNS.part(x,y,order=1)$dt
#'
#' ## Regression points
#' NNS.part(x,y,order=1)$regression.points
#'
#' ## Voronoi style plot
#' NNS.part(x,y,Voronoi=TRUE)
#'
#' ## Examine final counts by quadrant
#' DT=NNS.part(x,y)$dt
#' DT[,counts := .N,by=quadrant]
#' DT
#' @export

NNS.part = function(x, y,Voronoi=FALSE,type=NULL,order= NULL,min.obs=4,noise.reduction="mean"){

  if(!is.null(order)){if(order==0) {order=1} else {order=order}}


  PART = data.table(x=as.numeric(x),y=as.numeric(y),quadrant='q',prior.quadrant='q',key = c('quadrant','x','y'))

  mode=function(x) {
    if(length(x)>1){
      d <- density(x)
      d$x[which.max(d$y)]
    }else{x}
  }

  if(Voronoi==T){
    PART[,plot(x,y,col='steelblue',cex.lab=2,xlab = "X",ylab="Y")]
  }

  if(is.null(order)){order=Inf}

  if(!is.numeric(order)){
    min.obs=1
    type=type
  }else{
    min.obs=min.obs
    type=type
  }


  if(noise.reduction=='off'){min.obs=1}else{min.obs=min.obs}

  ### X and Y partition
  if(is.null(type)){
    i=0L
    while(i>=0){
      if(i==order){break}

      PART = PART[,counts := .N,quadrant]
      setkey(PART,counts,quadrant,x,y)
      l.PART = sum(PART$counts>min.obs)

      if(l.PART==0){break}

      #Segments for Voronoi...
      if(Voronoi==T){
        if(l.PART>0){
          PART[ counts>=min.obs , segments(min(x),mean(y),max(x),mean(y),lty=3),by=quadrant]
          PART[ counts>=min.obs , segments(mean(x),min(y),mean(x),max(y),lty=3),by=quadrant]
        }
      }



      if(noise.reduction=='mean' | noise.reduction=='off'){
        PART[ counts>=min.obs , prior.quadrant := (quadrant)]
        PART[ counts>=min.obs , quadrant :=
                ifelse( x <= mean(x) & y <= mean(y) , paste0(quadrant,4),
                        ifelse(x <= mean(x) & y >mean(y),paste0(quadrant,2),
                               ifelse(x > mean(x) & y <=mean(y),paste0(quadrant,3),                                    paste0(quadrant,1)))) ,
              by='quadrant']

        RP = PART[, .(x=mean(x),y=mean(y)),by=prior.quadrant]
      }

      if(noise.reduction=='median'){
        PART[ counts>=min.obs , prior.quadrant := (quadrant)]
        PART[ counts>=min.obs , quadrant :=
                ifelse( x <= median(x) & y <= median(y) , paste0(quadrant,4),
                        ifelse(x <= median(x) & y >median(y),paste0(quadrant,2),
                               ifelse(x > median(x) & y <=median(y),paste0(quadrant,3),                                    paste0(quadrant,1)))), by='quadrant']
        RP = PART[ , .(x=median(x),y=median(y)),by=prior.quadrant]
      }

      if(noise.reduction=='mode'){
        PART[ counts>=min.obs , prior.quadrant := (quadrant)]
        PART[ counts>=min.obs , quadrant :=
                ifelse( x <= mode(x) & y <= mode(y) , paste0(quadrant,4),
                        ifelse(x <= mode(x) & y > mode(y),paste0(quadrant,2),
                               ifelse(x > mode(x) & y <=mode(y),paste0(quadrant,3),                                    paste0(quadrant,1)))), by='quadrant']
        RP = PART[ , .(x=mode(x),y=mode(y)),by=prior.quadrant]
      }

      if(!is.numeric(order)){
        RP.new = PART[,.(x=mean(x),y=mean(y)),by=quadrant]
        if(length(RP.new$x)==length(RP$x)){break} }


      if(length(RP$x)==length(x)){break}
      i=i+1L
    }

    if(Voronoi==T){
      points(RP$x,RP$y,pch=15,lwd=2,col='red')
      title(main=paste0("NNS Order = ",i),cex.main=2)

    }

    setnames(RP,"prior.quadrant","quadrant")

    return(list("order"=i,"dt"=PART[, .(x, y, quadrant,prior.quadrant)],"regression.points"=RP))

  }


  ### X ONLY partition
  if(!is.null(type)){
    i=0L
    while(i>=0){
      if(i==order){break}

      PART = PART[,counts := .N,quadrant]
      setkey(PART,counts,quadrant,x,y)
      l.PART = sum(PART$counts>min.obs)

      if(l.PART==0){break}

      if(noise.reduction=='mean' | noise.reduction=='off'){
        RP = PART[,.(x=mean(x),y=mean(y)),by=quadrant]
        PART[counts>=min.obs, prior.quadrant := (quadrant)]
        PART[counts>=min.obs, quadrant :=
               ifelse( x <= mean(x) , paste0(quadrant,1),paste0(quadrant,2)), by='quadrant']
        if(!is.numeric(order)){
          RP.new = PART[,.(x=mean(x),y=mean(y)),by=quadrant]
          if(length(RP.new$x)==length(RP$x)){break}}
      }

      else {
        if(noise.reduction=='mode'){
          RP = PART[,.(x=mode(x),y=mode(y)),by=quadrant]
          PART[counts>=min.obs, prior.quadrant := (quadrant)]
          PART[counts>=min.obs, quadrant :=
                 ifelse( x <= mode(x), paste0(quadrant,1),paste0(quadrant,2)), by='quadrant']
          if(!is.numeric(order)){
            RP.new = PART[,.(x=mode(x),y=mode(y)),by=quadrant]
            if(length(RP.new$x)==length(RP$x)){break}}
        }

        else {
          if(noise.reduction=='median'){
            RP = PART[,.(x=median(x),y=median(y)),by=quadrant]
            PART[counts>=min.obs, prior.quadrant := (quadrant)]
            PART[counts>=min.obs, quadrant :=
                   ifelse( x <= median(x), paste0(quadrant,1),paste0(quadrant,2)), by='quadrant']
            if(!is.numeric(order)){
              RP.new = PART[,.(x=median(x),y=median(y)),by=quadrant]
              if(length(RP.new$x)==length(RP$x)){break}}
          }
        }}

      if(length(RP$x)==length(x)){break}
      i=i+1L
    }


    if(Voronoi==T){
      plot(x,y,col='steelblue',cex.lab=2,xlab = "X",ylab="Y")
      abline(v=RP$x,lty=3)
      points(RP$x,RP$y,pch=15,lwd=2,col='red')
      title(main=paste0("NNS Order = ",i),cex.main=2)
    }

    return(list("order"=i,"dt"=PART[, .(x, y, quadrant,prior.quadrant)],"regression.points"=RP))
  }

}
