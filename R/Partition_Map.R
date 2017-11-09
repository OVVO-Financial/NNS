#' NNS Partition Map
#'
#' Creates partitions based on partial moment quadrant means, iteratively assigning identifications to observations based on those quadrants (unsupervised partitional and hierarchial clustering method).  Basis for correlation \link{NNS.cor}, dependence \link{NNS.dep}, regression \link{NNS.reg} routines.
#' @param x a numeric vector.
#' @param y a numeric vector with compatible dimsensions to \code{x}.
#' @param Voronoi logical; \code{FALSE} (default) Displays a Voronoi type diagram using partial moment quadrants.
#' @param type \code{NULL} (default) Controls the partitioning basis.  Set to \code{(type="XONLY")} for X-axis based partitioning.  Defaults to \code{NULL} for both X and Y-axis partitioning.
#' @param order integer; Number of partial moment quadrants to be generated.  \code{(order="max")} will institute a perfect fit.
#' @param max.obs integer; (4 default) Desired observations per cluster where quadrants will not be further partitioned if observations are not greater than the entered value.  Reduces minimum number of necessary observations in a quadrant to 1 when \code{(max.obs=1)}.
#' @param min.obs.stop logical; \code{FALSE} (default) Stopping condition where quadrants will not be further partitioned if a single cluster contains less than the entered value of \code{max.obs}.
#' @param noise.reduction the method of determing regression points options: ("mean","median","mode","off"); \code{(noise.reduction="median")} uses medians instead of means for partitions, while \code{(noise.reduction="mode")} uses modes instead of means for partitions.  Defaults to \code{(noise.reduction="mean")}, while \code{(noise.reduction="off")} will partition quadrant to a single observation for a given \code{(order=...)}.
#' @return Returns:
#'  \itemize{
#'   \item{\code{"dt"}} a \link{data.table} of \code{x} and \code{y} observations with their partition assignment \code{"quadrant"} in the 3rd column and their prior partition assignment \code{"prior.quadrant"} in the 4th column.
#'   \item{\code{"regression.points"}} the \link{data.table} of regression points for that given \code{(order=...)}.
#'   \item{\code{"order"}}  the \code{order} of the final partition given \code{"min.obs.stop"} stopping condition.
#'   }
#'
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

NNS.part = function(x, y,Voronoi=FALSE,type=NULL,order= NULL,max.obs=4,min.obs.stop=FALSE,noise.reduction="mean"){
  if(is.null(max.obs)) max.obs=4
  if(!is.null(order)){if(order==0) {order=1} else {order=order}}
  x=as.numeric(x);y=as.numeric(y)

  PART = data.table(x,y,quadrant="q",prior.quadrant="pq")[, counts := .N]

  mode=function(x) {
    if(length(x)>1){
      d <- density(x)
      d$x[which.max(d$y)]
    }else{x}
  }

  if(Voronoi==T){
    plot(x,y,col='steelblue',cex.lab=2,xlab = "X",ylab="Y")
  }

  if(is.null(order)){order=Inf}

  if(!is.numeric(order)){
    max.obs=1
    type=type
  }else{
    max.obs=max.obs
    type=type
  }


  if(noise.reduction=='off') max.obs=1 else max.obs=max.obs

  ### X and Y partition
  if(is.null(type)){
    i=0L
    while(i>=0){
      if(i==order) break

      PART[counts >= max.obs, counts := .N, by=quadrant]
      l.PART=max(PART$counts)
      if(min.obs.stop && (min(PART$counts)<= max.obs)) break
      if (l.PART <= max.obs) break
      max.obs.rows = PART[counts >= max.obs, which=TRUE]
      #Segments for Voronoi...
      if(Voronoi==T){
        if(l.PART>max.obs){
          if(noise.reduction=='mean' | noise.reduction=='off'){
          PART[ max.obs.rows ,{
            segments(min(x),mean(y),max(x),mean(y),lty=3)
            segments(mean(x),min(y),mean(x),max(y),lty=3)
          },by=quadrant]
          }
          if(noise.reduction=='median'){
          PART[ max.obs.rows ,{
            segments(min(x),median(y),max(x),median(y),lty=3)
            segments(median(x),min(y),median(x),max(y),lty=3)
          },by=quadrant]
          }
          if(noise.reduction=='mode'){
          PART[ max.obs.rows ,{
            segments(min(x),mode(y),max(x),mode(y),lty=3)
            segments(mode(x),min(y),mode(x),max(y),lty=3)
          },by=quadrant]
          }
        }
      }


      if(noise.reduction=='mean' | noise.reduction=='off'){
        RP= PART[max.obs.rows, lapply(.SD, mean), by=quadrant, .SDcols = x:y]
        RP[, prior.quadrant := (quadrant)]

        PART[max.obs.rows , prior.quadrant := (quadrant)]
        old.parts=length(unique(PART$quadrant))
        PART[RP, on=.(quadrant), q_new := {
          lox = x.x <= i.x
          loy = x.y <= i.y
          1L + lox + loy*2L
        }]
        PART[max.obs.rows, quadrant := paste0(quadrant, q_new)]

        new.parts=length(unique(PART$quadrant))
      }

      if(noise.reduction=='median'){
        RP= PART[max.obs.rows, lapply(.SD, median), by=quadrant, .SDcols = x:y]
        RP[, prior.quadrant := (quadrant)]

        PART[max.obs.rows , prior.quadrant := (quadrant)]
        old.parts=length(unique(PART$quadrant))
        PART[RP, on=.(quadrant), q_new := {
          lox = x.x <= i.x
          loy = x.y <= i.y
          1L + lox + loy*2L
        }]
        PART[max.obs.rows, quadrant := paste0(quadrant, q_new)]

        new.parts=length(unique(PART$quadrant))
      }

      if(noise.reduction=='mode'){
        RP= PART[max.obs.rows, lapply(.SD, mode), by=quadrant, .SDcols = x:y]
        RP[, prior.quadrant := (quadrant)]

        PART[max.obs.rows , prior.quadrant := (quadrant)]
        old.parts=length(unique(PART$quadrant))
        PART[RP, on=.(quadrant), q_new := {
          lox = x.x <= i.x
          loy = x.y <= i.y
          1L + lox + loy*2L
        }]
        PART[max.obs.rows, quadrant := paste0(quadrant, q_new)]

        new.parts=length(unique(PART$quadrant))
      }

      if(old.parts==new.parts) break

      i=i+1L
    }

    if(Voronoi==T){
      points(RP$x,RP$y,pch=15,lwd=2,col='red')
      title(main=paste0("NNS Order = ",i),cex.main=2)
    }

    setnames(RP,"prior.quadrant","quadrant")
    PART[,`:=`(counts = NULL, q_new = NULL)]
    DT=setorder(PART[], quadrant,x, y)[]
    RP=setorder(RP[],quadrant)[]
    return(list("order"=i,"dt"=DT,"regression.points"=RP))

  }


  ### X ONLY partition
  if(!is.null(type)){
    i=0L
    while(i>=0){
      if(i==order) break

      PART[counts >= 2*max.obs, counts := .N, by=quadrant]

      if (max(PART$counts) <= 2*max.obs) break
      if(min.obs.stop && (min(PART$counts)<= 2*max.obs)) break
      max.obs.rows = PART[counts >= 2*max.obs, which=TRUE]


      if(noise.reduction=='mean' | noise.reduction=='off'){
        RP= PART[max.obs.rows, lapply(.SD, mean), by=quadrant, .SDcols = x:y]
        RP[, prior.quadrant := (quadrant)]

        PART[max.obs.rows , prior.quadrant := (quadrant)]
        old.parts=length(unique(PART$quadrant))
        PART[RP, on=.(quadrant), q_new := {
          lox = x.x > i.x
          1L + lox
        }]
        PART[max.obs.rows, quadrant := paste0(quadrant, q_new)]

        new.parts=length(unique(PART$quadrant))
      }

      if(noise.reduction=='mode'){
        RP= PART[max.obs.rows, lapply(.SD, mode), by=quadrant, .SDcols = x:y]
        RP[, prior.quadrant := (quadrant)]

        PART[max.obs.rows , prior.quadrant := (quadrant)]
        old.parts=length(unique(PART$quadrant))
        PART[RP, on=.(quadrant), q_new := {
          lox = x.x > i.x
          1L + lox
        }]
        PART[max.obs.rows, quadrant := paste0(quadrant, q_new)]

        new.parts=length(unique(PART$quadrant))
      }

      if(noise.reduction=='median'){
        RP= PART[max.obs.rows, lapply(.SD, median), by=quadrant, .SDcols = x:y]
        RP[, prior.quadrant := (quadrant)]

        PART[max.obs.rows , prior.quadrant := (quadrant)]
        old.parts=length(unique(PART$quadrant))
        PART[RP, on=.(quadrant), q_new := {
          lox = x.x > i.x
          1L + lox
        }]
        PART[max.obs.rows, quadrant := paste0(quadrant, q_new)]

        new.parts=length(unique(PART$quadrant))
      }

      if(old.parts==new.parts) break
      i=i+1L
    }


    if(Voronoi==T){
      abline(v=RP$x,lty=3)
      points(RP$x,RP$y,pch=15,lwd=2,col='red')
      title(main=paste0("NNS Order = ",i),cex.main=2)
    }

    setnames(RP,"prior.quadrant","quadrant")
    PART[,`:=`(counts = NULL, q_new = NULL)]
    DT=setorder(PART[], quadrant,x, y)[]
    RP=setorder(RP[],quadrant)[]
    return(list("order"=i,"dt"=DT,"regression.points"=RP))
  }

}
