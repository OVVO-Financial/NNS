Uni.caus <- function(x,y,tau,plot=TRUE,scale=FALSE,time.series=FALSE){

  xy=NNS.norm(cbind(x,y),linear = TRUE)
  NNS.x=xy[,1];NNS.y=xy[,2]

  if(mean(NNS.x)!=mean(NNS.y)){
    xy=apply(cbind(x,y),2,function(b) (b-min(b))/(max(b)-min(b)))
    x=xy[,1];y=xy[,2]
    if( time.series | scale |
        LPM.ratio(2,mean(x),x)<.25 | LPM.ratio(2,mean(y),y)<.25 |
        LPM.ratio(2,mean(x),x)>.75 | LPM.ratio(2,mean(y),y)>.75
    )
    {
      x=scale(x)[,1];y=scale(y)[,1]
    }
  }

  if(mean(NNS.x)==mean(NNS.y)){
    x=NNS.x;y=NNS.y
  if( time.series | scale |
      LPM.ratio(2,mean(x),x)<.25 | LPM.ratio(2,mean(y),y)<.25 |
      LPM.ratio(2,mean(x),x)>.75 | LPM.ratio(2,mean(y),y)>.75
      )
    {
    x=scale(x)[,1];y=scale(y)[,1]
  }
}

  min.length = min(length(x),length(y))

  x.vectors = list()
  y.vectors = list()

  ## Create tau vectors
  if(tau>0){
    for (i in 0:tau){
      x.vectors[[paste('x.tau.',i,sep="")]] <- numeric(0L)
      y.vectors[[paste('y.tau.',i,sep="")]] <- numeric(0L)
      start = tau-i+1
      end = min.length-i
      x.vectors[[i+1]] = x[start:end]
      y.vectors[[i+1]] = y[start:end]
    }

    x.vectors.tau = do.call(cbind,x.vectors)
    y.vectors.tau = do.call(cbind,y.vectors)


    ## Normalize x to x.tau
    x.norm.tau <- rowMeans(NNS.norm(x.vectors.tau))

    ## Normalize y to y.tau
    y.norm.tau <- rowMeans(NNS.norm(y.vectors.tau))
  } else
  {
    x.norm.tau <- x
    y.norm.tau <- y
  }


  ## Normalize x.norm.tau to y.norm.tau
  x.tau.y.tau = NNS.norm(cbind(x.norm.tau,y.norm.tau))
  x.norm.to.y = x.tau.y.tau[,1]
  y.norm.to.x = x.tau.y.tau[,2]


  ## Conditional Probability from Normalized Variables P(x.norm.to.y | y.norm.to.x)
  P.x.given.y = UPM.ratio(1,min(x.norm.to.y),y.norm.to.x) - UPM.ratio(1,max(x.norm.to.y),y.norm.to.x)


  ## Correlation of Normalized Variables
  rho.x.y = NNS.dep(x.norm.to.y,y.norm.to.x)$Dependence

  Causation.x.given.y= P.x.given.y*rho.x.y


  if(plot){
    par(mfrow=c(3,1))

    ## Raw Variable Plot
    ymin = min(c(min(x),min(y)))
    ymax = max(c(max(x),max(y)))
    par(mar=c(2, 4, 0, 1))
    plot(y,type = 'l', ylim=c(ymin, ymax),ylab='STANDARDIZED',col='red',lwd = 3)
    lines(x, col = 'steelblue',lwd = 3)
    legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
           col=c('steelblue','red'),ncol=2)



    ## Time Normalized Variables Plot
    ymin = min(c(min(x.norm.tau),min(y.norm.tau)))
    ymax = max(c(max(x.norm.tau),max(y.norm.tau)))
    par(mar=c(2, 4, 0, 1))
    plot(y.norm.tau,type = 'l', ylim=c(ymin, ymax),ylab='TIME NORMALIZED',col='red',lwd = 3)
    lines(x.norm.tau, col='steelblue',lwd=3)
    legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
           col=c('steelblue','red'),ncol=2)

    ## Time Normalized Variables Normalized to each other Plot
    ymin = min(c(min(x.norm.to.y),min(y.norm.to.x)))
    ymax = max(c(max(x.norm.to.y),max(y.norm.to.x)))
    par(mar=c(2, 4, 0, 1))
    plot(y.norm.to.x,type = 'l', ylim=c(ymin, ymax),ylab='X & Y NORMALIZED',col='red',lwd = 3)
    lines(x.norm.to.y,col='steelblue',lwd=3)
    legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
           col=c('steelblue','red'),ncol=2)

    par(mfrow=c(1,1))
  }

  return(Causation.x.given.y)

}
