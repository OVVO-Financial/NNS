#' NNS Causation uni-directional  (INTERNAL CALL FOR \link{NNS.caus})
#'
#' Returns the uni-directional causality from observational data between two variables.  Causation nets out the univariate effect.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param tau integer; Number of lagged observations to consider
#' @param plot logical; \code{TRUE} (default) Plots the raw variables, tau normalized, and cross-normalized variables.
#' @keywords causation
#' @author Fred Viole, OVVO Financial Systems

Uni.caus <- function(x,y,tau,plot=TRUE){


  min.length = min(length(x),length(y))

  x.vectors = list()
  y.vectors = list()

## Create tau vectors
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
  x.norm.tau <- NNS.norm(x.vectors.tau)[,1]

## Normalize y to y.tau
  y.norm.tau <- NNS.norm(y.vectors.tau)[,1]


  ## Normalize x.norm.tau to y.norm.tau
  x.tau.y.tau = NNS.norm(cbind(x.norm.tau,y.norm.tau))
  x.norm.to.y = x.tau.y.tau[,1]
  y.norm.to.x = x.tau.y.tau[,2]


  ## Conditional Probability from Normalized Variables P(x.norm.to.y | y.norm.to.x)
  P.x.given.y = UPM(0,min(x.norm.to.y),y.norm.to.x) - UPM(0,max(x.norm.to.y),y.norm.to.x)



  ## Correlation of Normalized Variables
  rho.x.y = NNS.dep(x.norm.to.y,y.norm.to.x)$Dependence

  Causation.x.given.y= P.x.given.y*rho.x.y


  if(plot==TRUE){
  par(mfrow=c(3,1))

  ## Raw Variable Plot
  ymin = min(c(min(x),min(y)))
  ymax = max(c(max(x),max(y)))
  par(mar=c(2, 4, 0, 1))
  plot(y,type = 'l', ylim=c(ymin, ymax),ylab='RAW',col='red',lwd = 3)
  lines(x, col = 'steelblue',lwd = 3)
  legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
         col=c('steelblue','red'),ncol=2,bty ="n")



  ## Time Normalized Variables Plot
  ymin = min(c(min(x.norm.tau),min(y.norm.tau)))
  ymax = max(c(max(x.norm.tau),max(y.norm.tau)))
  par(mar=c(2, 4, 0, 1))
  plot(y.norm.tau,type = 'l', ylim=c(ymin, ymax),ylab='TIME NORMALIZED',col='red',lwd = 3)
  lines(x.norm.tau, col='steelblue',lwd=3)
  legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
         col=c('steelblue','red'),ncol=2,bty ="n")

  ## Time Normalized Variables Normalized to each other Plot
  ymin = min(c(min(x.norm.to.y),min(y.norm.to.x)))
  ymax = max(c(max(x.norm.to.y),max(y.norm.to.x)))
  par(mar=c(2, 4, 0, 1))
  plot(y.norm.to.x,type = 'l', ylim=c(ymin, ymax),ylab='X & Y NORMALIZED',col='red',lwd = 3)
  lines(x.norm.to.y,col='steelblue',lwd=3)
  legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
         col=c('steelblue','red'),ncol=2,bty ="n")

  par(mfrow=c(1,1))
  }

  return(Causation.x.given.y)

}
