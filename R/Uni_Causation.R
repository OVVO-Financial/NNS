Uni.caus <- function(x, y, tau, plot = TRUE){

  if(tau=="cs") tau <- 0
  if(tau=="ts") tau <- 3

  xy <- NNS.norm(cbind(x, y), linear = FALSE, chart.type = NULL)

  NNS.x <- unlist(xy[ , 1])
  NNS.y <- unlist(xy[ , 2])

  min.length <- min(length(x), length(y))

  x.vectors <- list(tau+1)
  y.vectors <- list(tau+1)

  ## Create tau vectors
  if(tau > 0){
    for (i in 0:tau){
        x.vectors[[paste('x.tau.', i, sep = "")]] <- numeric(0L)
        y.vectors[[paste('y.tau.', i, sep = "")]] <- numeric(0L)
        start <- tau - i + 1
        end <- min.length - i
        x.vectors[[i + 1]] <- x[start : end]
        y.vectors[[i + 1]] <- y[start : end]
    }

      x.vectors.tau <- do.call(cbind, x.vectors)
      y.vectors.tau <- do.call(cbind, y.vectors)

      ## Normalize x to x.tau
      x.norm.tau <- unlist(NNS.norm(x.vectors.tau)[ , 1])

      ## Normalize y to y.tau
      y.norm.tau <- unlist(NNS.norm(y.vectors.tau)[ , 1])

  } else {
      x.norm.tau <- x
      y.norm.tau <- y
  }



  ## Normalize x.norm.tau to y.norm.tau
  x.tau.y.tau <- NNS.norm(cbind(x.norm.tau, y.norm.tau))
  x.norm.to.y <- as.vector(unlist(x.tau.y.tau[ , 1]))
  y.norm.to.x <- as.vector(unlist(x.tau.y.tau[ , 2]))


  ## Conditional Probability from Normalized Variables P(x.norm.to.y | y.norm.to.x)
  P.x.given.y <- 1 - (LPM.ratio(1, min(y.norm.to.x), x.norm.to.y) + UPM.ratio(1, max(y.norm.to.x), x.norm.to.y))


  ## Correlation of Normalized Variables
  rho.x.y <- NNS.dep(y.norm.to.x, x.norm.to.y, asym = TRUE)$Dependence

  Causation.x.given.y <- P.x.given.y * rho.x.y


  if(plot){
      original.par <- par(no.readonly = TRUE)
      par(mfrow = c(3, 1))

      ## Raw Variable Plot
      ymin <- min(c(min(x), min(y)))
      ymax <- max(c(max(x), max(y)))
      par(mar = c(2, 4, 0, 1))
      plot(y,type = 'l', ylim = c(ymin, ymax), ylab = 'STANDARDIZED', col = 'red', lwd = 3)
      lines(x, col = 'steelblue',lwd = 3)
      legend('top', c("X", "Y"), lty = 1,lwd = c(3, 3),
           col = c('steelblue', 'red'), ncol = 2)

      ## Time Normalized Variables Plot
      ymin <- min(c(min(x.norm.tau), min(y.norm.tau)))
      ymax <- max(c(max(x.norm.tau), max(y.norm.tau)))
      par(mar = c(2, 4, 0, 1))
      plot(y.norm.tau, type = 'l', ylim = c(ymin, ymax), ylab = 'TIME NORMALIZED', col = 'red', lwd = 3)
      lines(x.norm.tau, col = 'steelblue', lwd = 3)
      legend('top', c("X", "Y"), lty = 1, lwd = c(3, 3),
           col = c('steelblue', 'red'), ncol = 2)

      ## Time Normalized Variables Normalized to each other Plot
      ymin <- min(c(min(x.norm.to.y), min(y.norm.to.x)))
      ymax <- max(c(max(x.norm.to.y), max(y.norm.to.x)))
      par(mar = c(2, 4, 0, 1))
      plot(y.norm.to.x, type = 'l', ylim = c(ymin, ymax), ylab = 'X & Y NORMALIZED', col='red', lwd = 3)
      lines(x.norm.to.y, col = 'steelblue', lwd = 3)
      legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
           col = c('steelblue', 'red'), ncol = 2)

      par(original.par)
  }

  return(Causation.x.given.y)

}
