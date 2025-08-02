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
  dep.mtx <- NNS.dep(cbind(y.norm.to.x, x.norm.to.y), asym = TRUE)$Dependence
  rho.x.y <- dep.mtx[1, 2]
  rho.y.x <- dep.mtx[2, 1]

  Causation.x.given.y <- mean(c(P.x.given.y * rho.x.y, max(0, (rho.x.y - rho.y.x))))


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



# Internal helper: extract canonical signed causation from the third element of cp.
signed_from_cp <- function(cp){
  # Extract the third element (log-ratio in direction of stronger causation) and cap infinities
  raw_val <- as.numeric(cp[3])
  if(is.infinite(raw_val)) raw_val <- sign(raw_val) * 100
  return(raw_val)
}

# helper: cap infinite values at 100 (preserves sign) and apply tanh normalization; vectorized and scalar-safe.
cap_inf100_scalar <- function(v, cap = 100, scaling_factor = 2){  # switched to scaling_factor=2 to match log-ratio normalization
  v2 <- v
  v2[is.infinite(v2)] <- sign(v2[is.infinite(v2)]) * cap
  v2 <- ifelse(abs(v2) > cap, sign(v2) * cap, v2)
  return(v2)
#  return(tanh(v2 / scaling_factor))
}




# log-ratio logic now consolidated in core; helper removed.
# Core causation computation (no permutation logic) extracted from original implementation.
NNS.caus_core <- function(x, y = NULL,
                          factor.2.dummy = FALSE,
                          tau = 0,
                          plot = FALSE,
                          p.value = FALSE,
                          nperm = 100L,
                          permute = c("y","x","both"),
                          seed = NULL,
                          conf.int = 0.95){
  if(!is.null(y))  if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")
  if(is.null(y))  if(sum(is.na(x)) > 0) stop("You have some missing values, please address.")
  
  orig.tau <- tau
  orig.plot <- plot
  
  if(any(class(x)%in%c("tbl","data.table")) && dim(x)[2]==1) x <- as.vector(unlist(x))
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
  if(!is.null(y) && any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))
  
  if(factor.2.dummy){
    if(!is.null(dim(x))){
      if(!is.numeric(x)){
        x <- do.call(cbind, lapply(x, factor_2_dummy_FR))
      } else {
        x <- apply(x, 2, as.double)
      }
      if(is.list(x)){
        x <- do.call(cbind, x)
        x <- apply(x, 2, as.double)
      }
    } else {
      x <- factor_2_dummy(x)
      if(is.null(dim(x))){
        x <- as.double(x)
      } else {
        x <- apply(x, 2, as.double)
      }
    }
  }
  
  if(!is.null(y)){
    if(is.factor(y)) y <- as.numeric(y)
    
    if(is.numeric(tau)){
      Causation.x.given.y <- Uni.caus(x, y, tau = tau, plot = FALSE)
      Causation.y.given.x <- Uni.caus(y, x, tau = tau, plot = FALSE)
      Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
      Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
      if(Causation.x.given.y == Causation.y.given.x ||
         Causation.x.given.y == 0 || Causation.y.given.x == 0){
        Causation.x.given.y <- Uni.caus(x, y, tau = tau, plot = FALSE)
        Causation.y.given.x <- Uni.caus(y, x, tau = tau, plot = FALSE)
        Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
        Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
      }
    }
    
    if(identical(tau, "cs")){
      Causation.x.given.y <- Uni.caus(x, y, tau = 0, plot = FALSE)
      Causation.y.given.x <- Uni.caus(y, x, tau = 0, plot = FALSE)
      Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
      Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
      if(Causation.x.given.y == Causation.y.given.x ||
         Causation.x.given.y == 0 || Causation.y.given.x == 0){
        Causation.x.given.y <- Uni.caus(x, y, tau = 0, plot = FALSE)
        Causation.y.given.x <- Uni.caus(y, x, tau = 0, plot = FALSE)
        Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
        Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
      }
    }
    
    if(identical(tau, "ts")){
      l <- length(x)
      x_tau <- NNS.seas(x, plot = FALSE)$periods
      y_tau <- NNS.seas(y, plot = FALSE)$periods
      
      x_tau <- x_tau[x_tau <= (l)^(1/2)][1]
      y_tau <- y_tau[y_tau <= (l)^(1/2)][1]
      
      Causation.y.given.x <- Uni.caus(y, x, tau = x_tau, plot = FALSE)
      Causation.x.given.y <- Uni.caus(x, y, tau = y_tau, plot = FALSE)
      Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
      Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
    }
    
    # Choose plotting direction for diagnostics (basing on which direction is stronger) and correct direction naming according to updated semantics:
    # If Causation.y.given.x >= Causation.x.given.y then x causes y is stronger (so label C(x--->y)).
    if(abs(Causation.y.given.x) >= abs(Causation.x.given.y)){
      if(plot){
        if(identical(tau, "cs")) tau_plot <- 0 else if(identical(tau, "ts")) tau_plot <- mean(c(x_tau, y_tau)) else tau_plot <- tau
        Uni.caus(y, x, tau = tau_plot, plot = plot)
      }
      eps <- .Machine$double.eps
      net_log_ratio <- sign(Causation.y.given.x) * log((abs(Causation.y.given.x) + eps)/(abs(Causation.x.given.y) + eps))
      cp <- c(Causation.x.given.y = Causation.x.given.y,
              Causation.y.given.x = Causation.y.given.x,
              "C(x--->y)" = net_log_ratio)
    } else {
      if(plot){
        if(identical(tau, "cs")) tau_plot <- 0 else if(identical(tau, "ts")) tau_plot <- mean(c(x_tau, y_tau)) else tau_plot <- tau
        Uni.caus(x, y, tau = tau_plot, plot = plot)
      }
      eps <- .Machine$double.eps
      net_log_ratio <- sign(Causation.x.given.y) * log((abs(Causation.x.given.y) + eps)/(abs(Causation.y.given.x) + eps))
      cp <- c(Causation.x.given.y = Causation.x.given.y,
              Causation.y.given.x = Causation.y.given.x,
              "C(y--->x)" = net_log_ratio)
    }
    cp[3] <- cap_inf100_scalar(cp[3])
    return(cp)
  } else {
    return(NNS.caus.matrix(x, tau = orig.tau, factor.2.dummy = factor.2.dummy, plot = orig.plot, p.value = p.value, nperm = nperm, conf.int = conf.int, seed = seed))  
  }
}