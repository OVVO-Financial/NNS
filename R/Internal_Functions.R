### Continuous Mode of a distribution
mode <- function(x) NNS.mode(x, discrete = FALSE, multi = FALSE)


### Classification Mode of a distribution
mode_class <- function(x) NNS.mode(x, discrete = TRUE, multi = FALSE)


### Gravity of a distribution
gravity <- function(x, discrete = FALSE){
  l <- length(x)
  if(l <= 3) return(median(x))
  if(length(unique(x))==1) return(x[1])
  x_s <- x[order(x)]
  range <- abs(x_s[l]-x_s[1])

  if(range == 0) return(x[1])
  
  q1 <- sum(x_s[floor(l*.25)]+((l*.25)%%1 * (x_s[ceiling(l*.25)] - x_s[floor(l*.25)])))
  q2 <- (x_s[floor(l*.5)]+x_s[ceiling(l*.5)])/2
  q3 <- sum(x_s[floor(l*.75)]+((l*.75)%%1 * (x_s[ceiling(l*.75)] - x_s[floor(l*.75)])))

  z <- MESS::bin(x_s, range/128, origin = x_s[1], missinglast = FALSE)
  lz <- length(z$counts)
  max_z <- z$counts==max(z$counts)
  
  if(sum(max_z)>1){
    z_ind <- 1:lz
  } else {
    z_c <- which.max(z$counts)
    z_ind <- max(1, (z_c - 1)):min(lz,(z_c + 1))
  }

  z_names <- seq(x_s[1], x_s[l], z$width)

  m <- sum(z_names[z_ind] * z$counts[z_ind] )/sum(z$counts[z_ind])
  mu <- sum(x)/l
  
  res <- (q2 + m + mu + mean(c(q1, q2, q3)))/4
  if(is.na(res)) final <- q2 else final <- res
  if(discrete) return(ifelse(final%%1 < .5, floor(final), ceiling(final))) else return(final)
} 


gravity_class <- function(x) gravity(x, discrete = TRUE)


### Factor to dummy variable
factor_2_dummy <- function(x){
  if(is.factor(x) && length(unique(x)) > 1){
    x <- unlist(x)
    output <- model.matrix(~(x) -1, x)[,-1]
  } else {
    x <- unlist(x)
    output <- as.numeric(x)
  }
  output
}

### Factor to dummy variable FULL RANK
factor_2_dummy_FR <- function(x){
  if(is.factor(x) && length(unique(x)) > 1){
    x <- unlist(x)
    output <- model.matrix(~(x) -1, x)
  } else {
    x <- unlist(x)
    output <- as.numeric(x)
  }
  output
}



### Generator for 1:length(lag) vectors in NNS.ARMA
generate.vectors <- function(x, l){
  Component.series <- list()
  Component.index <- list()

  for (i in 1:length(l)){
    CS <- rev(x[seq(length(x)+1, 1, -l[i])])
    CS <- CS[!is.na(CS)]
    Component.series[[paste('Series.', i, sep = "")]] <- CS
    Component.index[[paste('Index.', i, sep = "")]] <- (1 : length(CS))
  }
  return(list(Component.index = Component.index, Component.series = Component.series))
}


### Weight and lag function for seasonality in NNS.ARMA
ARMA.seas.weighting <- function(sf,mat){
  M <- mat
  n <- ncol(M)
  if(is.null(n)){
    return(list(lag = M[1], Weights = 1))
  }

  if(n == 1){
    return(list(lag = 1, Weights = 1))
  }

  if(n > 1){
    if(sf){
      lag <- M$all.periods$Period[1]
      Weights <- 1
      return(list(lag = lag, Weights = Weights))
    }

    # Determine lag from seasonality test
    if(!sf){
      lag <- na.omit(unlist(M$Period))
      Observation.weighting <- (1 / sqrt(lag))
      if(is.na(M$Coefficient.of.Variation)  && length(M$Coefficient.of.Variation)==1){
        Lag.weighting <- 1
      } else {
        Lag.weighting <- (unlist(M$Variable.Coefficient.of.Variation) - unlist(M$Coefficient.of.Variation))
      }
      Weights <- (Lag.weighting * Observation.weighting) / sum(Lag.weighting * Observation.weighting)
      return(list(lag = lag, Weights = Weights))
    }
  }
}


### Lag matrix generator for NNS.VAR
### Vector of tau for single different tau per variables tau = c(1, 4)
### List of tau vectors for multiple different tau per variables tau = list(c(1,2,3), c(4,5,6))
lag.mtx <- function(x, tau){
  colheads <- NULL

  max_tau <- max(unlist(tau))

  if(is.null(dim(x)[2])) {
    colheads <- noquote(as.character(deparse(substitute(x))))
    x <- t(t(x))
  }

  j.vectors <- list()

  for(j in 1:ncol(x)){
    if(is.null(colheads)){
      colheads <- colnames(x)[j]

      colheads <- noquote(as.character(deparse(substitute(colheads))))
    }

    x.vectors <- list()
    heads <- paste0(colheads, "_tau_")
    heads <- gsub('"', '' ,heads)

    for (i in 0:max_tau){
      x.vectors[[paste(heads, i, sep = "")]] <- numeric(0L)
      start <- max_tau - i + 1
      end <- length(x[,j]) - i
      x.vectors[[i + 1]] <- x[start : end, j]
    }

    j.vectors[[j]] <- do.call(cbind, x.vectors)
    colheads <- NULL
  }
  mtx <- as.data.frame(do.call(cbind, j.vectors))

  relevant_lags <- list(length(tau))
  if(length(unlist(tau)) > 1){
    for(i in 1:(length(tau))){
        relevant_lags[[i]] <- c((i-1)*max_tau + i, (i-1)*max_tau + unlist(tau[[i]]) + i)
    }

    relevant_lags <- sort(unlist(relevant_lags))
    mtx <- mtx[ , relevant_lags]
  }
  vars <- which(grepl("tau_0", colnames(mtx)))

  everything_else <- seq_len(dim(mtx)[2])[-vars]
  mtx <- mtx[,c(vars, everything_else)]

  return(mtx)
}




### Refactored meboot::meboot.part function using tdigest
NNS.meboot.part <- function(x, n, z, xmin, xmax, desintxb, reachbnd)
{
  # Generate random numbers from the [0,1] uniform interval
  p <- runif(n, min=0, max=1)

  q <- quantile(x, p)

  ref1 <- which(p <= (1/n))
  if(length(ref1) > 0){
    qq <- approx(c(0,1/n), c(xmin,z[1]), p[ref1])$y
    q[ref1] <- qq
    if(!reachbnd)  q[ref1] <- qq + desintxb[1]-0.5*(z[1]+xmin)
  }

  ref4 <- which(p == ((n-1)/n))
  if(length(ref4) > 0)
    q[ref4] <- z[n-1]

  ref5 <- which(p > ((n-1)/n))
  if(length(ref5) > 0){
    # Right tail proportion p[i]
    qq <- approx(c((n-1)/n,1), c(z[n-1],xmax), p[ref5])$y
    q[ref5] <- qq   # this implicitly shifts xmax for algorithm
    if(!reachbnd)  q[ref5] <- qq + desintxb[n]-0.5*(z[n-1]+xmax)
    # such that the algorithm gives xmax when p[i]=1
    # this is the meaning of reaching the bounds xmax and xmin
  }

  q

}

### Refactored meboot::expand.sd function
NNS.meboot.expand.sd <- function(x, ensemble, fiv=5){
  sdx <- if (is.null(ncol(x))) sd(x) else apply(x, 2, sd)

  sdf <- c(sdx, apply(ensemble, 2, sd))

  sdfa <- sdf/sdf[1]  # ratio of actual sd to that of original data
  sdfd <- sdf[1]/sdf  # ratio of desired sd to actual sd

  # expansion is needed since some of these are <1 due to attenuation
  mx <- 1+(fiv/100)
  # following are expansion factors
  id <- which(sdfa < 1)
  if (length(id) > 0) sdfa[id] <- runif(n=length(id), min=1, max=mx)

  sdfdXsdfa <- sdfd[-1]*sdfa[-1]

  id <- which(floor(sdfdXsdfa) > 0)

  if (length(id) > 0) {
    if(length(id) > 1) ensemble[,id] <- ensemble[,id] %*% diag(sdfdXsdfa[id]) else ensemble[,id] <- ensemble[,id] * sdfdXsdfa[id]
  }

  if(is.ts(x)) ensemble <- ts(ensemble, frequency=frequency(x), start=start(x))


  ensemble
}


is.discrete <- function(x) is.factor(x) || is.character(x) || is.logical(x)

