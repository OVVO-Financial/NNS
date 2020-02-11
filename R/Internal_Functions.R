### Continuous Mode of a distribution
mode <- function(x) {
      d <-tryCatch(density(na.omit(as.numeric(x))), error = function(e) { median(x)})
      tryCatch(d$x[which.max(d$y)], error = function(e) {d})
  }

### Classification Mode of a distribution
mode_class <- function(x){
  x <- na.omit(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

### Central Tendency
gravity <- function(x){
  (mean(x) + median(x) + mode(x) + mean(range(x)) ) / 4
}


### cbind different length vectors
alt_cbind <- function(x,y,first = FALSE) {
  if(length(x)<length(y)) {
    if(first) x = c(rep(NA, length(y)-length(x)),x);y=y
    if(!first) x = c(x,rep(NA, length(y)-length(x)));y=y
  }
  if(length(y)<length(x)) {
    if(first) y = c(rep(NA, length(x)-length(y)),y);x=x
    if(!first) y = c(y,rep(NA, length(x)-length(y)));x=x
  }

  return(cbind(x,y))

}


### Factor to dummy variable
factor_2_dummy <- function(x){
  if(class(x) == "factor" & length(unique(x)) > 1){
    output <- model.matrix(~(x) -1, x)[,-1]
  } else {
    output <- as.numeric(x)
  }
  output
}

### Factor to dummy variable FULL RANK
factor_2_dummy_FR <- function(x){
  if(class(x) == "factor" & length(unique(x)) > 1){
    output <- model.matrix(~(x) -1, x)
  } else {
    output <- as.numeric(x)
  }
  output
}



### Generator for 1:length(lag) vectors in NNS.ARMA
generate.vectors <- function(x,l){
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
      lag <- na.omit(M$Period)
      Observation.weighting <- (1 / sqrt(lag))
      if(is.na(M$Coefficient.of.Variation)  && length(M$Coefficient.of.Variation)==1){
        Lag.weighting <- 1
      } else {
        Lag.weighting <- (M$Variable.Coefficient.of.Variation - M$Coefficient.of.Variation)
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

  relevant_lags <- list()
  if(length(unlist(tau)) > 1){
    for(i in 1:(length(tau))){
        relevant_lags[[i]] <- c((i-1)*max_tau + i, (i-1)*max_tau + unlist(tau[[i]]) + i)
    }

    relevant_lags <- sort(unlist(relevant_lags))
    mtx <- mtx[ , relevant_lags]
  }

  return(mtx)
}



### Row products for NNS.dep.hd
RP <- function(x, rows = NULL, cols = NULL, na.rm = FALSE) {

  if (!is.null(rows) && !is.null(cols)) x <- x[rows, cols, drop = FALSE]
  else if (!is.null(rows)) x <- x[rows, , drop = FALSE]
  else if (!is.null(cols)) x <- x[, cols, drop = FALSE]

  n <- nrow(x)
  y <- double(length = n)

  if (n == 0L) return(y)

  for (ii in seq_len(n)) {
    y[ii] <- prod(x[ii, , drop = TRUE], na.rm = na.rm)
  }

  y
}



