### Continuous Mode of a distribution
mode <- function(x) {
  if(length(x) > 1){
      d <- density(x)
      d$x[which.max(d$y)]
  } else {
    x
  }
}

### Classification Mode of a distribution
mode_class <- function(x){
  x <- na.omit(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

### Factor to dummy variable
factor_2_dummy <- function(x){
  if(class(x) == "factor"){
    output <- model.matrix(~(x) -1, x)[,-1]
  } else {
    output <- x
  }
  output
}

### Factor to dummy variable FULL RANK
factor_2_dummy_FR <- function(x){
  if(class(x) == "factor"){
    output <- model.matrix(~(x) -1, x)
  } else {
    output <- x
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
      if(is.na(M$Coefficient.of.Variance)  && length(M$Coefficient.of.Variance)==1){
        Lag.weighting <- 1
      } else {
        Lag.weighting <- (M$Variable.Coefficient.of.Variance - M$Coefficient.of.Variance)
      }
      Weights <- (Lag.weighting * Observation.weighting) / sum(Lag.weighting * Observation.weighting)
      return(list(lag = lag, Weights = Weights))
    }
  }
}


### Lag matrix generator for NNS.VAR
lag.mtx <- function(x, tau){
  colheads <- NULL

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
    heads <- paste0(colheads, ".tau.")
    heads <- gsub('"', '' ,heads)

    for (i in 0:tau){
      x.vectors[[paste(heads, i, sep = "")]] <- numeric(0L)
      start <- tau - i + 1
      end <- length(x[,j]) - i
      x.vectors[[i + 1]] <- x[start : end, j]
    }

    j.vectors[[j]] <- do.call(cbind, x.vectors)
    colheads <- NULL
  }

  return(as.data.frame(do.call(cbind, j.vectors)))
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



