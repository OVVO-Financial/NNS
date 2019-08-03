### Mode of a distribution
mode <- function(x){
  if(length(na.omit(x)) > 1){
      d <- density(na.omit(x))
      d$x[which.max(d$y)]
  } else {
      x
  }
}

### Factor to dummy variable
factor_2_dummy <- function(x){
  if(class(x) == "factor"){
      output <- model.matrix(~x -1, x)[,-1]
  } else {
      output <- x
  }
  output
}

### Factor to dummy variable FULL RANK
factor_2_dummy_FR <- function(x){
  if(class(x) == "factor"){
    output <- model.matrix(~x -1, x)
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
            lag <- M$Period
            Observation.weighting <- (1 / sqrt(lag))
            Lag.weighting <- (M$Variable.Coefficient.of.Variance - M$Coefficient.of.Variance)
            Weights <- (Lag.weighting * Observation.weighting) / sum(Lag.weighting * Observation.weighting)
            return(list(lag = lag, Weights = Weights))
        }
      }
    }
