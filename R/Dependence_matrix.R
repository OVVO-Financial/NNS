NNS.dep.matrix <- function(x, order = NULL, degree = NULL, asym = FALSE){

  n <- ncol(x)
  if(is.null(n)){
      stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }

  if(any(class(x)=="tbl")) x <- as.data.frame(x)

  x <- data.matrix(x)

  if(dim(x)[1] < 20 ) {
    order <- 2
    asym <- TRUE
  }

  raw.rhos_lower <- list()
  raw.deps_lower <- list()
  raw.both_lower <- list()

  for(i in 1 : (n-1)){
        raw.both_lower[[i]] <- sapply((i + 1) : n, function(b) NNS.dep(x[ , i], x[ , b], print.map = FALSE, asym = asym))

        raw.rhos_lower[[i]] <- unlist(raw.both_lower[[i]][row.names(raw.both_lower[[i]])=="Correlation"])
        raw.deps_lower[[i]] <- unlist(raw.both_lower[[i]][row.names(raw.both_lower[[i]])=="Dependence"])
  }


  rhos <- matrix(, n, n)
  rhos[lower.tri(rhos, diag = FALSE)] <- unlist(raw.rhos_lower)

  deps <- matrix(0, n, n)
  deps[lower.tri(deps, diag = FALSE)] <- unlist(raw.deps_lower)

    if(!asym){
        rhos <- pmax(rhos, t(rhos), na.rm = TRUE)
        deps <- pmax(deps, t(deps), na.rm = TRUE)
    } else {
        rhos_upper <- matrix(0, n, n)
        deps_upper <- matrix(0, n, n)

        raw.rhos_upper <- list()
        raw.deps_upper <- list()
        raw.both_upper <- list()

        for(i in 1 : (n-1)){
            raw.both_upper[[i]] <- sapply((i + 1) : n, function(b) NNS.dep(x[ , b], x[ , i], print.map = FALSE, asym = asym))

            raw.rhos_upper[[i]] <- unlist(raw.both_upper[[i]][row.names(raw.both_upper[[i]])=="Correlation"])
            raw.deps_upper[[i]] <- unlist(raw.both_upper[[i]][row.names(raw.both_upper[[i]])=="Dependence"])
        }

        rhos[is.na(rhos)] <- 0
        deps[is.na(deps)] <- 0

        rhos_upper[lower.tri(rhos_upper, diag=FALSE)] <- unlist(raw.rhos_upper)
        rhos_upper <- t(rhos_upper)

        deps_upper[lower.tri(deps_upper, diag=FALSE)] <- unlist(raw.deps_upper)
        deps_upper <- t(deps_upper)

    rhos <- rhos + rhos_upper
    deps <- deps + deps_upper
}

diag(rhos) <- 1
diag(deps) <- 1

colnames(rhos) <- colnames(x)
colnames(deps) <- colnames(x)
rownames(rhos) <- colnames(x)
rownames(deps) <- colnames(x)

return(list("Correlation" = rhos,
            "Dependence" = deps))

}



