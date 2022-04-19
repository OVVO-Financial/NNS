NNS.dep.matrix <- function(x, order = NULL, degree = NULL, asym = FALSE, ncores = NULL){

  n <- ncol(x)
  if(is.null(n)){
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }

  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)

  x <- data.matrix(x)

  if(dim(x)[1] < 20 ) order <- 2

  upper_lower <- function(x, y, asym){
    basic_dep <- NNS.dep(x, y, print.map = FALSE, asym = asym, ncores = 1)
    if(asym){
      asym_dep <- NNS.dep(y, x, print.map = FALSE, asym = asym, ncores = 1)
      return(list("Upper_cor" = basic_dep$Correlation,
                  "Upper_dep" = basic_dep$Dependence,
                  "Lower_cor" = asym_dep$Correlation,
                  "Lower_dep" = asym_dep$Dependence))
    } else {
      return(list("Upper_cor" = basic_dep$Correlation,
                  "Upper_dep" = basic_dep$Dependence,
                  "Lower_cor" = basic_dep$Correlation,
                  "Lower_dep" = basic_dep$Dependence))
    }
  }

  raw.both <- list((n-1))

  if(is.null(ncores)) {
    num_cores <- as.integer(parallel::detectCores()) - 1
  } else {
    num_cores <- ncores
  }

  if(num_cores>1){
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
  }


  i <- 1:(n-1)

  raw.both <- foreach(i = 1 : (n-1), .packages = c("NNS", "data.table"))%dopar%{
    raw.both[[i]] <-  sapply((i + 1) : n, function(b) upper_lower(x[ , i], x[ , b], asym = asym))
  }

  if(num_cores>1){
    parallel::stopCluster(cl)
    registerDoSEQ()
  }

  raw.both <- unlist(raw.both)
  l <- length(raw.both)

  raw.rhos_upper <- raw.both[seq(1, l, 4)]
  raw.deps_upper <- raw.both[seq(2, l, 4)]
  raw.rhos_lower <- raw.both[seq(3, l, 4)]
  raw.deps_lower <- raw.both[seq(4, l, 4)]

  rhos <- matrix(0, n, n)
  deps <- matrix(0, n, n)

  if(!asym){
    rhos[lower.tri(rhos, diag = FALSE)] <- (unlist(raw.rhos_upper) + unlist(raw.rhos_lower)) / 2
    deps[lower.tri(deps, diag = FALSE)] <- (unlist(raw.deps_upper) + unlist(raw.deps_lower)) / 2

    rhos <- pmax(rhos, t(rhos), na.rm = TRUE)
    deps <- pmax(deps, t(deps), na.rm = TRUE)
  } else {
    rhos[lower.tri(rhos, diag = FALSE)] <- unlist(raw.rhos_lower)
    deps[lower.tri(deps, diag = FALSE)] <- unlist(raw.deps_lower)

    rhos_upper <- matrix(0, n, n)
    deps_upper <- matrix(0, n, n)

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
