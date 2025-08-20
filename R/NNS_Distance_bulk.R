NNS.distance.bulk <- function(rpm, Xtest, k, class = NULL) {
  rpm <- data.table::as.data.table(rpm)
  stopifnot("y.hat" %in% names(rpm))
  
  # drop y.hat, align columns with Xtest by name if possible
  Xrpm <- as.data.frame(rpm[, !"y.hat"])
  if (!is.null(colnames(Xrpm)) && !is.null(colnames(Xtest))) {
    cmn <- intersect(colnames(Xrpm), colnames(Xtest))
    Xrpm  <- as.matrix(Xrpm[, cmn, drop = FALSE])
    Xtest <- as.matrix(as.data.frame(Xtest)[, cmn, drop = FALSE])
  } else {
    Xrpm  <- as.matrix(Xrpm)
    Xtest <- as.matrix(Xtest)
    if (ncol(Xtest) != ncol(Xrpm)) {
      stop("Column mismatch between RPM and Xtest and no names to align.")
    }
  }
  
  if (identical(k, "all")) k <- nrow(Xrpm)
  NNS_distance_bulk_cpp(Xrpm, as.numeric(rpm$y.hat), Xtest, as.integer(k), !is.null(class))
}



#' All-k bulk predictor (parallel if ncores > 1)
NNS.distance.path.bulk <- function(rpm, Xtest, kmax, class = NULL, ncores = 1L) {
  rpm   <- data.table::as.data.table(rpm)
  stopifnot("y.hat" %in% names(rpm))
  Xrpm  <- as.data.frame(rpm[, !"y.hat"])
  Xtest <- as.data.frame(Xtest)
  
  # Align by names if available
  if (!is.null(colnames(Xrpm)) && !is.null(colnames(Xtest))) {
    cmn <- intersect(colnames(Xrpm), colnames(Xtest))
    Xrpm  <- as.matrix(Xrpm[, cmn, drop = FALSE])
    Xtest <- as.matrix(Xtest[, cmn, drop = FALSE])
  } else {
    Xrpm  <- as.matrix(Xrpm)
    Xtest <- as.matrix(Xtest)
    if (ncol(Xrpm) != ncol(Xtest))
      stop("Column mismatch between RPM and Xtest and no names to align.")
  }
  
  if (identical(kmax, "all")) kmax <- nrow(Xrpm)
  kmax <- as.integer(kmax)
  is_class <- !is.null(class)
  
  if (ncores > 1L) {
    RcppParallel::setThreadOptions(numThreads = ncores)
    NNS_distance_path_parallel_cpp(Xrpm, as.numeric(rpm$y.hat), Xtest, kmax, is_class, ncores)
  } else {
    NNS_distance_path_cpp(Xrpm, as.numeric(rpm$y.hat), Xtest, kmax, is_class)
  }
}
