NNS.distance.bulk <- function(rpm, Xtest, k, class = NULL) {
  rpm <- data.table::as.data.table(rpm)
  stopifnot("y.hat" %in% names(rpm))
  
  # drop y.hat, align columns with Xtest by name if possible
  Xrpm <- as.data.frame(rpm[, !"y.hat"])
  if (!is.null(colnames(Xrpm)) && !is.null(colnames(Xtest))) {
    cmn <- intersect(colnames(Xrpm), colnames(Xtest))
    if (length(cmn) == 0L) {
      stop("No common feature columns between RPM and Xtest.")
    }
    Xrpm  <- as.matrix(Xrpm[, cmn, drop = FALSE])
    Xtest <- as.matrix(as.data.frame(Xtest)[, cmn, drop = FALSE])
  } else {
    Xrpm  <- as.matrix(Xrpm)
    Xtest <- as.matrix(Xtest)
    if (ncol(Xtest) != ncol(Xrpm)) {
      stop("Column mismatch between RPM and Xtest and no names to align.")
    }
  }
  
  if (identical(k, "all") ||
      (is.numeric(k) && length(k) == 1L && is.infinite(k))) {
    k <- nrow(Xrpm)
  } else {
    k <- suppressWarnings(as.integer(k[1L]))
    if (is.na(k)) k <- nrow(Xrpm)
    k <- max(1L, min(k, nrow(Xrpm)))
  }
  
  NNS_distance_bulk_cpp(Xrpm, as.numeric(rpm$y.hat), Xtest, as.integer(k), !is.null(class))
}

NNS.distance.path.bulk <- function(rpm, Xtest, kmax, class = NULL, ncores = 1L) {
  rpm   <- data.table::as.data.table(rpm)
  stopifnot("y.hat" %in% names(rpm))
  Xrpm  <- as.data.frame(rpm[, !"y.hat"])
  Xtest <- as.data.frame(Xtest)
  
  # Align by names if available
  if (!is.null(colnames(Xrpm)) && !is.null(colnames(Xtest))) {
    cmn <- intersect(colnames(Xrpm), colnames(Xtest))
    if (length(cmn) == 0L) {
      stop("No common feature columns between RPM and Xtest.")
    }
    Xrpm  <- as.matrix(Xrpm[, cmn, drop = FALSE])
    Xtest <- as.matrix(Xtest[, cmn, drop = FALSE])
  } else {
    Xrpm  <- as.matrix(Xrpm)
    Xtest <- as.matrix(Xtest)
    if (ncol(Xrpm) != ncol(Xtest)) {
      stop("Column mismatch between RPM and Xtest and no names to align.")
    }
  }
  
  if (identical(kmax, "all") ||
      (is.numeric(kmax) && length(kmax) == 1L && is.infinite(kmax))) {
    kmax <- nrow(Xrpm)
  } else {
    kmax <- suppressWarnings(as.integer(kmax[1L]))
    if (is.na(kmax)) kmax <- nrow(Xrpm)
    kmax <- max(1L, min(kmax, nrow(Xrpm)))
  }
  
  is_class <- !is.null(class)
  
  # Always use the parallel C++ routine to preserve the full multi-weight ensemble formulation matching NNS.distance
  RcppParallel::setThreadOptions(numThreads = as.integer(ncores))
  NNS_distance_path_parallel_cpp(Xrpm, as.numeric(rpm$y.hat), Xtest, kmax, is_class, as.integer(ncores))
}


NNS.distance.path.single.bulk <- function(rpm, Xtest, k, class = NULL, ncores = 1L) {
  rpm <- data.table::as.data.table(rpm)
  stopifnot("y.hat" %in% names(rpm))
  Xrpm <- as.data.frame(rpm[, !"y.hat"])
  
  if (is.null(dim(Xtest))) {
    Xtest <- as.data.frame(t(Xtest))
  } else {
    Xtest <- as.data.frame(Xtest)
  }
  
  # Align by names if available, matching NNS.distance.path.bulk.
  if (!is.null(colnames(Xrpm)) && !is.null(colnames(Xtest))) {
    cmn <- intersect(colnames(Xrpm), colnames(Xtest))
    if (length(cmn) == 0L) {
      stop("No common feature columns between RPM and Xtest.")
    }
    Xrpm  <- as.matrix(Xrpm[, cmn, drop = FALSE])
    Xtest <- as.matrix(Xtest[, cmn, drop = FALSE])
  } else {
    Xrpm  <- as.matrix(Xrpm)
    Xtest <- as.matrix(Xtest)
    if (ncol(Xrpm) != ncol(Xtest)) {
      stop("Column mismatch between RPM and Xtest and no names to align.")
    }
  }
  
  if (identical(k, "all") ||
      (is.numeric(k) && length(k) == 1L && is.infinite(k))) {
    k <- nrow(Xrpm)
  } else {
    k <- suppressWarnings(as.integer(k[1L]))
    if (is.na(k)) k <- nrow(Xrpm)
    k <- max(1L, min(k, nrow(Xrpm)))
  }
  
  is_class <- !is.null(class)
  
  RcppParallel::setThreadOptions(numThreads = as.integer(ncores))
  as.numeric(NNS_distance_path_single_parallel_cpp(
    Xrpm,
    as.numeric(rpm$y.hat),
    Xtest,
    as.integer(k),
    is_class,
    as.integer(ncores)
  ))
}
