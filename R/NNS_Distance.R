#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param dist.estimate Vector to generate distances from.
#' @param k \code{n.best} from \link{NNS.reg}
#' @param class if classification problem.
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export


NNS.distance <- function(rpm, dist.estimate, k = "all", class = NULL) {
  rpm <- data.table::as.data.table(rpm)
  if (!"y.hat" %in% names(rpm)) stop("rpm must contain column 'y.hat'")
  
  # 1) target vector
  dest <- unlist(dist.estimate, use.names = TRUE)
  n <- length(dest)
  y.hat <- as.numeric(rpm$y.hat)
  
  # 2) candidate feature columns (drop y.hat)
  feat_all <- setdiff(names(rpm), "y.hat")
  
  # 3) choose columns to match dist.estimate
  if (!is.null(names(dest)) && all(names(dest) %in% feat_all)) {
    # align by names (preferred)
    feat <- names(dest)
  } else {
    # fall back: take the first n *numeric* columns (like the original)
    numerics <- vapply(rpm[, ..feat_all], is.numeric, logical(1L))
    feat <- feat_all[numerics]
    if (length(feat) < n) stop("Not enough numeric feature columns in rpm")
    feat <- feat[seq_len(n)]
  }
  
  X <- as.matrix(rpm[, ..feat])
  if (ncol(X) != n) stop(sprintf("after alignment, ncol(X)=%d != length(dist.estimate)=%d", ncol(X), n))
  
  # 4) k handling
  if (identical(k, "all")) k <- nrow(X)
  k <- as.integer(k)
  
  # 5) call the C++ core
  NNS_distance_cpp(X, y.hat, as.numeric(dest), k, !is.null(class))
}
