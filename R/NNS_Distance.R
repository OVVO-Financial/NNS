#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param dist.estimate Query coordinates. Matching names select RPM features. If no names match, coordinates use the first numeric features in order. A partial match or a repeated matching name gives an error.
#' Each RPM row must share at least one finite coordinate with the query.
#' @param k \code{n.best} from \link{NNS.reg}
#' @param class if classification problem.
#'
#' @return The predicted response or class. Equal distances are ordered by RPM row.
#'
#'
#' @export


NNS.distance <- function(rpm, dist.estimate, k = "all", class = NULL) {
  rpm <- as.data.frame(rpm)
  if (!"y.hat" %in% names(rpm)) stop("rpm must contain column 'y.hat'")
  
  # 1) target vector
  dest <- unlist(dist.estimate, use.names = TRUE)
  n <- length(dest)
  y.hat <- as.numeric(rpm$y.hat)
  
  # 2) candidate feature columns, drop y.hat
  feat_all <- setdiff(names(rpm), "y.hat")
  
  # 3) choose columns to match dist.estimate
  feat <- .nns_distance_names(names(dest), feat_all)
  if (is.null(feat)) {
    # fall back: take the first n numeric columns, like the original
    numerics <- vapply(rpm[, feat_all, drop = FALSE], is.numeric, logical(1L))
    feat <- feat_all[numerics]
    if (length(feat) < n) stop("Not enough numeric feature columns in rpm")
    feat <- feat[seq_len(n)]
  } else {
    dest <- dest[feat]
  }

  X <- as.matrix(rpm[, feat, drop = FALSE])
  if (ncol(X) != n) {
    stop(sprintf(
      "after alignment, ncol(X)=%d != length(dist.estimate)=%d",
      ncol(X),
      n
    ))
  }
  
  # 4) k handling
  # Oversized k means use all available RPM rows.
  if (identical(k, "all") ||
      (is.numeric(k) && length(k) == 1L && is.infinite(k))) {
    k <- nrow(X)
  } else {
    k <- suppressWarnings(as.integer(k[1L]))
    if (is.na(k)) k <- nrow(X)
    k <- max(1L, min(k, nrow(X)))
  }
  
  # 5) call the C++ core
  NNS_distance_cpp(X, y.hat, as.numeric(dest), as.integer(k), !is.null(class))
}

.nns_distance_names <- function(supplied, available) {
  matched <- supplied %in% available
  if (!any(matched)) return(NULL)
  if (!all(matched))
    stop("Unknown distance feature: ", paste(supplied[!matched], collapse = ", "))
  if (anyDuplicated(supplied)) stop("Duplicate distance feature names")
  intersect(available, supplied)
}