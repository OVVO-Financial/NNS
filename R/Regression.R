# Internal validation helpers -------------------------------------------------

.nns_reg_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("[%s] must be a single TRUE or FALSE value.", name), call. = FALSE)
  }
  x
}


.nns_reg_response_vector <- function(y) {
  if (inherits(y, c("tbl", "data.table")) || is.data.frame(y)) {
    if (ncol(y) != 1L) stop("[y] must be a vector or one-column object.", call. = FALSE)
    return(y[[1L]])
  }
  if (is.matrix(y)) {
    if (ncol(y) != 1L) stop("[y] must be a vector or one-column object.", call. = FALSE)
    return(y[, 1L])
  }
  y
}

.nns_reg_type <- function(type, y) {
  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || is.na(type)) {
      stop("[type] must be NULL, 'REG', 'CLASS', or 'XONLY'.", call. = FALSE)
    }
    type <- tolower(type)
    if (!type %in% c("reg", "class", "xonly")) {
      stop("[type] must be NULL, 'REG', 'CLASS', or 'XONLY'.", call. = FALSE)
    }
  }
  if (identical(type, "reg") && !is.numeric(y)) {
    stop("[type = 'REG'] requires a numeric response.", call. = FALSE)
  }
  
  # The automatic test hashes y, so it only runs when no type was given.
  is.class <- identical(type, "class") ||
    (is.null(type) &&
       (is.factor(y) || is.character(y) || is.logical(y) ||
          (is.numeric(y) && length(unique(y)) == 2L &&
             setequal(sort(unique(as.numeric(y))), c(0, 1)))))
  is.xonly <- identical(type, "xonly")
  
  if (is.factor(y)) {
    class.levels <- levels(y)
    y.numeric <- as.numeric(y)
    class.values <- sort(unique(y.numeric))
  } else if (is.character(y)) {
    fy <- factor(y)
    class.levels <- levels(fy)
    y.numeric <- as.numeric(fy)
    class.values <- sort(unique(y.numeric))
  } else if (is.logical(y)) {
    y.numeric <- as.numeric(y)
    class.values <- sort(unique(y.numeric))
    class.levels <- as.character(class.values)
  } else {
    y.numeric <- as.numeric(y)
    class.values <- if (is.class) sort(unique(y.numeric)) else NULL
    class.levels <- if (is.class) as.character(class.values) else NULL
  }
  
  list(
    type = if (is.class) "class" else if (is.xonly) "xonly" else type,
    is.class = is.class,
    is.xonly = is.xonly,
    y = y.numeric,
    class.values = class.values,
    class.levels = class.levels
  )
}

.nns_reg_validate_order <- function(order) {
  if (is.null(order)) return(NULL)
  if (is.character(order) && length(order) == 1L &&
      !is.na(order) && tolower(order) == "max") return("max")
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
      order < 1 || order != floor(order)) {
    stop("[order] must be NULL, 'max', or a positive integer.", call. = FALSE)
  }
  as.integer(order)
}

.nns_reg_validate_nbest <- function(n.best) {
  # NA (as NNS.stack reports for a single predictor) selects the default.
  if (is.null(n.best) || (length(n.best) == 1L && is.na(n.best))) return(NULL)
  if (is.character(n.best) && length(n.best) == 1L &&
      !is.na(n.best) && tolower(n.best) == "all") return("all")
  if (!is.numeric(n.best) || length(n.best) != 1L || !is.finite(n.best) ||
      n.best < 1 || n.best != floor(n.best)) {
    stop("[n.best] must be NULL, 'all', or a positive integer.", call. = FALSE)
  }
  as.integer(n.best)
}

# NULL or NA means no screening (0), as does an NA entry of a vector. A vector
# holds one cutoff per class for classification, as NNS.stack reports.
.nns_reg_validate_threshold <- function(threshold) {
  if (is.null(threshold) || !length(threshold)) return(0)
  if (!is.numeric(threshold) && !all(is.na(threshold))) {
    stop("[threshold] must be NULL, NA, or finite nonnegative numbers.", call. = FALSE)
  }
  out <- as.numeric(threshold)
  names(out) <- names(threshold)
  out[is.na(out)] <- 0
  if (any(!is.finite(out)) || any(out < 0)) {
    stop("[threshold] must be NULL, NA, or finite nonnegative numbers.", call. = FALSE)
  }
  out
}

.nns_reg_validate_dist <- function(dist) {
  if (is.null(dist)) return("NNS")
  if (!is.character(dist) || length(dist) != 1L || is.na(dist)) {
    stop("[dist] must be NULL or one of 'NNS', 'L1', 'L2', or 'FACTOR'.", call. = FALSE)
  }
  dist <- toupper(dist)
  if (!dist %in% c("NNS", "L1", "L2", "FACTOR")) {
    stop("[dist] must be NULL or one of 'NNS', 'L1', 'L2', or 'FACTOR'.", call. = FALSE)
  }
  dist
}

.nns_reg_validate_noise <- function(noise.reduction) {
  if (!is.character(noise.reduction) || length(noise.reduction) != 1L ||
      is.na(noise.reduction)) {
    stop("[noise.reduction] must be one of 'mean', 'median', 'mode', or 'off'.",
         call. = FALSE)
  }
  noise.reduction <- tolower(noise.reduction)
  if (!noise.reduction %in% c("mean", "median", "mode", "off")) {
    stop("[noise.reduction] must be one of 'mean', 'median', 'mode', or 'off'.",
         call. = FALSE)
  }
  noise.reduction
}

.nns_reg_validate_ci <- function(confidence.interval) {
  if (is.null(confidence.interval)) return(NULL)
  if (!is.numeric(confidence.interval) || length(confidence.interval) != 1L ||
      !is.finite(confidence.interval) || confidence.interval <= 0 ||
      confidence.interval >= 1) {
    stop("[confidence.interval] must be NULL or a scalar strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(confidence.interval)
}

.nns_reg_as_frame <- function(x) {
  if (inherits(x, c("tbl", "data.table"))) x <- as.data.frame(x)
  if (is.null(dim(x))) {
    out <- data.frame(x = x, check.names = FALSE)
  } else {
    out <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (ncol(out) < 1L) stop("[x] must contain at least one predictor.", call. = FALSE)
  nm <- names(out)
  missing.names <- is.na(nm) | !nzchar(nm)
  nm[missing.names] <- paste0("x", which(missing.names))
  names(out) <- make.unique(nm, sep = ".")
  out
}

# Reconcile supplied prediction / test column names with the training predictor
# names.  Names only carry alignment information when both sides refer to the
# same predictors, so three cases are distinguished:
#
#   * same set of names  -> reorder the supplied columns into the training
#     order, which catches columns handed over in a different order;
#   * no names in common -> the names cannot express an alignment at all, as in
#     cbind(test.x_1, test.x_2) against cbind(x_1, x_2), so keep the supplied
#     order and match positionally per the documented "compatible dimensions";
#   * CSV BOM first header -> normalize an Excel/CSV UTF-8 import artifact
#     such as "X...Male" (or a literal BOM) before name-based matching;
#   * other overlap      -> ambiguous, and in practice a wrong, misspelled or
#     misplaced column, so error.
#
# [supplied] and [train.names] must already be make.unique() normalized and of
# equal length.  Returns the permutation to apply to the supplied columns, or
# NULL when they are to be taken in the order supplied.
.nns_match_predictor_names <- function(supplied, train.names, label) {
  if (setequal(supplied, train.names)) return(match(train.names, supplied))
  if (!any(supplied %in% train.names)) return(NULL)
  # A CSV UTF-8 byte-order mark can corrupt only the *first* header.
  # read.csv() sometimes turns "Male" into "X...Male"; other import paths
  # can retain the BOM itself. Repair this recognizable artifact only.
  # An arbitrary in-place rename is ambiguous and must not silently pass.
  repaired <- supplied
  if (startsWith(repaired[1L], "X...")) {
    repaired[1L] <- substring(repaired[1L], 5L)
  } else if (startsWith(repaired[1L], "\ufeff")) {
    repaired[1L] <- substring(repaired[1L], 2L)
  }
  if (setequal(repaired, train.names)) return(match(train.names, repaired))

  missing.names <- setdiff(train.names, supplied)
  extra.names <- setdiff(supplied, train.names)
  stop(sprintf(
    paste0("[%s] columns must exactly match the training predictors, or ",
           "share no names with them to be matched positionally. ",
           "Missing: %s; extra: %s."),
    label,
    if (length(missing.names)) paste(missing.names, collapse = ", ") else "none",
    if (length(extra.names)) paste(extra.names, collapse = ", ") else "none"
  ), call. = FALSE)
}

.nns_reg_prepare_points <- function(point.est, train.names) {
  if (is.null(point.est)) return(NULL)
  p <- length(train.names)
  
  if (p == 1L) {
    if (is.null(dim(point.est))) {
      out <- data.frame(point.est, check.names = FALSE)
      names(out) <- train.names
      return(out)
    }
    out <- as.data.frame(point.est, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(out) != 1L) {
      stop("[point.est] must contain exactly one predictor column.", call. = FALSE)
    }
    # A lone column can only align one way, so its name never selects between
    # predictors: cbind(test.x) against a single training predictor is matched
    # positionally.
    names(out) <- train.names
    return(out)
  }
  
  if (is.null(dim(point.est))) {
    if (length(point.est) != p) {
      stop(sprintf("A vector [point.est] must have exactly %d values.", p), call. = FALSE)
    }
    supplied <- names(point.est)
    if (!is.null(supplied) && all(nzchar(supplied))) {
      # Apply the same duplicate-name normalization used for the training frame.
      # Example: c(x = ..., x = ...) becomes c("x", "x.1") on both sides.
      supplied <- make.unique(supplied, sep = ".")
      ordering <- .nns_match_predictor_names(supplied, train.names, "point.est")
      if (!is.null(ordering)) point.est <- point.est[ordering]
    }
    out <- as.data.frame(as.list(point.est), check.names = FALSE,
                         stringsAsFactors = FALSE)
    names(out) <- train.names
    return(out)
  }
  
  supplied.names <- colnames(point.est)
  has.supplied.names <- !is.null(supplied.names) && all(nzchar(supplied.names))
  out <- as.data.frame(point.est, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(out) != p) {
    stop(sprintf("[point.est] must contain exactly %d predictor columns.", p), call. = FALSE)
  }
  
  if (has.supplied.names) {
    # Training names are normalized with make.unique() in .nns_reg_as_frame().
    # Normalize prediction names identically before validating/reordering so
    # cbind(x, x) matches training columns c("x", "x.1") by position.
    supplied.names <- make.unique(supplied.names, sep = ".")
    ordering <- .nns_match_predictor_names(supplied.names, train.names, "point.est")
    if (!is.null(ordering)) out <- out[, ordering, drop = FALSE]
  }
  names(out) <- train.names
  out
}

.nns_reg_encode_predictors <- function(x, point.est = NULL,
                                       factor.2.dummy = TRUE) {
  factor.2.dummy <- .nns_reg_scalar_logical(factor.2.dummy, "factor.2.dummy")
  # A plain numeric vector is the single numeric predictor "x", as
  # data.frame(x = x) below names it; building that data frame was most of the
  # encoding time for short series. Plain numeric points are its values.
  if (is.numeric(x) && is.null(dim(x)) && is.null(oldClass(x)) &&
      length(x) && all(is.finite(x)) &&
      (is.null(point.est) ||
         (is.numeric(point.est) && is.null(dim(point.est)) &&
            is.null(oldClass(point.est))))) {
    point.matrix <- NULL
    if (!is.null(point.est)) {
      if (any(!is.finite(point.est))) {
        stop("[point.est] predictor 'x' must contain only finite numeric values.",
             call. = FALSE)
      }
      point.matrix <- matrix(as.numeric(point.est), ncol = 1L,
                             dimnames = list(NULL, "x"))
    }
    return(list(
      x = matrix(as.numeric(x), ncol = 1L, dimnames = list(NULL, "x")),
      point.est = point.matrix,
      raw.names = "x",
      encoded.names = "x",
      metadata = list(list(name = "x", categorical = FALSE, levels = NULL))
    ))
  }
  train <- .nns_reg_as_frame(x)
  points <- .nns_reg_prepare_points(point.est, names(train))
  
  train.parts <- list()
  point.parts <- if (is.null(points)) NULL else list()
  meta <- vector("list", ncol(train))
  
  for (j in seq_along(train)) {
    nm <- names(train)[j]
    z <- train[[j]]
    zp <- if (is.null(points)) NULL else points[[j]]
    
    categorical <- is.factor(z) || is.character(z) || is.logical(z)
    
    if (categorical) {
      # A factor's declared levels are part of its training support even when
      # a resampled fold or subset has no row at one of them; sorting the
      # unordered levels keeps the encoding independent of level order.
      levels.train <- if (is.ordered(z)) levels(z) else
        if (is.factor(z)) sort(levels(z)) else sort(unique(as.character(z)))
      values.train <- as.character(z)
      if (anyNA(values.train)) {
        stop(sprintf("Predictor '%s' contains missing values.", nm), call. = FALSE)
      }
      if (!is.null(zp)) {
        values.point <- as.character(zp)
        if (anyNA(values.point)) {
          stop(sprintf("[point.est] predictor '%s' contains missing values.", nm),
               call. = FALSE)
        }
        unseen <- setdiff(unique(values.point), levels.train)
        if (length(unseen)) {
          stop(sprintf("[point.est] predictor '%s' contains unseen level(s): %s",
                       nm, paste(unseen, collapse = ", ")), call. = FALSE)
        }
      } else {
        values.point <- NULL
      }
      
      if (factor.2.dummy) {
        tr <- outer(values.train, levels.train, FUN = function(a, b) as.numeric(a == b))
        if (is.null(dim(tr))) tr <- matrix(tr, nrow = length(values.train),
                                           ncol = length(levels.train))
        colnames(tr) <- paste0(nm, "_", make.names(levels.train, unique = TRUE))
        train.parts[[length(train.parts) + 1L]] <- tr
        
        if (!is.null(points)) {
          pt <- outer(values.point, levels.train, FUN = function(a, b) as.numeric(a == b))
          if (is.null(dim(pt))) pt <- matrix(pt, nrow = length(values.point),
                                             ncol = length(levels.train))
          colnames(pt) <- colnames(tr)
          point.parts[[length(point.parts) + 1L]] <- pt
        }
      } else {
        tr <- matrix(match(values.train, levels.train), ncol = 1L,
                     dimnames = list(NULL, nm))
        train.parts[[length(train.parts) + 1L]] <- tr
        if (!is.null(points)) {
          pt <- matrix(match(values.point, levels.train), ncol = 1L,
                       dimnames = list(NULL, nm))
          point.parts[[length(point.parts) + 1L]] <- pt
        }
      }
      
      meta[[j]] <- list(name = nm, categorical = TRUE, levels = levels.train)
    } else {
      tr <- suppressWarnings(as.numeric(z))
      if (length(tr) != length(z) || any(!is.finite(tr))) {
        stop(sprintf("Predictor '%s' must contain only finite numeric values.", nm),
             call. = FALSE)
      }
      tr <- matrix(tr, ncol = 1L, dimnames = list(NULL, nm))
      train.parts[[length(train.parts) + 1L]] <- tr
      
      if (!is.null(points)) {
        pt <- if (is.factor(zp) || is.character(zp)) {
          suppressWarnings(as.numeric(as.character(zp)))
        } else {
          suppressWarnings(as.numeric(zp))
        }
        if (length(pt) != length(zp) || any(!is.finite(pt))) {
          stop(sprintf("[point.est] predictor '%s' must contain only finite numeric values.",
                       nm), call. = FALSE)
        }
        pt <- matrix(pt, ncol = 1L, dimnames = list(NULL, nm))
        point.parts[[length(point.parts) + 1L]] <- pt
      }
      meta[[j]] <- list(name = nm, categorical = FALSE, levels = NULL)
    }
  }
  
  x.matrix <- do.call(cbind, train.parts)
  storage.mode(x.matrix) <- "double"
  point.matrix <- if (is.null(points)) NULL else do.call(cbind, point.parts)
  if (!is.null(point.matrix)) storage.mode(point.matrix) <- "double"
  
  list(
    x = x.matrix,
    point.est = point.matrix,
    raw.names = names(train),
    encoded.names = colnames(x.matrix),
    metadata = meta
  )
}

.nns_reg_reduce_value <- function(z, noise.reduction, is.class = FALSE) {
  z <- z[is.finite(z)]
  if (!length(z)) return(NA_real_)
  if (is.class) {
    tab <- table(z)
    return(as.numeric(names(tab)[which.max(tab)]))
  }
  switch(noise.reduction,
         mean = mean(z),
         median = stats::median(z),
         mode = mode(z),
         off = gravity(z))
}

.nns_reg_snap_class <- function(x, class.values) {
  if (is.null(class.values) || !length(class.values)) return(x)
  vapply(x, function(z) {
    if (!is.finite(z)) return(NA_real_)
    class.values[which.min(abs(class.values - z))]
  }, numeric(1L))
}

# A vector with fewer than two distinct values, as length(unique(z)) < 2L but
# without hashing every value.
.nns_constant <- function(z) {
  if (length(z) < 2L) return(TRUE)
  if (anyNA(z)) return(length(unique(z)) < 2L)
  all(z == z[1L])
}

.nns_reg_dependence <- function(x, y) {
  if (.nns_constant(x) || .nns_constant(y)) return(0.1)
  # The order heuristic needs a fine-scale signal-to-noise reading, which the
  # local measure gives; NNS.dep's multi-level measure would also count
  # dependence between coarse cells and over-partition noisy data (#95).
  d1 <- tryCatch(.nns_dep_local(x, y, asym = TRUE),
                 error = function(e) NA_real_)
  d2 <- tryCatch({
    rescaled.x <- NNS.rescale(x, 0, 1)
    m <- cbind(rescaled.x, rescaled.x, NNS.rescale(y, 0, 1), deparse.level = 0)
    reflected <- m
    reflected[, 3L] <- 1 - reflected[, 3L]
    # Model complexity should not depend on the sign chosen for the response.
    # The directional copula itself is retained; both orientations inform this
    # unsigned signal-to-noise summary.
    mean(c(NNS.copula(m), NNS.copula(reflected)))
  }, error = function(e) NA_real_)
  d <- mean(c(d1, d2), na.rm = TRUE)
  if (!is.finite(d)) d <- 0.1
  min(1, max(0, d))
}

.nns_reg_default_order <- function(x, y) {
  dep <- .nns_reg_dependence(x, y)
  ord <- max(1L, as.integer(floor(dep * 10 + 0.5)))
  if (length(y) < 100L) ord <- max(1L, as.integer(floor(ord / 2)))
  ord
}

.nns_reg_build_points <- function(x, y, order, noise.reduction,
                                  is.class, xonly = FALSE) {
  reducer <- function(z) .nns_reg_reduce_value(z, noise.reduction, is.class)
  actual.order <- NA_integer_
  
  if (identical(order, "max")) {
    split.y <- split(y, x)
    rx <- as.numeric(names(split.y))
    ry <- unname(vapply(split.y, reducer, numeric(1L)))
    actual.order <- as.integer(length(rx))
  } else {
    ord <- if (is.null(order)) .nns_reg_default_order(x, y) else order
    nr <- if (is.class) "mode_class" else noise.reduction
    part <- NNS.part(x, y, type = "XONLY", order = ord, obs.req = 0,
                     min.obs.stop = TRUE, noise.reduction = nr)
    actual.order <- as.integer(part$order)
    rx <- part$regression.points$x
    ry <- part$regression.points$y
  }

  # The points are kept as two plain vectors: data frame subsetting and rbind()
  # cost more than the partition itself on short series, and NNS.ARMA calls
  # this once per lag and forecast step.
  keep <- is.finite(rx) & is.finite(ry)
  rx <- rx[keep]
  ry <- ry[keep]
  if (!length(rx)) stop("NNS regression produced no finite regression points.", call. = FALSE)
  # The order of strictly increasing values is the identity.
  if (is.unsorted(rx, strictly = TRUE)) {
    o <- order(rx, method = "radix")
    rx <- rx[o]
    ry <- ry[o]
  }
  # Equal x values (after the character rounding split() applies) share one
  # reduced point.
  merge.ties <- function() {
    groups <- split(ry, rx)
    rx <<- as.numeric(names(groups))
    ry <<- unname(vapply(groups, reducer, numeric(1L)))
    o <- order(rx, method = "radix")
    rx <<- rx[o]
    ry <<- ry[o]
  }
  if (anyDuplicated(rx)) merge.ties()

  # Always include training boundaries using the same response reducer.
  xmin <- min(x)
  xmax <- max(x)
  ymin <- reducer(y[x == xmin])
  ymax <- reducer(y[x == xmax])
  rx <- c(rx, xmin, xmax)
  ry <- c(ry, ymin, ymax)
  o <- order(rx, method = "radix")
  rx <- rx[o]
  ry <- ry[o]
  if (anyDuplicated(rx)) merge.ties()

  rp <- structure(list(x = unname(rx), y = unname(ry)),
                  row.names = c(NA_integer_, -length(rx)), class = "data.frame")
  # Preserve the order actually completed by NNS.part(), which can differ from
  # the requested order when stopping rules apply.  Plot titles must report
  # this realized order rather than the placeholder "auto".
  attr(rp, "nns.order") <- actual.order
  rp
}

.nns_reg_derivative <- function(rp) {
  if (nrow(rp) < 2L) {
    return(data.frame(Coefficient = 0,
                      X.Lower.Range = rp$x[1L],
                      X.Upper.Range = rp$x[1L]))
  }
  run <- diff(rp$x)
  slope <- ifelse(run == 0, 0, diff(rp$y) / run)
  .nns_plain_df(list(
    Coefficient = slope,
    X.Lower.Range = rp$x[-nrow(rp)],
    X.Upper.Range = rp$x[-1L]
  ))
}

# A data frame of equal-length, unnamed plain vectors (at least one row),
# exactly as data.frame() builds it, without data.frame()'s per-column
# coercion and name checks.
.nns_plain_df <- function(columns) {
  structure(columns, class = "data.frame",
            row.names = c(NA_integer_, -length(columns[[1L]])))
}

.nns_reg_predict_univariate <- function(xout, rp, smooth = FALSE,
                                        is.class = FALSE,
                                        class.values = NULL,
                                        smooth.fit = NULL) {
  if (!length(xout)) return(numeric())
  if (nrow(rp) == 1L) {
    pred <- rep(rp$y[1L], length(xout))
  } else if (isTRUE(smooth) && nrow(rp) >= 4L && !is.class) {
    if (is.null(smooth.fit)) {
      smooth.fit <- stats::smooth.spline(rp$x, rp$y)
    }
    pred <- stats::predict(smooth.fit, xout)$y
  } else {
    pred <- stats::approx(rp$x, rp$y, xout = xout, method = "linear",
                          rule = 2, ties = "ordered")$y
    left <- xout < min(rp$x)
    right <- xout > max(rp$x)
    if (any(left)) {
      slope <- (rp$y[2L] - rp$y[1L]) / (rp$x[2L] - rp$x[1L])
      pred[left] <- rp$y[1L] + (xout[left] - rp$x[1L]) * slope
    }
    if (any(right)) {
      n <- nrow(rp)
      slope <- (rp$y[n] - rp$y[n - 1L]) / (rp$x[n] - rp$x[n - 1L])
      pred[right] <- rp$y[n] + (xout[right] - rp$x[n]) * slope
    }
  }
  if (is.class) pred <- .nns_reg_snap_class(pred, class.values)
  as.numeric(pred)
}

.nns_reg_r2 <- function(actual, predicted) {
  # Racine-Hastie nonparametric goodness-of-fit: squared sample correlation
  # between observed and fitted values. Unlike predictive 1 - SSE / SST, this
  # within-sample descriptive measure is bounded in [0, 1].
  actual.centered <- actual - mean(actual)
  predicted.centered <- predicted - mean(predicted)
  denominator <- sum(actual.centered^2) * sum(predicted.centered^2)

  if (denominator == 0) {
    return(if (isTRUE(all.equal(actual, predicted))) 1 else 0)
  }

  r2 <- sum(actual.centered * predicted.centered)^2 / denominator
  # Protect the mathematical bounds against floating-point roundoff.
  min(1, max(0, r2))
}

.nns_reg_intervals <- function(actual, fitted, point.pred, confidence.interval,
                               is.class = FALSE, class.values = NULL) {
  if (is.null(confidence.interval)) {
    return(list(conf.lower = NULL, conf.upper = NULL, pred.int = NULL))
  }
  alpha <- 1 - confidence.interval
  errors <- actual - fitted
  q <- as.numeric(stats::quantile(errors, probs = c(alpha / 2, 1 - alpha / 2),
                                  na.rm = TRUE, names = FALSE, type = 8))
  conf.lower <- fitted + q[1L]
  conf.upper <- fitted + q[2L]
  
  pred.int <- NULL
  if (!is.null(point.pred)) {
    lower <- point.pred + q[1L]
    upper <- point.pred + q[2L]
    if (is.class) {
      lower <- .nns_reg_snap_class(lower, class.values)
      upper <- .nns_reg_snap_class(upper, class.values)
      lo <- pmin(lower, upper)
      hi <- pmax(lower, upper)
      lower <- lo
      upper <- hi
    }
    pred.int <- data.frame(pred.int.neg = lower, pred.int.pos = upper)
  }
  list(conf.lower = conf.lower, conf.upper = conf.upper, pred.int = pred.int)
}

.nns_reg_dimred_coefficients <- function(x, y, dim.red.method, tau,
                                         threshold) {
  p <- ncol(x)
  if (!is.numeric(dim.red.method)) {
    if (!is.character(dim.red.method) || length(dim.red.method) != 1L ||
        is.na(dim.red.method)) {
      stop("[dim.red.method] must be NULL, a supported method, or a numeric vector.",
           call. = FALSE)
    }
    method <- tolower(dim.red.method)
    if (!method %in% c("cor", "nns.dep", "nns.caus", "all", "equal")) {
      stop("Unsupported [dim.red.method].", call. = FALSE)
    }
  } else {
    if (length(dim.red.method) != p || any(!is.finite(dim.red.method))) {
      stop(sprintf("A numeric [dim.red.method] must contain exactly %d finite coefficients.", p),
           call. = FALSE)
    }
    if (!is.null(names(dim.red.method))) {
      if (anyDuplicated(names(dim.red.method)) ||
          !setequal(names(dim.red.method), colnames(x))) {
        stop(paste0("Named numeric [dim.red.method] coefficients must exactly ",
                    "match the encoded design columns."), call. = FALSE)
      }
      dim.red.method <- dim.red.method[colnames(x)]
    }
    method <- "numeric"
  }
  
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("[threshold] must be a single finite nonnegative number.", call. = FALSE)
  }
  
  cor.coef <- vapply(seq_len(p), function(j) {
    z <- suppressWarnings(stats::cor(x[, j], y, method = "spearman"))
    if (is.finite(z)) z else 0
  }, numeric(1L))
  
  dep.coef <- function() vapply(seq_len(p), function(j) {
    z <- tryCatch(NNS.dep(x[, j], y, print.map = FALSE, asym = TRUE)$Dependence,
                  error = function(e) 0)
    if (is.finite(z)) z else 0
  }, numeric(1L))
  
  caus.coef <- function() {
    if (is.null(tau)) tau.use <- "cs" else {
      if (!is.character(tau) || length(tau) != 1L || is.na(tau) ||
          !tolower(tau) %in% c("cs", "ts")) {
        stop("[tau] must be NULL, 'cs', or 'ts'.", call. = FALSE)
      }
      tau.use <- tolower(tau)
    }
    vapply(seq_len(p), function(j) {
      z <- tryCatch(Uni.caus(y, x[, j], tau = tau.use, plot = FALSE),
                    error = function(e) 0)
      if (is.finite(z)) z else 0
    }, numeric(1L))
  }
  
  coef <- switch(method,
                 cor = cor.coef,
                 nns.dep = dep.coef(),
                 nns.caus = caus.coef(),
                 all = rowMeans(cbind(caus.coef(), dep.coef(), cor.coef,
                                      rep(1, p))),
                 equal = rep(1, p),
                 numeric = as.numeric(dim.red.method))
  
  preserved <- coef
  coef[abs(coef) < threshold] <- 0
  if (!any(abs(coef) > 0)) {
    coef <- preserved
    if (!any(abs(coef) > 0)) coef <- rep(1, p)
  }
  coef
}

.nns_reg_dimreduce <- function(x, point.est, y, dim.red.method, tau, threshold) {
  coef <- .nns_reg_dimred_coefficients(x, y, dim.red.method, tau, threshold)
  mins <- apply(x, 2L, min)
  maxs <- apply(x, 2L, max)
  ranges <- maxs - mins
  
  normalize <- function(m) {
    out <- matrix(0.5, nrow = nrow(m), ncol = ncol(m),
                  dimnames = dimnames(m))
    active <- ranges > 0
    if (any(active)) {
      out[, active] <- sweep(sweep(m[, active, drop = FALSE], 2L,
                                   mins[active], "-"), 2L, ranges[active], "/")
    }
    out
  }
  
  nx <- normalize(x)
  np <- if (is.null(point.est)) NULL else normalize(point.est)
  denominator <- sum(abs(coef) > 0)
  if (denominator < 1L) stop("Dimension reduction retained no predictors.", call. = FALSE)
  
  x.star <- as.numeric((nx %*% coef) / denominator)
  point.star <- if (is.null(np)) NULL else as.numeric((np %*% coef) / denominator)
  equation <- data.frame(
    Variable = c(colnames(x), "DENOMINATOR"),
    Coefficient = c(coef, denominator),
    stringsAsFactors = FALSE
  )
  list(x.star = x.star, point.star = point.star,
       equation = equation, coefficient = coef,
       minimums = mins, maximums = maxs)
}



.nns_reg_partition_points_fast <- function(x, y, order = NULL,
                                           noise.reduction = "off",
                                           is.class = FALSE) {
  rp <- .nns_reg_build_points(as.numeric(x), as.numeric(y), order,
                              noise.reduction, is.class, xonly = TRUE)
  out <- sort(unique(as.numeric(rp[, 1L])))
  if (!length(out)) out <- sort(unique(as.numeric(x)))
  out
}

.nns_reg_univariate_fast <- function(train_x, train_y, test_x, order = NULL,
                                     noise.reduction = "off",
                                     is.class = FALSE,
                                     class.values = NULL) {
  rp <- .nns_reg_build_points(as.numeric(train_x), as.numeric(train_y),
                              order, noise.reduction, is.class)
  pred <- .nns_reg_predict_univariate(as.numeric(test_x), rp, smooth = FALSE,
                                      is.class = is.class,
                                      class.values = class.values)
  list(prediction = pred, regression.points = rp)
}


#' NNS Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x a vector, matrix or data frame of variables of numeric or factor data types.
#' @param y a numeric or factor vector with compatible dimensions to \code{x}.
#' @param factor.2.dummy logical; \code{TRUE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.
#' @param order integer; Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{(order = ...)} integer settings with \code{(noise.reduction = "off")}.  \code{(order = "max")} retains distinct predictor locations. A perfect training fit requires compatible responses at tied locations and interpolation at those locations; duplicate predictors with different responses cannot all be interpolated.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "all", "equal", \code{numeric vector}, NULL) method for determining synthetic X* coefficients (per Dana and Dawes (2004)).  Selection of a method automatically engages the dimension reduction regression.  The default is \code{NULL} for full multivariate regression.  \code{(dim.red.method = "NNS.dep")} uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for directional association weights.  \code{(dim.red.method = "cor")} uses Spearman rank correlation for weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.  \code{(dim.red.method = "equal")} uses unit weights.  Alternatively, user can specify a numeric vector of coefficients. Named vectors must exactly match the encoded predictor columns and are matched by name; unnamed vectors follow encoded column order.
#' @param tau options("ts", NULL); \code{NULL}(default) To be used in conjunction with \code{(dim.red.method = "NNS.caus")} or \code{(dim.red.method = "all")}.  If the regression is using time-series data, set \code{(tau = "ts")} to use the time-series lag convention.
#' @param type \code{NULL} (default).  To perform a classification, set to \code{(type = "CLASS")}.  Like a logistic regression, it is not necessary for target variable of two classes e.g. [0, 1]: a numeric response with exactly the values 0 and 1, or a factor, character or logical response, is classified automatically.  Use \code{type = "REG"} to fit a numeric response as regression, including a two-valued response, with one or more predictors. \code{(type = "XONLY")} is the single-predictor regression alias and keeps the default x-only partition.
#' @param point.est a numeric or factor vector with compatible dimensions to \code{x}.  Returns the fitted value \code{y.hat} for any value of \code{x}.  Columns are matched to \code{x} by name when the two share the same predictor names, and positionally when they share no names at all (as in \code{cbind(test.x_1, test.x_2)} against \code{cbind(x_1, x_2)}).  A CSV UTF-8 byte-order mark in the first header (for example, \code{Male} becoming \code{X...Male}) is repaired before matching names. Other partial overlaps, including misspelled predictors, are ambiguous and return an error.
#' @param location Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.
#' @param return.values logical; \code{TRUE} (default), set to \code{FALSE} in order to only display a regression plot and call values as needed.
#' @param plot logical; \code{TRUE} (default) To plot regression.
#' @param plot.regions logical; \code{FALSE} (default).  Generates 3d regions associated with each regression point for multivariate regressions.  Note, adds significant time to routine.
#' @param residual.plot logical; \code{TRUE} (default) To plot \code{y.hat} and \code{Y}.
#' @param confidence.interval numeric [0, 1]; \code{NULL} (default) Plots the associated confidence interval with the estimate and adds \code{conf.int.pos} / \code{conf.int.neg} to \code{Fitted.xy}.  The intervals are global: they are built from the type-8 quantiles of the in-sample residuals (no per-segment standard errors are returned).  Also applies the same level for the prediction intervals.
#' @param threshold  numeric [0, 1]; \code{(threshold = 0)} (default) Sets the threshold for dimension reduction of independent variables when \code{(dim.red.method)} is not \code{NULL}. Predictors whose coefficient magnitude is below \code{threshold} are dropped; \code{NULL} or \code{NA} means 0. For classification, a vector with one cutoff per class (named by class, or in class order) applies each cutoff to that class's index, as \link{NNS.stack} can report in \code{NNS.dim.red.threshold}.
#' @param n.best integer; \code{NULL} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at \code{floor((1 - D) * sqrt(nrow(x)))} (at least 1, at most the number of regression points), where D averages the copula summaries with response y and -y so the smoothing choice is unchanged by response reflection.  \code{(n.best = "all")} will select and weight all generated regression points.  Analogous to \code{k} in a
#' \code{k Nearest Neighbors} algorithm.  Different values of \code{n.best} are tested using cross-validation in \link{NNS.stack}. \code{NA} selects the default, as \code{NULL} does.
#' @param smooth logical; \code{FALSE} (default) Applies a smoothing spline instead of local linear fit to regression points.
#' @param noise.reduction the method of determining regression points options: ("mean", "median", "mode", "off"); In low signal:noise situations,\code{(noise.reduction = "mean")}  uses means for \link{NNS.dep} restricted partitions, \code{(noise.reduction = "median")} uses medians instead of means for \link{NNS.dep} restricted partitions, while \code{(noise.reduction = "mode")}  uses modes instead of means for \link{NNS.dep} restricted partitions.  \code{(noise.reduction = "off")} uses an overall central tendency measure for partitions.
#' @param dist options:(NULL, "NNS", "L1", "L2", "FACTOR") the method of distance calculation. \code{dist = NULL} is the default and selects the native blended NNS distance, \eqn{sum(abs(z) + z^2)}, over range-normalized coordinates. \code{dist = "NNS"} is an explicit alias for the default. \code{dist = "L2"} selects Euclidean distance; \code{dist = "L1"} selects Manhattan distance; \code{dist = "FACTOR"} uses a frequency.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#' @param multivariate.call Internal argument for multivariate regressions.
#' @param point.only Internal argument for abbreviated output.
#' @return UNIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:
#' \itemize{
#'  \item{\code{"R2"}} is squared Pearson correlation between fitted and observed responses (Racine--Hastie goodness of fit), or fitted classification accuracy;
#'
#'  \item{\code{"SE"}} returns the overall standard error of the estimate between \code{y} and \code{y.hat};
#'
#'  \item{\code{"Prediction.Accuracy"}} returns the proportion of fitted class decisions equal to the observed class;
#'
#'  \item{\code{"derivative"}} for the coefficient of the \code{x} and its applicable range;
#'
#'  \item{\code{"Point.est"}} for the predicted value generated;
#'  
#'  \item{\code{"pred.int"}} lower and upper prediction intervals for the \code{"Point.est"} returned using the \code{"confidence.interval"} provided;
#'  
#'  \item{\code{"regression.points"}} provides the points used in the regression equation for the given order of partitions;
#'
#'  \item{\code{"Fitted.xy"}} returns a \code{data.frame} of \code{x}, \code{y}, \code{y.hat}, \code{resid}, \code{NNS.ID}, \code{gradient};
#' }
#'
#'
#' MULTIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:
#' \itemize{
#'  \item{\code{"R2"}} is squared Pearson correlation between fitted and observed responses (Racine--Hastie goodness of fit), or fitted classification accuracy;
#'
#'  \item{\code{"equation"}} returns the numerator of the synthetic X* dimension reduction equation as a \code{data.frame} consisting of regressor and its coefficient.  Denominator is the number of non-zero coefficients (\code{abs(coefficient) > 0}), returned in last row of \code{equation} \code{data.frame}.
#'
#'  \item{\code{"x.star"}} returns the synthetic X* as a vector;
#'
#'  \item{\code{"rhs.partitions"}} returns the partition points for each regressor \code{x};
#'
#'  \item{\code{"RPM"}} provides the Regression Point Matrix, the points for each \code{x} used in the regression equation for the given order of partitions;
#'
#'  \item{\code{"Point.est"}} returns the predicted value generated;
#'  
#'  \item{\code{"pred.int"}} lower and upper prediction intervals for the \code{"Point.est"} returned using the \code{"confidence.interval"} provided;
#'
#'  \item{\code{"Fitted.xy"}} returns a \code{data.frame} of the regressors, \code{y}, \code{y.hat}, \code{NNS.ID} and \code{residuals} for full multivariate regression, or of the synthetic \code{x}, \code{y}, \code{y.hat}, \code{NNS.ID}, \code{gradient} and \code{residuals} for dimension reduction regression.  \code{confidence.interval} adds \code{conf.int.neg} and \code{conf.int.pos}.
#' }
#'
#' @details
#' For regression, \code{R2} measures association between fitted and observed responses: it is \eqn{cor(y, \hat y)^2}, with the package's finite constant-vector convention. It is not \eqn{1-SSE/SST}; a biased or incorrectly scaled prediction can have \code{R2 = 1}. Evaluate held-out prediction errors for predictive accuracy.
#' 
#' The default partition order combines a finite-partition dependence score with the average of the duplicated-predictor copula score under both response orientations. This complexity heuristic is invariant to reflecting the response; the underlying directional dependence and copula definitions are unchanged. A supplied order bypasses this heuristic.
#' 
#' Univariate regression linearly interpolates the selected regression points and extends the first or last segment beyond their range (unless a smoothing spline is requested). Endpoint slopes can be unstable with sparse or noisy boundary points. Multivariate distance weighting is a weighted average of regression-point responses, so it does not automatically continue a trend outside the training support. Neither construction supplies a general extrapolation guarantee. With the default \code{n.best = NULL}, purely categorical predictors use exact observed joint cells, and an unobserved joint combination uses its nearest cells. An explicit \code{n.best} weights that many nearest cells, as for numeric predictors. A declared factor level with no training rows is accepted in \code{point.est} and is placed by the same distance rule. A coordinate with exactly two observed values retains both terminal cells, including a dummy column in a mixed design. Other finite-order numeric coordinates retain interval cells; \code{order = "max"} uses exact joint support.
#' 
#' Classification forms an indicator for each observed class, retains the empirical class frequencies in each predictor cell, and interpolates or distance-weights those probability vectors. Class labels never enter numeric interpolation. \code{Point.probabilities}, \code{Fitted.probabilities}, and \code{RPM.probabilities} expose the corresponding vectors. Tied class decisions follow first appearance in the training response. Factor and character responses retain numeric output codes with \code{class.levels}; numeric responses retain their supplied values. For dimension reduction, one signed NNS index is fitted per class indicator, and the resulting probabilities are normalized across classes; \code{equation} is then a list indexed by class. \code{confidence.interval} returns \code{Prediction.set}, the smallest leading set of classes attaining the requested estimated probability mass. These are descriptive probability sets, not calibrated coverage guarantees; nominal classes have no numeric prediction interval.
#'
#' @note
#' \itemize{
#'  \item Please ensure \code{point.est} is of compatible dimensions to \code{x}, error message will ensue if not compatible.
#'
#'  \item Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].  For other numeric responses set \code{type = "CLASS"}; class values can be any distinct numbers.
#'
#'  \item For low signal:noise instances, increasing the dimension may yield better results using \code{NNS.stack(cbind(x,x), y, method = 1, ...)}.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#'
#' Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"  \doi{10.1007/s10614-017-9713-5}
#'
#' Vinod, H. and Viole, F. (2018) "Clustering and Curve Fitting by Line Segments"  \doi{10.20944/preprints201801.0090.v1}
#' 
#' Viole, F. (2020) "Partitional Estimation Using Partial Moments" \doi{10.2139/ssrn.3592491}
#' 
#' Dana, J., and Dawes, R. M. (2004). The Superiority of Simple Alternatives to Regression for Social Science Predictions. Journal of Educational and Behavioral Statistics, 29(3), 317–331.
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.reg(x, y)
#'
#' ## Manual {order} selection
#' NNS.reg(x, y, order = 2)
#'
#' ## Maximum {order} selection
#' NNS.reg(x, y, order = "max")
#'
#'
#' ## For Multiple Regression:
#' x <- cbind(rnorm(100), rnorm(100), rnorm(100)) ; y <- rnorm(100)
#' NNS.reg(x, y, point.est = c(.25, .5, .75))
#'
#' ## For Multiple Regression based on Synthetic X* (Dimension Reduction):
#' x <- cbind(rnorm(100), rnorm(100), rnorm(100)) ; y <- rnorm(100)
#' NNS.reg(x, y, point.est = c(.25, .5, .75), dim.red.method = "cor", ncores = 1)
#'
#' ## IRIS dataset examples:
#' # Dimension Reduction:
#' NNS.reg(iris[,1:4], iris[,5], dim.red.method = "cor", order = 5, ncores = 1)
#'
#' # Dimension Reduction using causal weights:
#' NNS.reg(iris[,1:4], iris[,5], dim.red.method = "NNS.caus", order = 5, ncores = 1)
#'
#' # Multiple Regression:
#' NNS.reg(iris[,1:4], iris[,5], order = 2, noise.reduction = "off")
#'
#' # Classification:
#' NNS.reg(iris[,1:4], iris[,5], point.est = iris[1:10, 1:4], type = "CLASS")$Point.est
#'
#' ## To call fitted values:
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.reg(x, y)$Fitted
#'
#' ## To call partial derivative (univariate regression only):
#' NNS.reg(x, y)$derivative
#' }
#' @export


NNS.reg <- function(x, y,
                    factor.2.dummy = TRUE, order = NULL,
                    dim.red.method = NULL, tau = NULL,
                    type = NULL,
                    point.est = NULL,
                    location = "top",
                    return.values = TRUE,
                    plot = TRUE, plot.regions = FALSE, residual.plot = TRUE,
                    confidence.interval = NULL,
                    threshold = 0,
                    n.best = NULL,
                    smooth = FALSE,
                    noise.reduction = "off",
                    dist = NULL,
                    ncores = NULL,
                    point.only = FALSE,
                    multivariate.call = FALSE) {
  
  # Capture the original calls before x and y are validated/coerced/reassigned.
  # Calling substitute(y) after y has been overwritten deparses the full numeric
  # response vector and produces the wall of vertical axis text seen in vignettes.
  # They are only deparsed when a plot is drawn.
  x.expr <- substitute(x)
  y.expr <- substitute(y)
  
  factor.2.dummy <- .nns_reg_scalar_logical(factor.2.dummy, "factor.2.dummy")
  return.values <- .nns_reg_scalar_logical(return.values, "return.values")
  plot <- .nns_reg_scalar_logical(plot, "plot")
  plot.regions <- .nns_reg_scalar_logical(plot.regions, "plot.regions")
  residual.plot <- .nns_reg_scalar_logical(residual.plot, "residual.plot")
  smooth <- .nns_reg_scalar_logical(smooth, "smooth")
  point.only <- .nns_reg_scalar_logical(point.only, "point.only")
  multivariate.call <- .nns_reg_scalar_logical(multivariate.call, "multivariate.call")
  order <- .nns_reg_validate_order(order)
  n.best <- .nns_reg_validate_nbest(n.best)
  threshold <- .nns_reg_validate_threshold(threshold)
  dist <- .nns_reg_validate_dist(dist)
  noise.reduction <- .nns_reg_validate_noise(noise.reduction)
  confidence.interval <- .nns_reg_validate_ci(confidence.interval)
  if (!plot) residual.plot <- FALSE
  
  y <- .nns_reg_response_vector(y)
  if (length(y) < 2L) stop("[y] must contain at least two observations.", call. = FALSE)
  if (anyNA(y)) stop("[y] contains missing values.", call. = FALSE)
  
  task <- .nns_reg_type(type, y)
  y.numeric <- task$y
  if (any(!is.finite(y.numeric))) stop("[y] must contain finite values.", call. = FALSE)
  
  encoded <- .nns_reg_encode_predictors(x, point.est, factor.2.dummy)
  if (nrow(encoded$x) != length(y.numeric)) {
    stop(sprintf("[x] has %d rows but [y] has %d values.",
                 nrow(encoded$x), length(y.numeric)), call. = FALSE)
  }
  
  if (task$is.xonly && ncol(encoded$x) != 1L) {
    stop("[type = 'XONLY'] is only valid for a single encoded predictor.", call. = FALSE)
  }
  
  if (task$is.class) {
    ans <- .nns_reg_classify(encoded, task, order, n.best, dist,
                             .nns_mreg_validate_cores(ncores), dim.red.method,
                             tau, threshold, plot = plot, point.only = point.only,
                             confidence.interval = confidence.interval)
    if (return.values) return(ans) else return(invisible(ans))
  }

  # Full multivariate regression.
  if (ncol(encoded$x) > 1L && is.null(dim.red.method)) {
    ans <- NNS.M.reg(
      X_n = encoded$x,
      Y = y.numeric,
      factor.2.dummy = FALSE,
      order = order,
      n.best = n.best,
      type = "REG",
      point.est = encoded$point.est,
      point.only = point.only,
      plot = plot,
      residual.plot = residual.plot,
      location = location,
      noise.reduction = noise.reduction,
      dist = dist,
      return.values = return.values,
      plot.regions = plot.regions,
      ncores = ncores,
      confidence.interval = confidence.interval,
      categorical = all(vapply(encoded$metadata, `[[`, logical(1), "categorical"))
    )
    ans$class.levels <- task$class.levels
    if (return.values) return(ans) else return(invisible(ans))
  }
  
  synthetic.x.equation <- NULL
  x.star <- NULL
  
  # Dimension reduction converts the encoded matrix to one training-fitted X*.
  if (ncol(encoded$x) > 1L) {
    if (length(threshold) > 1L) {
      stop(paste0("[threshold] must be a single number for regression; one cutoff ",
                  "per class applies to classification."), call. = FALSE)
    }
    dr <- .nns_reg_dimreduce(encoded$x, encoded$point.est, y.numeric,
                             dim.red.method, tau, threshold)
    ux <- dr$x.star
    up <- dr$point.star
    synthetic.x.equation <- dr$equation
    x.star <- data.frame(x = ux)
  } else {
    ux <- as.numeric(encoded$x[, 1L])
    up <- if (is.null(encoded$point.est)) NULL else
      as.numeric(encoded$point.est[, 1L])
  }
  
  rp <- .nns_reg_build_points(ux, y.numeric, order, noise.reduction,
                              task$is.class, task$is.xonly)
  actual.order <- attr(rp, "nns.order", exact = TRUE)
  if (is.null(actual.order) || length(actual.order) != 1L ||
      !is.finite(actual.order)) {
    actual.order <- if (identical(order, "max")) {
      length(unique(ux))
    } else if (is.numeric(order)) {
      as.integer(order)
    } else {
      .nns_reg_default_order(ux, y.numeric)
    }
  }
  actual.order <- as.integer(actual.order)
  
  if (multivariate.call) {
    return(.NNS.df(.nns_plain_df(list(x = rp$x, y = rp$y))))
  }
  
  # Restore the original NNS smoothing rule: the smoothing parameter is tied
  # to dependence rather than delegated to smooth.spline's generic default.
  smooth.condition <- isTRUE(smooth) && nrow(rp) >= 4L && !task$is.class
  smooth.fit <- NULL
  if (smooth.condition) {
    dependence <- .nns_reg_dependence(ux, y.numeric)
    smooth.fit <- stats::smooth.spline(
      x = rp$x,
      y = rp$y,
      spar = (dependence + 0.5) / 2
    )
  }
  
  # The fitted values are only needed for the fit table and for intervals.
  fitted.pred <- if (point.only && is.null(confidence.interval)) NULL else
    .nns_reg_predict_univariate(
      ux, rp, smooth = smooth.condition, is.class = task$is.class,
      class.values = task$class.values, smooth.fit = smooth.fit
    )
  point.pred <- if (is.null(up)) NULL else .nns_reg_predict_univariate(
    up, rp, smooth = smooth.condition, is.class = task$is.class,
    class.values = task$class.values, smooth.fit = smooth.fit
  )
  
  derivative <- .nns_reg_derivative(rp)
  intervals <- .nns_reg_intervals(
    y.numeric, fitted.pred, point.pred, confidence.interval,
    task$is.class, task$class.values
  )
  
  if (point.only) {
    out <- list(
      R2 = NULL,
      SE = NULL,
      Prediction.Accuracy = NULL,
      equation = .NNS.df(synthetic.x.equation),
      x.star = .NNS.df(x.star),
      derivative = .NNS.df(derivative),
      Point.est = point.pred,
      pred.int = .NNS.df(intervals$pred.int),
      regression.points = .NNS.df(rp),
      Fitted.xy = NULL,
      class.levels = task$class.levels
    )
    return(out)
  }
  
  ids <- findInterval(ux, rp$x, left.open = FALSE, rightmost.closed = TRUE)
  ids <- pmax(1L, pmin(ids, nrow(rp)))
  grad.idx <- findInterval(ux, derivative$X.Lower.Range,
                           left.open = FALSE, rightmost.closed = TRUE)
  grad.idx <- pmax(1L, pmin(grad.idx, nrow(derivative)))
  
  fitted <- .nns_plain_df(list(
    x = ux,
    y = y.numeric,
    y.hat = fitted.pred,
    NNS.ID = ids,
    gradient = derivative$Coefficient[grad.idx],
    residuals = fitted.pred - y.numeric
  ))
  
  metric <- if (task$is.class) mean(fitted.pred == y.numeric) else
    .nns_reg_r2(y.numeric, fitted.pred)
  se <- sqrt(mean((fitted.pred - y.numeric)^2))
  if (!is.null(intervals$conf.lower)) {
    fitted$conf.int.neg <- intervals$conf.lower
    fitted$conf.int.pos <- intervals$conf.upper
  }
  
  if (plot) {
    xlim <- range(c(ux, up), finite = TRUE)
    ylim <- range(c(y.numeric, fitted.pred, point.pred, rp$y,
                    intervals$conf.lower, intervals$conf.upper), finite = TRUE)
    
    # Report the order actually carried out by NNS.part(), not "auto" and
    # not merely the originally requested ceiling.
    plot.order <- actual.order
    
    graphics::plot(
      ux, y.numeric,
      pch = 1, lwd = 2, col = "steelblue",
      xlim = xlim, ylim = ylim,
      xlab = if (is.null(dim.red.method)) paste(deparse(x.expr), collapse = " ") else "Synthetic X*",
      ylab = paste(deparse(y.expr), collapse = " "),
      main = paste0("NNS Order = ", plot.order),
      mgp = c(2.5, 0.5, 0),
      cex.lab = 1.5,
      cex.main = 2
    )
    
    if (!is.null(intervals$conf.lower)) {
      o <- order(ux, method = "radix")
      graphics::polygon(
        c(ux[o], rev(ux[o])),
        c(intervals$conf.upper[o], rev(intervals$conf.lower[o])),
        col = grDevices::rgb(1, 192 / 255, 203 / 255, alpha = 0.375),
        border = NA
      )
    }
    
    graphics::points(rp$x, rp$y, col = "red", pch = 15)
    
    if (smooth.condition) {
      o <- order(ux, method = "radix")
      graphics::lines(ux[o], fitted.pred[o], col = "red", lwd = 2)
    } else {
      graphics::lines(rp$x, rp$y, col = "red", lwd = 2, lty = 2)
    }
    
    if (!is.null(up)) {
      graphics::points(up, point.pred, col = "green", pch = 18, cex = 1.5)
    }
    
    label <- if (task$is.class) {
      paste("Accuracy:", format(metric, digits = 4))
    } else {
      bquote(bold(R^2 == .(format(metric, digits = 4))))
    }
    graphics::legend(location, legend = label, bty = "n", y.intersp = 0.75)
  }
  
  out <- list(
    R2 = metric,
    SE = se,
    Prediction.Accuracy = if (task$is.class) metric else NULL,
    equation = .NNS.df(synthetic.x.equation),
    x.star = .NNS.df(x.star),
    derivative = .NNS.df(derivative),
    Point.est = point.pred,
    pred.int = .NNS.df(intervals$pred.int),
    regression.points = .NNS.df(rp),
    Fitted.xy = .NNS.df(fitted),
    class.levels = task$class.levels
  )
  
  if (return.values) out else invisible(out)
}
