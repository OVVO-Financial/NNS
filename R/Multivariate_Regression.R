.nns_mreg_weighted_mode <- function(values, weights) {
  keep <- is.finite(values) & is.finite(weights) & weights >= 0
  values <- values[keep]
  weights <- weights[keep]
  if (!length(values)) return(NA_real_)
  totals <- tapply(weights, values, sum)
  winners <- as.numeric(names(totals)[totals == max(totals)])
  min(winners)
}

.nns_mreg_normalize_weights <- function(w) {
  w[!is.finite(w) | w < 0] <- 0
  s <- sum(w)
  if (s <= 0) rep(1 / length(w), length(w)) else w / s
}

.nns_mreg_ensemble_weights <- function(distances) {
  k <- length(distances)
  if (k == 1L) return(1)
  
  ranks <- seq_len(k)
  uniform <- rep(1 / k, k)
  student <- .nns_mreg_normalize_weights(stats::dt(distances, df = k))
  inverse <- .nns_mreg_normalize_weights(1 / pmax(distances, 1e-12))
  exponential <- .nns_mreg_normalize_weights(stats::dexp(ranks, rate = 1 / k))
  
  rank.sd <- stats::sd(ranks)
  lognormal <- if (is.finite(rank.sd) && rank.sd > 0) {
    rev(.nns_mreg_normalize_weights(abs(stats::dlnorm(ranks, 0, rank.sd, log = TRUE))))
  } else rep(0, k)
  
  power <- .nns_mreg_normalize_weights(ranks^(-2))
  distance.sd <- stats::sd(distances)
  normal <- if (is.finite(distance.sd) && distance.sd > 0) {
    .nns_mreg_normalize_weights(stats::dnorm(distances, 0, distance.sd))
  } else rep(0, k)
  
  distance.var <- stats::var(distances)
  rbf <- if (is.finite(distance.var) && distance.var > 0) {
    .nns_mreg_normalize_weights(exp(-distances / (2 * distance.var)))
  } else rep(0, k)
  
  .nns_mreg_normalize_weights(
    uniform + student + inverse + exponential + lognormal + power + normal + rbf
  )
}

.nns_mreg_distance_vector <- function(rpm.x, destination, dist,
                                      minimums, maximums) {
  if (any(!is.finite(destination))) {
    stop("[point.est] contains missing or nonfinite values.", call. = FALSE)
  }
  
  if (dist == "FACTOR") {
    # Encoded categorical predictors are represented by stable training-fitted
    # dummy/code columns. Hamming distance is therefore the coherent categorical
    # metric and does not invent an ordering among unequal categories.
    return(rowMeans(sweep(rpm.x, 2L, destination, "!=") * 1))
  }
  
  ranges <- maximums - minimums
  active <- is.finite(ranges) & ranges > 0
  if (!any(active)) return(rep(0, nrow(rpm.x)))
  z <- sweep(rpm.x[, active, drop = FALSE], 2L, destination[active], "-")
  z <- sweep(z, 2L, ranges[active], "/")
  
  if (dist == "L1") rowSums(abs(z)) else sqrt(rowSums(z^2))
}

.nns_mreg_predict_one <- function(destination, rpm, k, dist,
                                  minimums, maximums, is.class) {
  rpm.x <- as.matrix(rpm[, setdiff(names(rpm), "y.hat"), drop = FALSE])
  d <- .nns_mreg_distance_vector(rpm.x, destination, dist, minimums, maximums)
  ord <- order(d, seq_along(d), method = "radix")
  k <- min(k, length(ord))
  
  # For k = 1, aggregate all exact nearest-distance ties deterministically.
  if (k == 1L) {
    tied <- which(d == min(d))
    vals <- rpm$y.hat[tied]
    if (is.class) return(.nns_mreg_weighted_mode(vals, rep(1, length(vals))))
    return(gravity(vals))
  }
  
  idx <- ord[seq_len(k)]
  w <- .nns_mreg_ensemble_weights(d[idx])
  if (is.class) .nns_mreg_weighted_mode(rpm$y.hat[idx], w) else
    sum(rpm$y.hat[idx] * w)
}

.nns_mreg_predict <- function(Xtest, rpm, k, dist,
                              minimums, maximums, is.class,
                              ncores = 1L) {
  if (is.null(Xtest)) return(NULL)
  Xtest <- as.matrix(Xtest)
  if (!nrow(Xtest)) return(numeric())
  if (ncol(Xtest) != length(minimums)) {
    stop("Prediction data and the fitted RPM have incompatible dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(Xtest))) {
    stop("[point.est] contains missing or nonfinite values.", call. = FALSE)
  }

  rpm.x <- as.matrix(rpm[, setdiff(names(rpm), "y.hat"), drop = FALSE])
  storage.mode(rpm.x) <- "double"
  storage.mode(Xtest) <- "double"
  dist.code <- match(dist, c("L2", "L1", "FACTOR")) - 1L

  as.numeric(if (isTRUE(getOption("NNS.native.mreg", TRUE))) {
    NNS_mreg_predict_v2_cpp(
      rpm.x, as.numeric(rpm$y.hat), Xtest, as.integer(k), dist.code,
      as.numeric(minimums), as.numeric(maximums), isTRUE(is.class),
      as.integer(ncores)
    )
  } else {
    NNS_mreg_predict_cpp(
      rpm.x, as.numeric(rpm$y.hat), Xtest, as.integer(k), dist.code,
      as.numeric(minimums), as.numeric(maximums), isTRUE(is.class)
    )
  })
}

# Reference pure-R implementation of the prediction rule, retained for
# equivalence tests against NNS_mreg_predict_cpp.
.nns_mreg_predict_reference <- function(Xtest, rpm, k, dist,
                                        minimums, maximums, is.class) {
  vapply(seq_len(nrow(Xtest)), function(i) .nns_mreg_predict_one(
    Xtest[i, ], rpm, k, dist, minimums, maximums, is.class
  ), numeric(1L))
}

.nns_mreg_group_reduce <- function(z, noise.reduction, is.class) {
  .nns_reg_reduce_value(z, noise.reduction, is.class)
}

.nns_mreg_build_rpm <- function(X, y, ids, noise.reduction, is.class) {
  # Fast path: every observation in its own cell (the common continuous case).
  # Each reducer returns a singleton's own value, so the RPM is the data
  # itself, ordered by interval ID exactly as split() would order the groups.
  if (!anyDuplicated(ids)) {
    o <- order(ids)
    rpm <- as.data.frame(X[o, , drop = FALSE], stringsAsFactors = FALSE)
    rpm$y.hat <- y[o]
    names(rpm) <- c(colnames(X), "y.hat")
    rownames(rpm) <- NULL
    return(rpm)
  }
  groups <- split(seq_len(nrow(X)), ids)
  rows <- lapply(groups, function(idx) {
    c(vapply(seq_len(ncol(X)), function(j) {
      .nns_mreg_group_reduce(X[idx, j], noise.reduction, FALSE)
    }, numeric(1L)),
    y.hat = .nns_mreg_group_reduce(y[idx], noise.reduction, is.class))
  })
  rpm <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(rpm) <- c(colnames(X), "y.hat")
  rownames(rpm) <- NULL
  rpm
}

.nns_mreg_partition_matrix <- function(X, y, order, noise.reduction, is.class) {
  p <- ncol(X)
  boundaries <- vector("list", p)

  # Exact duplicate predictor columns share one partition calculation.
  reps <- if (isTRUE(getOption("NNS.native.mreg", TRUE))) {
    tryCatch(
      as.integer(NNS_duplicate_column_map_cpp(as.matrix(X))),
      error = function(e) seq_len(p)
    )
  } else {
    seq_len(p)
  }

  for (j in seq_len(p)) {
    r <- reps[j]
    if (r < j && !is.null(boundaries[[r]])) {
      boundaries[[j]] <- boundaries[[r]]
      next
    }

    if (identical(order, "max")) {
      # Maximum order is the observed coordinate support itself.  No recursive
      # partitioning is required.
      boundaries[[j]] <- sort(unique(X[, j]))
    } else {
      boundaries[[j]] <- .nns_reg_partition_points_fast(
        X[, j], y,
        order = order,
        noise.reduction = noise.reduction,
        is.class = is.class
      )
    }
  }

  max.length <- max(lengths(boundaries))
  rhs <- do.call(cbind, lapply(boundaries, function(z) {
    length(z) <- max.length
    z
  }))
  colnames(rhs) <- colnames(X)

  # order = "max" is not an interval-boundary problem.  Each distinct observed
  # coordinate is a regression-point coordinate, so use its exact rank.  This
  # preserves a separate ID for the maximum value rather than folding it into
  # the preceding interval via rightmost.closed = TRUE.
  id.parts <- if (identical(order, "max")) {
    lapply(seq_len(p), function(j) {
      match(X[, j], boundaries[[j]])
    })
  } else {
    lapply(seq_len(p), function(j) {
      findInterval(
        X[, j], boundaries[[j]],
        left.open = FALSE,
        rightmost.closed = TRUE
      )
    })
  }

  ids <- do.call(paste, c(as.data.frame(id.parts), sep = "."))
  list(rhs = rhs, ids = ids, boundaries = boundaries)
}

.nns_mreg_prepare_model <- function(
    X,
    y,
    order = NULL,
    noise.reduction = "off",
    is.class = FALSE,
    use.native = isTRUE(getOption("NNS.native.mreg", TRUE))
) {
  X <- as.data.frame(X, check.names = FALSE)
  y <- as.numeric(y)
  minimums <- vapply(X, min, numeric(1L))
  maximums <- vapply(X, max, numeric(1L))
  partition <- .nns_mreg_partition_matrix(
    X, y, order, noise.reduction, is.class
  )

  # Maximum order already supplies the exact observed coordinates.  Bypass the
  # generic rightmost-closed native interval setup and build the RPM directly
  # from the exact rank IDs.
  if (identical(order, "max")) {
    rpm <- .nns_mreg_build_rpm(
      X, y, partition$ids, noise.reduction, is.class
    )
    row.ids <- if (!anyDuplicated(partition$ids)) {
      partition$ids[order(partition$ids, method = "radix")]
    } else {
      names(split(seq_len(nrow(X)), partition$ids))
    }

    return(list(
      RPM = rpm,
      rhs.partitions = partition$rhs,
      ids = partition$ids,
      row.ids = row.ids,
      boundaries = partition$boundaries,
      minimums = minimums,
      maximums = maximums
    ))
  }

  reducer.code <- switch(
    noise.reduction,
    mean = 0L,
    median = 1L,
    mode = 2L,
    off = 3L
  )

  if (use.native) {
    native <- tryCatch(
      NNS_mreg_setup_cpp(
        as.matrix(X), y, partition$boundaries,
        reducer.code, is.class
      ),
      error = function(e) NULL
    )

    if (!is.null(native)) {
      rpm <- as.data.frame(native$RPM, stringsAsFactors = FALSE)
      names(rpm) <- c(colnames(X), "y.hat")
      rownames(rpm) <- NULL
      row.ids <- if (!is.null(native$row_ids)) {
        as.character(native$row_ids)
      } else if (!anyDuplicated(partition$ids)) {
        partition$ids[order(partition$ids, method = "radix")]
      } else {
        names(split(seq_len(nrow(X)), partition$ids))
      }

      return(list(
        RPM = rpm,
        rhs.partitions = partition$rhs,
        ids = partition$ids,
        row.ids = row.ids,
        boundaries = partition$boundaries,
        minimums = minimums,
        maximums = maximums
      ))
    }
  }

  rpm <- .nns_mreg_build_rpm(
    X, y, partition$ids, noise.reduction, is.class
  )
  row.ids <- if (!anyDuplicated(partition$ids)) {
    partition$ids[order(partition$ids, method = "radix")]
  } else {
    names(split(seq_len(nrow(X)), partition$ids))
  }

  list(
    RPM = rpm,
    rhs.partitions = partition$rhs,
    ids = partition$ids,
    row.ids = row.ids,
    boundaries = partition$boundaries,
    minimums = minimums,
    maximums = maximums
  )
}

.nns_mreg_default_nbest <- function(X, y, rpm.rows) {
  dep <- tryCatch(NNS.copula(cbind(X, y)), error = function(e) NA_real_)
  if (!is.finite(dep)) dep <- 0.5
  k <- max(1L, as.integer(floor((1 - dep) * sqrt(nrow(X)))))
  min(k, rpm.rows)
}

.nns_mreg_validate_cores <- function(ncores) {
  if (is.null(ncores)) return(1L)
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || ncores != floor(ncores)) {
    stop("[ncores] must be NULL or a positive integer.", call. = FALSE)
  }
  as.integer(ncores)
}

# Internal multivariate regression engine used by NNS.reg.
NNS.M.reg <- function(X_n, Y, factor.2.dummy = TRUE, order = NULL,
                      n.best = NULL, type = NULL, point.est = NULL,
                      point.only = FALSE, plot = FALSE,
                      residual.plot = TRUE, location = NULL,
                      noise.reduction = "off", dist = "L2",
                      return.values = FALSE, plot.regions = FALSE,
                      ncores = NULL, confidence.interval = NULL) {
  
  factor.2.dummy <- .nns_reg_scalar_logical(factor.2.dummy, "factor.2.dummy")
  point.only <- .nns_reg_scalar_logical(point.only, "point.only")
  plot <- .nns_reg_scalar_logical(plot, "plot")
  residual.plot <- .nns_reg_scalar_logical(residual.plot, "residual.plot")
  return.values <- .nns_reg_scalar_logical(return.values, "return.values")
  plot.regions <- .nns_reg_scalar_logical(plot.regions, "plot.regions")
  order <- .nns_reg_validate_order(order)
  n.best <- .nns_reg_validate_nbest(n.best)
  dist <- .nns_reg_validate_dist(dist)
  noise.reduction <- .nns_reg_validate_noise(noise.reduction)
  confidence.interval <- .nns_reg_validate_ci(confidence.interval)
  ncores <- .nns_mreg_validate_cores(ncores)
  
  if (!is.null(type) && (!is.character(type) || length(type) != 1L ||
                         is.na(type) || tolower(type) != "class")) {
    stop("NNS.M.reg [type] must be NULL or 'CLASS'.", call. = FALSE)
  }
  Y <- .nns_reg_response_vector(Y)
  task <- .nns_reg_type(if (is.null(type)) NULL else "CLASS", Y)
  is.class <- task$is.class
  
  encoded <- .nns_reg_encode_predictors(X_n, point.est, factor.2.dummy)
  X <- encoded$x
  Xpoint <- encoded$point.est
  y <- task$y
  
  if (nrow(X) != length(y)) {
    stop(sprintf("[X_n] has %d rows but [Y] has %d values.",
                 nrow(X), length(y)), call. = FALSE)
  }
  if (ncol(X) < 2L) {
    stop("NNS.M.reg requires at least two encoded predictors.", call. = FALSE)
  }
  if (length(y) < 2L || any(!is.finite(y))) {
    stop("[Y] must contain at least two finite values.", call. = FALSE)
  }
  
  prepared <- .nns_mreg_prepare_model(
    X, y, order, noise.reduction, is.class
  )
  minimums <- prepared$minimums
  maximums <- prepared$maximums
  partition <- list(
    rhs = prepared$rhs.partitions,
    ids = prepared$ids,
    boundaries = prepared$boundaries
  )
  rpm <- prepared$RPM
  if (!nrow(rpm)) stop("NNS.M.reg produced an empty RPM.", call. = FALSE)
  
  if (is.null(n.best) && identical(order, "max")) {
    k <- 1L
  } else if (is.null(n.best)) {
    k <- .nns_mreg_default_nbest(X, y, nrow(rpm))
  } else if (identical(n.best, "all")) {
    k <- nrow(rpm)
  } else {
    k <- min(as.integer(n.best), nrow(rpm))
  }
  k <- max(1L, k)
  
  # Maximum order with the default k = 1 fit is an exact limit condition.
  # Every unique training row is already its own regression point; duplicate
  # rows share one reduced regression point.  Reuse those values directly
  # instead of performing nrow(X) x nrow(RPM) distance comparisons merely to
  # rediscover the identity mapping.
  if (identical(order, "max") && k == 1L) {
    if (!anyDuplicated(prepared$ids)) {
      fitted.pred <- y
    } else {
      rpm.index <- match(prepared$ids, prepared$row.ids)
      if (anyNA(rpm.index)) {
        stop(
          "Maximum-order RPM IDs could not be mapped back to training rows.",
          call. = FALSE
        )
      }
      fitted.pred <- rpm$y.hat[rpm.index]
    }
  } else {
    fitted.pred <- .nns_mreg_predict(
      X, rpm, k, dist, minimums, maximums, is.class, ncores
    )
  }

  # External point estimates still require the repaired distance rule.  With
  # point.est = NULL this returns immediately.
  point.pred <- .nns_mreg_predict(
    Xpoint, rpm, k, dist, minimums, maximums, is.class, ncores
  )
  
  if (is.class) {
    observed <- sort(unique(y))
    fitted.pred <- .nns_reg_snap_class(fitted.pred, observed)
    if (!is.null(point.pred)) point.pred <- .nns_reg_snap_class(point.pred, observed)
  }
  
  fitted <- as.data.frame(X, stringsAsFactors = FALSE)
  fitted$y <- y
  fitted$y.hat <- fitted.pred
  fitted$NNS.ID <- partition$ids
  fitted$residuals <- fitted.pred - y
  
  metric <- if (is.class) mean(fitted.pred == y) else
    .nns_reg_r2(y, fitted.pred)
  intervals <- .nns_reg_intervals(
    y, fitted.pred, point.pred, confidence.interval,
    is.class, if (is.class) sort(unique(y)) else NULL
  )
  if (!is.null(intervals$conf.lower)) {
    fitted$conf.int.neg <- intervals$conf.lower
    fitted$conf.int.pos <- intervals$conf.upper
  }
  
  if (point.only) {
    return(list(
      R2 = NULL,
      rhs.partitions = .NNS.df(as.data.frame(partition$rhs)),
      RPM = .NNS.df(rpm),
      Point.est = point.pred,
      pred.int = .NNS.df(intervals$pred.int),
      Fitted.xy = NULL,
      n.best = k,
      dist = dist,
      class.levels = task$class.levels
    ))
  }
  
  if (plot && ncol(X) == 2L) {
    .nns_require_rgl()
    rgl::plot3d(X[, 1L], X[, 2L], y, box = FALSE, size = 3,
                col = "steelblue", xlab = colnames(X)[1L],
                ylab = colnames(X)[2L], zlab = "Y")
    rgl::points3d(rpm[[1L]], rpm[[2L]], rpm$y.hat,
                  col = "red", size = 5)
    if (!is.null(Xpoint)) {
      rgl::points3d(Xpoint[, 1L], Xpoint[, 2L], point.pred,
                    col = "green", size = 5)
    }
    
    if (plot.regions) {
      for (id in unique(partition$ids)) {
        idx <- partition$ids == id
        x1 <- range(X[idx, 1L])
        x2 <- range(X[idx, 2L])
        z <- .nns_mreg_group_reduce(y[idx], noise.reduction, is.class)
        rgl::quads3d(
          x = c(x1[1L], x1[1L], x1[2L], x1[2L]),
          y = c(x2[1L], x2[2L], x2[2L], x2[1L]),
          z = rep(z, 4L), col = "pink", alpha = 0.6
        )
      }
    }
  }
  
  if (residual.plot) {
    graphics::plot(seq_along(y), y, pch = 1, lwd = 2,
                   col = "steelblue", xlab = "Index",
                   ylab = expression(paste("y (blue)   ", hat(y), " (red)")))
    graphics::lines(seq_along(fitted.pred), fitted.pred,
                    col = "red", lwd = 2)
    if (!is.null(intervals$conf.lower)) {
      idx <- seq_along(y)
      graphics::polygon(c(idx, rev(idx)),
                        c(intervals$conf.upper, rev(intervals$conf.lower)),
                        col = grDevices::rgb(1, 192 / 255, 203 / 255, alpha = 0.375),
                        border = NA)
    }
    label <- if (is.class) paste("Accuracy:", format(metric, digits = 4)) else
      bquote(bold(R^2 == .(format(metric, digits = 4))))
    graphics::legend(if (is.null(location)) "top" else location,
                     legend = label, bty = "n")
  }
  
  out <- list(
    R2 = metric,
    rhs.partitions = .NNS.df(as.data.frame(partition$rhs)),
    RPM = .NNS.df(rpm),
    Point.est = point.pred,
    pred.int = .NNS.df(intervals$pred.int),
    Fitted.xy = .NNS.df(fitted),
    n.best = k,
    dist = dist,
    class.levels = task$class.levels
  )
  
  if (return.values) out else invisible(out)
}
