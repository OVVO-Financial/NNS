# Nominal classification uses empirical class frequencies in NNS cells.
# Class numbers only identify columns; arithmetic is performed on indicators.

.nns_class_normalize <- function(p, prior) {
  p <- as.matrix(p)
  p[!is.finite(p) | p < 0] <- 0
  p[p > 1] <- 1
  total <- rowSums(p)
  good <- total > 0
  if (any(good)) p[good, ] <- p[good, , drop = FALSE] / total[good]
  if (any(!good)) p[!good, ] <- matrix(prior, sum(!good), length(prior), byrow = TRUE)
  p
}

.nns_class_select <- function(p, values, tie.order = seq_along(values)) {
  if (is.null(p)) return(NULL)
  p <- as.matrix(p)
  # A tied decision uses first appearance in the training response, not the
  # smallest class code. Relabelling or reordering factor levels preserves it.
  values[tie.order[max.col(p[, tie.order, drop = FALSE], ties.method = "first")]]
}

.nns_class_groups <- function(y) {
  lapply(unique(y), function(value) which(y == value))
}

.nns_class_strength <- function(x, probabilities) {
  scores <- vapply(seq_len(ncol(probabilities)), function(j) {
    .nns_reg_dependence(x, probabilities[, j])
  }, numeric(1))
  max(scores, 0)
}

.nns_class_base <- function(X, probabilities, point.est, order, n.best,
                            dist, ncores, categorical = FALSE,
                            need.fitted = TRUE, k.path = NULL) {
  X <- as.matrix(X)
  p <- ncol(X)
  x.min <- apply(X, 2, min)
  x.max <- apply(X, 2, max)
  exact.coordinate <- vapply(seq_len(p), function(j) {
    categorical || identical(order, "max") || length(unique(X[, j])) <= 2L
  }, logical(1))
  strength <- vapply(seq_len(p), function(j) {
    .nns_class_strength(X[, j], probabilities)
  }, numeric(1))
  boundaries <- lapply(seq_len(p), function(j) {
    if (exact.coordinate[j]) return(sort(unique(X[, j])))
    depth <- if (is.null(order)) {
      d <- max(1L, as.integer(floor(strength[j] * 10 + .5)))
      if (nrow(X) < 100L) d <- max(1L, as.integer(floor(d / 2)))
      d
    } else order
    .nns_reg_partition_points_fast(X[, j], X[, j], order = depth,
                                   noise.reduction = "mean", is.class = FALSE)
  })
  id.parts <- lapply(seq_len(p), function(j) {
    # A two-valued coordinate already has two terminal cells. Treating those
    # coordinates as interval endpoints would merge both into the same cell.
    if (exact.coordinate[j]) match(X[, j], boundaries[[j]]) else
      findInterval(X[, j], boundaries[[j]], rightmost.closed = TRUE)
  })
  ids <- do.call(paste, c(id.parts, sep = "."))
  groups <- split(seq_len(nrow(X)), ids)
  centres <- do.call(rbind, lapply(groups, function(g) colMeans(X[g, , drop = FALSE])))
  shares <- do.call(rbind, lapply(groups, function(g) colMeans(probabilities[g, , drop = FALSE])))
  colnames(centres) <- colnames(X)
  colnames(shares) <- colnames(probabilities)
  prior <- colMeans(probabilities)
  k <- if (identical(n.best, "all")) nrow(centres) else if (!is.null(n.best)) n.best else
    if (identical(order, "max") || categorical) 1L else
      max(1L, floor((1 - mean(strength)) * sqrt(nrow(centres))))
  k <- max(1L, min(as.integer(k), nrow(centres)))
  neighbours <- function(kk) {
    kk <- if (identical(kk, "all")) nrow(centres) else as.integer(kk)
    max(1L, min(kk, nrow(centres)))
  }
  # One probability matrix per neighbour count in ks. Each point's distances
  # are computed and sorted once and reused for every count.
  predict <- function(points, ks = k) {
    if (is.null(points)) return(NULL)
    points <- as.matrix(points)
    ks <- vapply(ks, neighbours, integer(1))
    # Counts above the number of cells all select every cell.
    unique.ks <- unique(ks)
    if (p == 1L && !categorical) {
      o <- order(centres[, 1])
      ans <- vapply(seq_len(ncol(shares)), function(j) {
        if (length(o) == 1L) rep(shares[o, j], nrow(points)) else
          stats::approx(centres[o, 1], shares[o, j], xout = points[, 1],
                        rule = 2, ties = mean)$y
      }, numeric(nrow(points)))
      ans <- matrix(ans, ncol = ncol(shares))
      colnames(ans) <- colnames(shares)
      return(rep(list(.nns_class_normalize(ans, prior)), length(ks)))
    }
    out <- lapply(unique.ks, function(kk) {
      matrix(0, nrow(points), ncol(shares), dimnames = list(NULL, colnames(shares)))
    })
    rank.terms <- vector("list", nrow(centres))
    for (i in seq_len(nrow(points))) {
      d <- .nns_mreg_distance_vector(centres, points[i, ], dist, x.min, x.max)
      # Exact matches retain their empirical class frequencies.
      exact <- which(d == 0)
      if (length(exact)) {
        v <- colMeans(shares[exact, , drop = FALSE])
        for (a in seq_along(unique.ks)) out[[a]][i, ] <- v
        next
      }
      # The k nearest cells, with every cell tied at the k-th distance.
      # Equal-distance cells share their rank weight.
      o <- order(d)
      ds <- d[o]
      for (a in seq_along(unique.ks)) {
        m <- sum(ds <= ds[unique.ks[a]])
        near <- ds[seq_len(m)]
        if (m > 1L && is.null(rank.terms[[m]])) rank.terms[[m]] <- .nns_mreg_rank_weights(m)
        weights <- .nns_mreg_ensemble_weights(near, if (m > 1L) rank.terms[[m]])
        for (v in unique(near[duplicated(near)])) {
          tied <- which(near == v)
          weights[tied] <- mean(weights[tied])
        }
        out[[a]][i, ] <- as.numeric(crossprod(weights, shares[o[seq_len(m)], , drop = FALSE]))
      }
    }
    lapply(out, .nns_class_normalize, prior = prior)[match(ks, unique.ks)]
  }
  list(fitted = if (need.fitted) predict(X)[[1L]] else NULL,
       predicted = if (is.null(k.path)) predict(point.est)[[1L]] else NULL,
       predicted.path = if (is.null(k.path)) NULL else predict(point.est, k.path),
       centres = centres, shares = shares, ids = ids, boundaries = boundaries,
       k = k, strength = strength)
}

.nns_reg_classify <- function(encoded, task, order, n.best, dist, ncores,
                              dim.red.method = NULL, tau = NULL, threshold = 0,
                              dim.count = NULL, plot = FALSE,
                              point.only = FALSE, confidence.interval = NULL) {
  X <- encoded$x
  Q <- encoded$point.est
  values <- task$class.values
  labels <- if (length(task$class.levels) == length(values)) task$class.levels else
    task$class.levels[values]
  indicators <- outer(task$y, values, "==") * 1
  colnames(indicators) <- labels
  tie.order <- match(unique(task$y), values)
  categorical <- all(vapply(encoded$metadata, `[[`, logical(1), "categorical"))
  equations <- NULL
  if (!is.null(dim.red.method) && ncol(X) > 1L) {
    # Each nominal class has its own signed NNS indicator index. A single
    # index fitted to numeric class codes would impose an arbitrary ordering.
    models <- lapply(seq_along(values), function(j) {
      coef <- .nns_reg_dimred_coefficients(X, indicators[, j], dim.red.method,
                                          tau, threshold)
      if (!is.null(dim.count) && dim.count < length(coef)) {
        keep <- order(abs(coef), decreasing = TRUE)[seq_len(dim.count)]
        coef[-keep] <- 0
      }
      if (!any(abs(coef) > 0)) coef[] <- 1
      reduced <- .nns_reg_dimreduce(X, Q, indicators[, j], coef, tau, 0)
      fit <- .nns_class_base(matrix(reduced$x.star, ncol = 1),
                             cbind(1 - indicators[, j], indicators[, j]),
                             if (is.null(Q)) NULL else matrix(reduced$point.star, ncol = 1),
                             order, n.best, dist, ncores, need.fitted = !point.only)
      list(fitted = if (point.only) NULL else fit$fitted[, 2],
           predicted = if (is.null(Q)) NULL else fit$predicted[, 2],
           equation = reduced$equation)
    })
    fitted.prob <- if (point.only) NULL else
      .nns_class_normalize(do.call(cbind, lapply(models, `[[`, "fitted")), colMeans(indicators))
    point.prob <- if (is.null(Q)) NULL else
      .nns_class_normalize(do.call(cbind, lapply(models, `[[`, "predicted")), colMeans(indicators))
    equations <- setNames(lapply(models, `[[`, "equation"), labels)
    centres <- shares <- boundaries <- NULL
    ids <- seq_len(nrow(X))
    k <- NA_integer_
  } else {
    fit <- .nns_class_base(X, indicators, Q, order, n.best, dist, ncores, categorical,
                           need.fitted = !point.only)
    fitted.prob <- fit$fitted; point.prob <- fit$predicted
    centres <- fit$centres; shares <- fit$shares
    boundaries <- fit$boundaries; ids <- fit$ids; k <- fit$k
  }
  if (!is.null(point.prob)) colnames(point.prob) <- labels
  point.pred <- .nns_class_select(point.prob, values, tie.order)
  # Fitted classes are not computed for point-only calls, which return none.
  fitted <- fitted.pred <- accuracy <- NULL
  if (!point.only) {
    colnames(fitted.prob) <- labels
    fitted.pred <- .nns_class_select(fitted.prob, values, tie.order)
    fitted <- as.data.frame(X)
    fitted$y <- task$y; fitted$y.hat <- fitted.pred
    fitted$NNS.ID <- ids; fitted$error <- as.numeric(fitted.pred != task$y)
    accuracy <- mean(fitted.pred == task$y)
  }
  if (plot && !point.only && ncol(X) == 1L) {
    graphics::plot(X[, 1], task$y, xlab = colnames(X)[1], ylab = "Class", pch = 1)
    graphics::points(X[, 1], fitted.pred, col = "red", pch = 20)
    if (!is.null(Q)) graphics::points(Q[, 1], point.pred, col = "green", pch = 18)
  }
  prediction.set <- if (is.null(point.prob) || is.null(confidence.interval)) NULL else
    lapply(seq_len(nrow(point.prob)), function(i) {
      o <- tie.order[order(point.prob[i, tie.order], decreasing = TRUE)]
      keep <- seq_len(which(cumsum(point.prob[i, o]) >= confidence.interval - 1e-12)[1])
      values[o[keep]]
    })
  .NNS.out(list(R2 = if (point.only) NULL else accuracy, SE = NULL,
    Prediction.Accuracy = if (point.only) NULL else accuracy,
    equation = equations, x.star = NULL, derivative = NULL,
    Point.est = point.pred, pred.int = NULL, Prediction.set = prediction.set,
    regression.points = if (is.null(centres)) NULL else as.data.frame(centres),
    RPM = if (is.null(centres)) NULL else as.data.frame(centres),
    RPM.probabilities = shares, rhs.partitions = boundaries,
    Fitted.xy = if (point.only) NULL else fitted, Fitted = if (point.only) NULL else fitted,
    Fitted.probabilities = if (point.only) NULL else fitted.prob,
    Point.probabilities = point.prob, n.best = k,
    class.levels = task$class.levels, class.values = values))
}
