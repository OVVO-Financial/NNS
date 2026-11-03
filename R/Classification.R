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
                            dist, ncores, categorical = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  strength <- vapply(seq_len(p), function(j) {
    .nns_class_strength(X[, j], probabilities)
  }, numeric(1))
  boundaries <- lapply(seq_len(p), function(j) {
    if (identical(order, "max") || categorical) return(sort(unique(X[, j])))
    depth <- if (is.null(order)) {
      d <- max(1L, as.integer(floor(strength[j] * 10 + .5)))
      if (nrow(X) < 100L) d <- max(1L, as.integer(floor(d / 2)))
      d
    } else order
    .nns_reg_partition_points_fast(X[, j], X[, j], order = depth,
                                   noise.reduction = "mean", is.class = FALSE)
  })
  id.parts <- lapply(seq_len(p), function(j) {
    if (identical(order, "max") || categorical) match(X[, j], boundaries[[j]]) else
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
      max(1L, floor((1 - mean(strength)) * sqrt(nrow(X))))
  k <- max(1L, min(as.integer(k), nrow(centres)))
  predict <- function(points) {
    if (is.null(points)) return(NULL)
    points <- as.matrix(points)
    if (p == 1L && !categorical) {
      o <- order(centres[, 1])
      ans <- vapply(seq_len(ncol(shares)), function(j) {
        if (length(o) == 1L) rep(shares[o, j], nrow(points)) else
          stats::approx(centres[o, 1], shares[o, j], xout = points[, 1],
                        rule = 2, ties = mean)$y
      }, numeric(nrow(points)))
    } else {
      # Equal-distance cells share their rank weight, including a tie at k.
      # Exact matches retain their empirical class frequencies.
      ans <- t(vapply(seq_len(nrow(points)), function(i) {
        d <- .nns_mreg_distance_vector(centres, points[i, ], dist,
                                       apply(X, 2, min), apply(X, 2, max))
        exact <- which(d == 0)
        if (length(exact)) return(colMeans(shares[exact, , drop = FALSE]))
        cutoff <- sort(d, partial = k)[k]
        selected <- which(d <= cutoff)
        selected <- selected[order(d[selected])]
        weights <- .nns_mreg_ensemble_weights(d[selected])
        for (v in unique(d[selected])) {
          tied <- which(d[selected] == v)
          weights[tied] <- mean(weights[tied])
        }
        as.numeric(crossprod(weights, shares[selected, , drop = FALSE]))
      }, numeric(ncol(shares))))
    }
    ans <- matrix(ans, ncol = ncol(shares))
    colnames(ans) <- colnames(shares)
    .nns_class_normalize(ans, prior)
  }
  list(fitted = predict(X), predicted = predict(point.est),
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
                             order, n.best, dist, ncores)
      list(fitted = fit$fitted[, 2], predicted = if (is.null(Q)) NULL else fit$predicted[, 2],
           equation = reduced$equation)
    })
    fitted.prob <- .nns_class_normalize(do.call(cbind, lapply(models, `[[`, "fitted")), colMeans(indicators))
    point.prob <- if (is.null(Q)) NULL else
      .nns_class_normalize(do.call(cbind, lapply(models, `[[`, "predicted")), colMeans(indicators))
    equations <- setNames(lapply(models, `[[`, "equation"), labels)
    centres <- shares <- boundaries <- NULL
    ids <- seq_len(nrow(X))
    k <- NA_integer_
  } else {
    fit <- .nns_class_base(X, indicators, Q, order, n.best, dist, ncores, categorical)
    fitted.prob <- fit$fitted; point.prob <- fit$predicted
    centres <- fit$centres; shares <- fit$shares
    boundaries <- fit$boundaries; ids <- fit$ids; k <- fit$k
  }
  colnames(fitted.prob) <- labels
  if (!is.null(point.prob)) colnames(point.prob) <- labels
  fitted.pred <- .nns_class_select(fitted.prob, values, tie.order)
  point.pred <- .nns_class_select(point.prob, values, tie.order)
  fitted <- as.data.frame(X)
  fitted$y <- task$y; fitted$y.hat <- fitted.pred
  fitted$NNS.ID <- ids; fitted$error <- as.numeric(fitted.pred != task$y)
  accuracy <- mean(fitted.pred == task$y)
  if (plot && ncol(X) == 1L) {
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
