# Cross-validated NNS classification aggregates class probability vectors.
# Regression and its existing native search are kept in NNS.stack itself.
.nns_stack_classify <- function(x, y, z, class.values, splits, method, stack,
                                dim.red.method, order, dist, ncores,
                                optimize.threshold, pred.int, score, objective,
                                balance.indices, numeric.response, status,
                                numeric.design, tau = NULL) {
  values <- sort(unique(y))
  tie.order <- match(unique(y), values)
  labels <- as.character(class.values)
  n <- length(y); K <- length(values)
  encode <- function(train, test) {
    # Keep Stack's predictor names and factor columns in every fold. Scaling
    # remains in the NNS geometry and dimension-reduction routines below.
    design <- numeric.design(train, test, scale = FALSE)
    list(x = as.matrix(design$train), point.est = as.matrix(design$test),
         metadata = lapply(train, function(column) {
           list(categorical = is.factor(column) || is.logical(column))
         }))
  }
  full <- encode(x, z)
  p <- ncol(full$x)
  balanced <- lapply(splits, function(s) balance.indices(y[s$train]))
  final.rows <- balance.indices(y)
  decode <- function(code) {
    if (numeric.response) class.values[code] else as.numeric(code)
  }
  choose <- function(prob, threshold = .5) {
    if (K == 2L) {
      ifelse(prob[, tie.order[2L]] > threshold,
             values[tie.order[2L]], values[tie.order[1L]])
    } else .nns_class_select(prob, values, tie.order)
  }
  evaluate <- function(prob, actual) {
    thresholds <- if (K == 2L && optimize.threshold) seq(.01, .99, .01) else .5
    loss <- vapply(thresholds, function(t) score(decode(choose(prob, t)), decode(actual)), numeric(1))
    usable <- which(is.finite(loss))
    if (!length(usable)) stop("No classification threshold has a finite objective.")
    best <- if (objective == "min") min(loss[usable]) else max(loss[usable])
    tied <- usable[loss[usable] == best]
    take <- tied[which.min(abs(thresholds[tied] - .5))]
    list(score = loss[take], threshold = thresholds[take])
  }
  design <- function(train, response, test, component, select.count = NULL) {
    encoded <- encode(train, test)
    if (component == 1L && stack && !is.null(select.count) && select.count < ncol(encoded$x)) {
      ind <- outer(response, values, "==") * 1
      importance <- vapply(seq_len(ncol(encoded$x)), function(j) {
        .nns_class_strength(encoded$x[, j], ind)
      }, numeric(1))
      # Retain actual predictor coordinates, not two copies of the index.
      keep <- order(importance, decreasing = TRUE)[seq_len(select.count)]
      encoded$x <- encoded$x[, keep, drop = FALSE]
      encoded$point.est <- encoded$point.est[, keep, drop = FALSE]
    }
    encoded
  }
  fit <- function(train, response, test, component, candidate, select.count = NULL) {
    encoded <- design(train, response, test, component, select.count)
    task <- list(y = response, class.values = values, class.levels = labels)
    .nns_reg_classify(encoded, task, order,
                      if (component == 1L) candidate else NULL, dist, ncores,
                      if (component == 2L) dim.red.method else NULL, tau = tau,
                      dim.count = if (component == 2L) candidate else NULL,
                      point.only = TRUE)
  }
  # The method 1 cells, class shares and distances do not depend on the
  # neighbour count, so each fold builds them once for every candidate. The
  # probabilities equal separate fits at each count.
  fit.path <- function(train, response, test, candidates, select.count = NULL) {
    encoded <- design(train, response, test, 1L, select.count)
    categorical <- all(vapply(encoded$metadata, `[[`, logical(1), "categorical"))
    .nns_class_base(encoded$x, outer(response, values, "==") * 1, encoded$point.est,
                    order, NULL, dist, ncores, categorical,
                    need.fitted = FALSE, k.path = candidates)$predicted.path
  }
  record <- function(total, count, candidate) {
    valid <- count > 0
    probability <- matrix(NA_real_, n, K)
    probability[valid, ] <- total[valid, , drop = FALSE] / count[valid]
    assessment <- evaluate(probability[valid, , drop = FALSE], y[valid])
    c(assessment, list(candidate = candidate, probability = probability, valid = valid))
  }
  tune <- function(component, candidates, select.count = NULL) {
    if (component == 1L) {
      totals <- lapply(candidates, function(candidate) matrix(0, n, K))
      count <- integer(n)
      for (b in seq_along(splits)) {
        s <- splits[[b]]; tr <- s$train[balanced[[b]]]; va <- s$validation
        path <- fit.path(x[tr, , drop = FALSE], y[tr], x[va, , drop = FALSE],
                         candidates, select.count)
        for (a in seq_along(candidates))
          totals[[a]][va, ] <- totals[[a]][va, , drop = FALSE] + path[[a]]
        count[va] <- count[va] + 1L
      }
      records <- lapply(seq_along(candidates), function(a) {
        record(totals[[a]], count, candidates[[a]])
      })
      losses <- vapply(records, `[[`, numeric(1), "score")
      best <- if (objective == "min") which.min(losses) else which.max(losses)
      return(records[[best]])
    }
    records <- lapply(candidates, function(candidate) {
      total <- matrix(0, n, K); count <- integer(n)
      for (b in seq_along(splits)) {
        s <- splits[[b]]; tr <- s$train[balanced[[b]]]; va <- s$validation
        fitted <- fit(x[tr, , drop = FALSE], y[tr], x[va, , drop = FALSE],
                      component, candidate, select.count)
        total[va, ] <- total[va, , drop = FALSE] + as.matrix(fitted$Point.probabilities)
        count[va] <- count[va] + 1L
      }
      record(total, count, candidate)
    })
    losses <- vapply(records, `[[`, numeric(1), "score")
    best <- if (objective == "min") which.min(losses) else which.max(losses)
    records[[best]]
  }
  if (status) message("Selecting NNS class-frequency models on held-out probability vectors")
  dim.fit <- if (2L %in% method) tune(2L, as.list(seq_len(p))) else NULL
  reg.fit <- if (1L %in% method) {
    candidates <- if (p == 1L) list("all") else c(as.list(seq_len(max(1L, floor(sqrt(n))))), list("all"))
    tune(1L, candidates, if (is.null(dim.fit)) NULL else dim.fit$candidate)
  } else NULL
  weights <- c(reg = as.numeric(!is.null(reg.fit)), dim.red = as.numeric(!is.null(dim.fit)))
  threshold <- if (is.null(reg.fit)) dim.fit$threshold else reg.fit$threshold
  if (!is.null(reg.fit) && !is.null(dim.fit)) {
    # As in NNS 13.0, the models are weighted by their squared objective
    # scores, so nearly equal scores give a nearly equal blend. Choosing the
    # blend weight by its own out-of-fold accuracy gave one model all the
    # weight whenever it was ahead by a single out-of-fold row.
    weights <- .nns_stack_objective_weights(reg.fit$score, dim.fit$score, objective)
    valid <- reg.fit$valid & dim.fit$valid
    threshold <- evaluate(
      weights[["reg"]] * reg.fit$probability[valid, , drop = FALSE] +
        weights[["dim.red"]] * dim.fit$probability[valid, , drop = FALSE],
      y[valid])$threshold
  }
  # The coefficient cutoff at which NNS.reg(x, y, dim.red.method = ..., threshold
  # = cutoff) keeps the tuned number of predictors in every class index,
  # computed as for regression on all training rows: the smallest kept
  # |coefficient|, or 0 when every predictor is kept. Ties at the cutoff keep
  # every tied predictor. Two-class indices are mirror images, so one number
  # serves both; with more classes each index can need its own cutoff, and the
  # cutoffs are returned per class, which NNS.reg also accepts.
  dim.cutoff <- function(count) {
    cutoffs <- vapply(seq_along(values), function(j) {
      coef <- .nns_reg_dimred_coefficients(full$x, as.numeric(y == values[j]),
                                           dim.red.method, tau, 0)
      if (count >= length(coef)) return(0)
      cut <- sort(abs(coef), decreasing = TRUE)[count]
      if (is.finite(cut)) cut else 0
    }, numeric(1))
    if (isTRUE(all.equal(min(cutoffs), max(cutoffs)))) min(cutoffs) else
      setNames(cutoffs, labels)
  }
  final <- function(component, record) {
    if (is.null(record)) return(NULL)
    fit(x[final.rows, , drop = FALSE], y[final.rows], z, component,
        record$candidate, if (is.null(dim.fit)) NULL else dim.fit$candidate)
  }
  reg.final <- final(1L, reg.fit); dim.final <- final(2L, dim.fit)
  rp <- if (is.null(reg.final)) NULL else as.matrix(reg.final$Point.probabilities)
  dp <- if (is.null(dim.final)) NULL else as.matrix(dim.final$Point.probabilities)
  probability <- if (is.null(rp)) dp else if (is.null(dp)) rp else weights[1] * rp + weights[2] * dp
  colnames(probability) <- labels
  prediction.set <- function(prob) {
    if (is.null(prob) || is.null(pred.int)) return(NULL)
    lapply(seq_len(nrow(prob)), function(i) {
      o <- tie.order[order(prob[i, tie.order], decreasing = TRUE)]
      count <- which(cumsum(prob[i, o]) >= pred.int - 1e-12)[1L]
      decode(values[o[seq_len(count)]])
    })
  }
  interval <- function(prob) {
    sets <- prediction.set(prob)
    if (is.null(sets) || K != 2L) return(NULL)
    setNames(as.data.frame(t(vapply(sets, function(s) range(as.numeric(s)), numeric(2)))),
             c("pred.int.neg", "pred.int.pos"))
  }
  result <- list(
    OBJfn.reg = if (is.null(reg.fit)) NA_real_ else reg.fit$score,
    # As for regression, a win for all cells reports the final fit's cell count.
    NNS.reg.n.best = if (is.null(reg.fit)) NA_integer_ else if (identical(reg.fit$candidate, "all"))
      as.integer(reg.final$n.best) else reg.fit$candidate,
    probability.threshold = threshold,
    OBJfn.dim.red = if (is.null(dim.fit)) NA_real_ else dim.fit$score,
    NNS.dim.red.threshold = if (is.null(dim.fit)) NA_real_ else dim.cutoff(dim.fit$candidate),
    NNS.dim.red.count = if (is.null(dim.fit)) NA_integer_ else dim.fit$candidate,
    NNS.dim.red.variables = if (is.null(dim.final)) NULL else unique(unlist(lapply(
      dim.final$equation, function(eq) eq$Variable[eq$Variable != "DENOMINATOR" & abs(eq$Coefficient) > 0]))),
    NNS.dim.red.coefficients = if (is.null(dim.final)) NULL else dim.final$equation,
    reg = if (is.null(rp)) NA_real_ else decode(choose(rp, reg.fit$threshold)),
    reg.pred.int = interval(rp), dim.red = if (is.null(dp)) NA_real_ else decode(choose(dp, dim.fit$threshold)),
    dim.red.pred.int = interval(dp), stack = decode(choose(probability, threshold)),
    pred.int = interval(probability), Prediction.set = prediction.set(probability),
    probabilities = probability, reg.probabilities = rp, dim.red.probabilities = dp,
    weights = weights, class.levels = class.values)
  .NNS.out(result)
}

# NNS 13.0's stack weights: each model's squared objective score, or its inverse
# square when the objective is minimized, normalized to sum to one. A zero score
# counts as 1e-10, a non-finite weight as 0, and two zero weights give an equal
# blend.
.nns_stack_objective_weights <- function(reg.score, dim.score, objective) {
  s <- c(reg.score, dim.score)
  s[!is.na(s) & s == 0] <- 1e-10
  w <- pmax(if (identical(objective, "min")) 1 / s^2 else s^2, 1e-10)
  w[!is.finite(w)] <- 0
  w <- if (sum(w) > 0) w / sum(w) else c(.5, .5)
  c(reg = w[1], dim.red = w[2])
}
