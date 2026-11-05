# NNS.VAR completes a variable whose latest values are missing with the fit NNS
# 13.0 obtained from NNS.stack(cbind(index, index), y, method = 1) on the
# observed periods.
#
# That NNS.stack call always selected one neighbour (n.best = 1): its
# cross-validation scored larger neighbour counts with a k-NN path that read the
# predictor matrix in row-major order (since corrected), so they never won. The
# result was therefore NNS.reg's multivariate regression with one neighbour on
# the duplicated index. That regression has since changed in ways that matter
# here: its last cell absorbs the newest observation, estimates beyond the
# observed range are no longer extended along the boundary gradient, and its
# default order partitions an index more coarsely. The functions below
# reproduce the 13.0 computation. Run with NNS 13.0's own
# NNS.dep, NNS.reg and gravity, they return 13.0's estimates exactly; with the
# current versions of those functions the regression points can differ
# slightly.

# Regression points on the index: NNS.reg's univariate partition at the order
# NNS 13.0 used for order = NULL. That order is ten times the mean of
# NNS.dep(asym = TRUE) and the copula of (x, x, y), rounded down, halved and
# rounded down again below 100 observations, and at least 1. (13.0 wrote the
# rounding as dependence*10 %% 1 < .5, which R reads as dependence*(10 %% 1),
# so it always rounded down.) At dependence 1 every observed index is a point.
.nns_var_index_points <- function(x, y) {
  # 13.0 computed this inside NNS.reg, which silenced warnings (for example
  # NNS.rescale's on a constant series).
  dep <- tryCatch(suppressWarnings(NNS.dep(x, y, print.map = FALSE, asym = TRUE)$Dependence),
                  error = function(e) .1)
  dep <- tryCatch(suppressWarnings(mean(c(dep, NNS.copula(cbind(apply(cbind(x, x, y), 2,
                                                      function(z) NNS.rescale(z, 0, 1))))))),
                  error = function(e) dep)
  dep[is.na(dep)] <- .1
  ord <- floor(dep * 10)
  if (length(y) < 100) ord <- floor(ord / 2)
  ord <- max(1, ord)
  if (dep == 1) ord <- "max"
  # type = "REG": 13.0 treated every response as numeric here, including a
  # 0/1 series that NNS.reg would now read as classes.
  rp <- NNS.reg(x, y, order = ord, type = "REG", plot = FALSE, ncores = 1,
                multivariate.call = TRUE)$x
  sort(unique(as.numeric(rp)))
}

# Estimates at x.new from observed values y at index positions x. `points`
# supplies the regression points directly (tests only).
.nns_var_index_fit <- function(x, y, x.new, points = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x.new <- as.numeric(x.new)
  if (length(y) < 2L)
    return(rep(if (length(y)) y[length(y)] else NA_real_, length(x.new)))

  rp <- if (is.null(points)) .nns_var_index_points(x, y) else
    sort(unique(as.numeric(points)))

  # Cells are left-closed intervals between regression points. The newest
  # observation sits on the last regression point, so it forms its own cell
  # and is retained as a regression point. Each cell is summarised by the
  # gravity of its index and of its values; cells are ordered by their
  # identifiers as 13.0 sorted them, which decides ties.
  cell <- findInterval(x, rp)
  ids <- paste(cell, cell, sep = ".")
  grp <- factor(ids, levels = sort(unique(ids), method = "radix"))
  rx <- vapply(split(x, grp), gravity, numeric(1), USE.NAMES = FALSE)
  ry <- vapply(split(y, grp), gravity, numeric(1), USE.NAMES = FALSE)

  # One neighbour under 13.0's distance on the two identical index columns,
  # scaled by the range of the regression points and the query point.
  nearest <- function(q) {
    rng <- max(rx, q) - min(rx, q)
    if (is.finite(rng) && rng > 0) {
      d <- (rx - q) * (1 / rng)
      s <- 2 * (d * d + abs(d))
      s[s == 0] <- 1e-10
    } else {
      s <- rep(1e-10, length(rx))
    }
    ry[which.min(s)]
  }

  est <- vapply(x.new, nearest, numeric(1))

  # Beyond the observed index: start from the estimate at the boundary and
  # extend it along the gradient toward the centre of the regression points,
  # averaging the gradients to the centre and to two points between the
  # boundary and the centre with weights 3:2:1. Distances are measured on the
  # two index columns, as in 13.0.
  lo <- min(x)
  hi <- max(x)
  outside <- which(x.new < lo | x.new > hi)
  if (length(outside)) {
    centre <- gravity(rx)
    e.centre <- nearest(centre)
    for (i in outside) {
      b <- min(max(x.new[i], lo), hi)
      m1 <- (b + centre) / 2
      m2 <- (b + m1) / 2
      e.b <- nearest(b)
      g1 <- (e.b - e.centre) / max(sqrt((b - centre)^2 + (b - centre)^2), 1e-10)
      g2 <- (e.b - nearest(m1)) / max(sqrt((b - m1)^2 + (b - m1)^2), 1e-10)
      g3 <- (e.b - nearest(m2)) / max(sqrt((b - m2)^2 + (b - m2)^2), 1e-10)
      gradient <- (g1 * 3 + g2 * 2 + g3 * 1) / 6
      est[i] <- sqrt((x.new[i] - b)^2 + (x.new[i] - b)^2) * gradient + e.b
    }
  }
  est
}
