#' Partial Derivative dy/d_[wrt]
#'
#' Returns the numerical partial derivative of \code{y} with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{NNS.reg} estimates as \code{f(x + h)} and \code{f(x - h)} values.
#'
#' @param x a numeric matrix or data frame.
#' @param y a numeric vector with compatible dimensions to \code{x}.
#' @param wrt integer; Selects the regressor to differentiate with respect to (vectorized).
#' @param eval.points numeric or options: ("obs", "apd", "mean", "median", "last"); Regressor points to be evaluated.
#' \itemize{
#' \item Numeric values must be in matrix or data.frame form to be evaluated for each regressor, otherwise, a vector of points will evaluate only at the \code{wrt} regressor.  See examples for use cases.
#' \item Set to \code{(eval.points = "obs")} (default) to find the average partial derivative at every observation of the variable with respect to \emph{for specific tuples of given observations.}
#' \item Set to \code{(eval.points = "apd")} to find the average partial derivative at every observation of the variable with respect to \emph{over the entire distribution of other regressors.}
#' \item Set to \code{(eval.points = "mean")} to find the partial derivative at the mean of value of every variable.
#' \item Set to \code{(eval.points = "median")} to find the partial derivative at the median value of every variable.
#' \item Set to \code{(eval.points = "last")} to find the partial derivative at the last observation of every value (relevant for time-series data).
#' }
#' @param mixed logical; \code{FALSE} (default) If mixed derivative is to be evaluated, set \code{(mixed = TRUE)}.
#' @param messages logical; \code{TRUE} (default) Prints status messages.
#' @return Returns column-wise matrix of wrt regressors:
#' \itemize{
#' \item{\code{dy.d_(...)[, wrt]$First}} the 1st derivative
#' \item{\code{dy.d_(...)[, wrt]$Second}} the 2nd derivative
#' \item{\code{dy.d_(...)[, wrt]$Mixed}} the mixed derivative (for two independent variables only).
#' }
#'
#'
#' @note For binary regressors, it is suggested to use \code{eval.points = seq(0, 1, .05)} for a better resolution around the midpoint.
#'
#' @author Fred Viole, OVVO Financial Systems
#'
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#'
#' Vinod, H. and Viole, F. (2020) "Comparing Old and New Partial Derivative Estimates from Nonlinear Nonparametric Regressions"  \doi{10.2139/ssrn.3681104}
#'
#' @examples
#' \dontrun{
#' set.seed(123) ; x_1 <- runif(1000) ; x_2 <- runif(1000) ; y <- x_1 ^ 2 * x_2 ^ 2
#' B <- cbind(x_1, x_2)
#'
#' ## To find derivatives of y wrt 1st regressor for specific points of both regressors
#' dy.d_(B, y, wrt = 1, eval.points = t(c(.5, 1)))
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' only supply 1 value in [eval.points], or a vector of [eval.points]:
#' dy.d_(B, y, wrt = 1, eval.points = .5)
#'
#' dy.d_(B, y, wrt = 1, eval.points = fivenum(B[,1]))
#'
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' for every observation of 1st regressor:
#' apd <- dy.d_(B, y, wrt = 1, eval.points = "apd")
#' plot(B[,1], apd[,1]$First)
#'
#' ## 95% Confidence Interval to test if 0 is within
#' ### Lower CI
#' LPM.VaR(.025, 0, apd[,1]$First)
#'
#' ### Upper CI
#' UPM.VaR(.025, 0, apd[,1]$First)
#' }
#' @export


# -----------------------------------------------------------------------------
# Reconciled dy.d_ : the original NNS 0.5.7 finite-difference design
# (Vinod & Viole 2020, SSRN 3681436) reconciled with the current NNS.reg engine.
# Two changes make the estimates uniform across identically-distributed
# regressors and independent of the retired data.table machinery:
#   * h_step shares the dy.dx() logic - a locally-adaptive step centred on the
#     evaluation point's percentile:
#         p      <- LPM.ratio(1, eval, x[, wrt])
#         h_step <- LPM.VaR(p + H, 1, x[, wrt]) - LPM.VaR(p - H, 1, x[, wrt])
#     (no cumulative window). This removes the cross-regressor scatter.
#   * estimates come from NNS.stack on the equal-weight synthetic regressor X*
#     via the increased-dimension trick cbind(X*, X*), with method = c(1, 2),
#     dim.red.method = "equal", order = "max", folds = 5.  The cross-validated
#     n.best regularises the (sharper) current engine back toward the paper's
#     regime.
# Bandwidths follow v0.5.7: h_s = 1/log(length(x), c(2, 10)); c(h_s, 10*h_s);
# doubled when NNS.dep(x[, wrt], y) < 0.5.  First = (upper - lower)/(2*h_step);
# Second = (upper - 2*f(x) + lower)/h_step^2 (matching dy.dx); rowMeans(na.rm)
# blend across bandwidths.
# -----------------------------------------------------------------------------

dy.d_ <- function(x, y, wrt,
                  eval.points = "obs",
                  mixed = FALSE,
                  messages = TRUE){

  n <- nrow(x)
  l <- ncol(x)

  if (is.null(l)) stop("Please ensure (x) is a matrix or data.frame type object.")
  if (l < 2) stop("Please use NNS::dy.dx(...) for univariate partial derivatives.")
  if (anyNA(cbind(x, y))) stop("You have some missing values, please address.")

  dummies <- list()
  for (i in seq_len(l)) {
    dummies[[i]] <- factor_2_dummy_FR(x[, i])
    if (!is.null(ncol(dummies[[i]]))) {
      base_name <- if (is.null(colnames(x))) paste0("X", i) else colnames(x)[i]
      colnames(dummies[[i]]) <- paste0(base_name, "_", colnames(dummies[[i]]))
    }
  }
  x <- do.call(cbind, dummies)

  if (messages) {
    message("Currently generating NNS.reg finite difference estimates...Regressor ",
            wrt, "\r", appendLF = TRUE)
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", seq_len(ncol(x)))
  }

  if (any(class(x) %in% c("tbl", "data.table"))) x <- as.data.frame(x)
  if (!is.null(y) && any(class(y) %in% c("tbl", "data.table"))) {
    y <- as.vector(unlist(y))
  }

  x <- as.matrix(x)
  l <- ncol(x)
  n <- nrow(x)

  if (l != 2) mixed <- FALSE

  if (is.character(eval.points)) {
    eval.points <- tolower(eval.points)
    eval.points <- switch(eval.points,
                          "median" = t(apply(x, 2, median)),
                          "last"   = tail(x, 1),
                          "mean"   = t(apply(x, 2, mean)),
                          "apd"    = as.vector(x[ , wrt, drop = FALSE]),
                          x)
  }

  # ---- Estimates: NNS.stack on the equal-weight synthetic regressor X* -------
  stack.est <- function(train, resp, pts) {
    xs <- as.numeric(rowMeans(as.matrix(train)))
    ps <- as.numeric(rowMeans(as.matrix(pts)))
    as.numeric(NNS.stack(
      IVs.train = cbind(Xstar = xs, Xstar2 = xs),
      DV.train = resp,
      IVs.test = cbind(Xstar = ps, Xstar2 = ps),
      method = c(1, 2),
      dim.red.method = "equal",
      order = "max",
      folds = 5,
      status = FALSE,
      ncores = 1
    )$stack)
  }

  # ---- dy.dx() locally-adaptive step -----------------------------------------
  dydx.step <- function(col, ev, H) {
    p <- LPM.ratio(1, ev, col)
    LPM.VaR(min(1, p + H), 1, col) - LPM.VaR(max(0, p - H), 1, col)
  }

  col <- x[, wrt]

  # v0.5.7 bandwidths.
  h_s <- 1/log(length(x), c(2, 10))
  h_s <- c(h_s, 10 * h_s)
  if (NNS.dep(col, y)$Dependence < .5) h_s <- 2 * h_s

  is_vector <- is.vector(eval.points) || (!is.null(ncol(eval.points)) && ncol(eval.points) == 1)

  # The NNS.stack fit (CV n.best, dimension-reduction coefficients, blend weight)
  # depends only on (x, y), never on the evaluation points, so every bandwidth -
  # and the mixed-derivative corners - are predicted from a single fit. Gather
  # every test block, run ONE NNS.stack, then slice it back.
  chunks <- list(); csize <- integer(0)
  push <- function(block) {
    i <- length(chunks) + 1L
    chunks[[i]] <<- block
    csize[i] <<- nrow(block)
    i
  }

  steps_per_band <- vector("list", length(h_s))
  main_idx <- integer(length(h_s))
  position <- NULL; id <- NULL

  if (is_vector) {
    eval_vec <- as.numeric(unlist(eval.points))
    k <- length(eval_vec)
    grid <- apply(x, 2, function(z) LPM.VaR(seq(0, 1, .05), 0, z))
    if (is.null(dim(grid)) || ncol(grid) != l) grid <- matrix(grid, ncol = l, byrow = FALSE)
    sampsize <- nrow(grid)
    position <- rep(rep(c("l", "m", "u"), each = sampsize), times = k)
    id <- rep(seq_len(k), each = 3L * sampsize)

    for (bi in seq_along(h_s)) {
      H <- h_s[bi]
      steps <- vapply(eval_vec, function(ev) dydx.step(col, ev, H), numeric(1))
      blocks <- vector("list", 3L * k)
      for (g in seq_len(k)) {
        lower_g <- grid; middle_g <- grid; upper_g <- grid
        lower_g[, wrt]  <- eval_vec[g] - steps[g]
        middle_g[, wrt] <- eval_vec[g]
        upper_g[, wrt]  <- eval_vec[g] + steps[g]
        blocks[[3L * g - 2L]] <- lower_g
        blocks[[3L * g - 1L]] <- middle_g
        blocks[[3L * g]]      <- upper_g
      }
      steps_per_band[[bi]] <- steps
      main_idx[bi] <- push(do.call(rbind, blocks))
    }
    mixed_eval <- if (k == 2L) matrix(eval_vec, nrow = 1L) else NULL

  } else {
    eval_mat <- as.matrix(eval.points)
    if (ncol(eval_mat) != l) stop("Matrix/data-frame `eval.points` must have one column per expanded predictor.")
    n_eval <- nrow(eval_mat)

    for (bi in seq_along(h_s)) {
      H <- h_s[bi]
      steps <- vapply(seq_len(n_eval), function(i) dydx.step(col, eval_mat[i, wrt], H), numeric(1))
      lower <- eval_mat; upper <- eval_mat
      lower[, wrt] <- eval_mat[, wrt] - steps
      upper[, wrt] <- eval_mat[, wrt] + steps
      steps_per_band[[bi]] <- steps
      main_idx[bi] <- push(rbind(lower, eval_mat, upper))
    }
    mixed_eval <- eval_mat
  }

  # Mixed-derivative corners (also predicted from the same single fit).
  mixed_meta <- vector("list", length(h_s))
  if (mixed) {
    if (is.null(mixed_eval) || ncol(mixed_eval) != 2) stop("Mixed Derivatives are only for 2 IV")
    for (bi in seq_along(h_s)) {
      H <- h_s[bi]
      corner_blocks <- list(); scales <- numeric(nrow(mixed_eval))
      for (m in seq_len(nrow(mixed_eval))) {
        p <- mixed_eval[m, ]
        s1 <- dydx.step(x[, 1], p[1], H); s2 <- dydx.step(x[, 2], p[2], H)
        if (is.finite(s1) && is.finite(s2) && s1 != 0 && s2 != 0) {
          corner_blocks[[length(corner_blocks) + 1L]] <- rbind(
            c(p[1] + s1, p[2] + s2), c(p[1] - s1, p[2] + s2),
            c(p[1] + s1, p[2] - s2), c(p[1] - s1, p[2] - s2))
          scales[m] <- 4 * s1 * s2
        } else scales[m] <- NA_real_
      }
      if (length(corner_blocks)) {
        mixed_meta[[bi]] <- list(idx = push(do.call(rbind, corner_blocks)), scales = scales)
      } else {
        mixed_meta[[bi]] <- list(idx = NA_integer_, scales = scales)
      }
    }
  }

  # ---- the single NNS.stack fit + prediction --------------------------------
  big <- do.call(rbind, chunks)
  colnames(big) <- colnames(x)
  preds <- stack.est(x, y, big)
  offs <- c(0L, cumsum(csize))
  parts <- lapply(seq_along(csize), function(i) preds[(offs[i] + 1L):offs[i + 1L]])

  band_first <- vector("list", length(h_s))
  band_second <- vector("list", length(h_s))
  for (bi in seq_along(h_s)) {
    steps <- steps_per_band[[bi]]; block <- parts[[main_idx[bi]]]
    if (is_vector) {
      kk <- length(steps); f1 <- numeric(kk); f2 <- numeric(kk)
      for (g in seq_len(kk)) {
        lo <- mean(block[position == "l" & id == g])
        mm <- mean(block[position == "m" & id == g])
        up <- mean(block[position == "u" & id == g])
        h <- steps[g]
        if (is.finite(h) && h != 0) {
          f1[g] <- (up - lo) / (2 * h)
          f2[g] <- (up - 2 * mm + lo) / (h ^ 2)
        } else {
          f1[g] <- NA_real_; f2[g] <- NA_real_
        }
      }
    } else {
      n_eval <- length(steps)
      lo <- block[seq_len(n_eval)]
      mm <- block[n_eval + seq_len(n_eval)]
      up <- block[2L * n_eval + seq_len(n_eval)]
      finite <- is.finite(steps) & steps != 0
      f1 <- (up - lo) / (2 * steps)
      f2 <- (up - 2 * mm + lo) / (steps ^ 2)
      f1[!finite] <- NA_real_; f2[!finite] <- NA_real_
    }
    band_first[[bi]] <- f1
    band_second[[bi]] <- f2
  }

  row_nanmean <- function(bands) {
    m <- do.call(cbind, bands)
    if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
    rowMeans(m, na.rm = TRUE)
  }

  final_results <- list("First"  = row_nanmean(band_first),
                        "Second" = row_nanmean(band_second))

  if (mixed) {
    band_mixed <- vector("list", length(h_s))
    for (bi in seq_along(h_s)) {
      meta <- mixed_meta[[bi]]
      vals <- rep(NA_real_, length(meta$scales))
      if (!is.na(meta$idx)) {
        z <- parts[[meta$idx]]; pos <- 0L
        for (m in seq_along(meta$scales)) {
          if (is.finite(meta$scales[m])) {
            c4 <- z[(pos + 1L):(pos + 4L)]
            vals[m] <- (c4[1] + c4[4] - c4[2] - c4[3]) / meta$scales[m]
            pos <- pos + 4L
          }
        }
      }
      band_mixed[[bi]] <- vals
    }
    final_results$Mixed <- row_nanmean(band_mixed)
  }

  if (messages) message("", "\r", appendLF = TRUE)
  return(final_results)
}

dy.d_ <- Vectorize(dy.d_, vectorize.args = c("wrt"))
