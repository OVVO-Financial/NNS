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

  # v0.5.7 bandwidths.
  h_s <- 1/log(length(x), c(2, 10))
  h_s <- c(h_s, 10 * h_s)
  if (NNS.dep(x[, wrt], y)$Dependence < .5) h_s <- 2 * h_s

  col <- x[, wrt]

  band_first <- list()
  band_second <- list()
  band_mixed <- list()

  is_vector <- is.vector(eval.points) || (!is.null(ncol(eval.points)) && ncol(eval.points) == 1)

  if (is_vector) {
    eval_vec <- as.numeric(unlist(eval.points))
    k <- length(eval_vec)
    grid <- apply(x, 2, function(z) LPM.VaR(seq(0, 1, .05), 0, z))
    if (is.null(dim(grid)) || ncol(grid) != l) grid <- matrix(grid, ncol = l, byrow = FALSE)
    sampsize <- nrow(grid)

    for (bi in seq_along(h_s)) {
      H <- h_s[bi]
      steps <- vapply(eval_vec, function(ev) dydx.step(col, ev, H), numeric(1))

      blocks <- vector("list", 3L * k)
      position <- character(0)
      id <- integer(0)
      for (g in seq_len(k)) {
        lower_g <- grid; middle_g <- grid; upper_g <- grid
        lower_g[, wrt]  <- eval_vec[g] - steps[g]
        middle_g[, wrt] <- eval_vec[g]
        upper_g[, wrt]  <- eval_vec[g] + steps[g]
        blocks[[3L * g - 2L]] <- lower_g
        blocks[[3L * g - 1L]] <- middle_g
        blocks[[3L * g]]      <- upper_g
        position <- c(position, rep(c("l", "m", "u"), each = sampsize))
        id <- c(id, rep(g, 3L * sampsize))
      }
      deriv.points <- do.call(rbind, blocks)
      colnames(deriv.points) <- colnames(x)
      estimates <- stack.est(x, y, deriv.points)

      f1 <- numeric(k); f2 <- numeric(k)
      for (g in seq_len(k)) {
        lo <- mean(estimates[position == "l" & id == g])
        mm <- mean(estimates[position == "m" & id == g])
        up <- mean(estimates[position == "u" & id == g])
        h <- steps[g]
        if (is.finite(h) && h != 0) {
          f1[g] <- (up - lo) / (2 * h)
          f2[g] <- (up - 2 * mm + lo) / (h ^ 2)
        } else {
          f1[g] <- NA_real_; f2[g] <- NA_real_
        }
      }
      band_first[[bi]] <- f1
      band_second[[bi]] <- f2
    }
    mixed_eval <- if (k == 2L) matrix(eval_vec, nrow = 1L) else NULL

  } else {
    eval_mat <- as.matrix(eval.points)
    if (ncol(eval_mat) != l) stop("Matrix/data-frame `eval.points` must have one column per expanded predictor.")
    n_eval <- nrow(eval_mat)

    for (bi in seq_along(h_s)) {
      H <- h_s[bi]
      steps <- vapply(seq_len(n_eval), function(i) dydx.step(col, eval_mat[i, wrt], H), numeric(1))
      finite <- is.finite(steps) & steps != 0

      lower <- eval_mat; upper <- eval_mat
      lower[, wrt] <- eval_mat[, wrt] - steps
      upper[, wrt] <- eval_mat[, wrt] + steps
      deriv.points <- rbind(lower, eval_mat, upper)
      colnames(deriv.points) <- colnames(x)
      estimates <- stack.est(x, y, deriv.points)

      lo <- estimates[seq_len(n_eval)]
      mm <- estimates[n_eval + seq_len(n_eval)]
      up <- estimates[2L * n_eval + seq_len(n_eval)]
      f1 <- (up - lo) / (2 * steps)
      f2 <- (up - 2 * mm + lo) / (steps ^ 2)
      f1[!finite] <- NA_real_; f2[!finite] <- NA_real_
      band_first[[bi]] <- f1
      band_second[[bi]] <- f2
    }
    mixed_eval <- eval_mat
  }

  row_nanmean <- function(bands) {
    m <- do.call(cbind, bands)
    if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
    rowMeans(m, na.rm = TRUE)
  }

  if (mixed) {
    if (is.null(mixed_eval) || ncol(mixed_eval) != 2) stop("Mixed Derivatives are only for 2 IV")
    for (bi in seq_along(h_s)) {
      H <- h_s[bi]
      vals <- vapply(seq_len(nrow(mixed_eval)), function(i) {
        p <- mixed_eval[i, ]
        s1 <- dydx.step(x[, 1], p[1], H)
        s2 <- dydx.step(x[, 2], p[2], H)
        if (!is.finite(s1) || !is.finite(s2) || s1 == 0 || s2 == 0) return(NA_real_)
        corners <- rbind(c(p[1] + s1, p[2] + s2),
                         c(p[1] - s1, p[2] + s2),
                         c(p[1] + s1, p[2] - s2),
                         c(p[1] - s1, p[2] - s2))
        colnames(corners) <- colnames(x)
        z <- stack.est(x, y, corners)
        (z[1] + z[4] - z[2] - z[3]) / (4 * s1 * s2)
      }, numeric(1))
      band_mixed[[bi]] <- vals
    }
  }

  final_results <- list("First"  = row_nanmean(band_first),
                        "Second" = row_nanmean(band_second))
  if (mixed) final_results$Mixed <- row_nanmean(band_mixed)

  if (messages) message("", "\r", appendLF = TRUE)
  return(final_results)
}

dy.d_ <- Vectorize(dy.d_, vectorize.args = c("wrt"))
