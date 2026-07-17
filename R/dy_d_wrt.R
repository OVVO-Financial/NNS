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
# This restores the ORIGINAL NNS 0.5.7 finite-difference design used to produce
# the results in Vinod & Viole (2020), "Comparing Old and New Partial Derivative
# Estimates from Nonlinear Nonparametric Regressions" (SSRN 3681436).  It is a
# straight base-R port of the original data.table implementation (which no longer
# runs now that data.table has been dropped from Imports), and reproduces that
# version bit-for-bit when run against the same NNS.reg engine.  Key design:
#   * bandwidths h_s = 1/log(length(x), c(2, 10)); c(h_s, 10*h_s); doubled if
#     NNS.dep(x[, wrt], y) < 0.5;
#   * quantile-spacing step h_step = |mean(diff(LPM.VaR(seq(.01, 1, h), 0, x)))|;
#   * a degree-0 quantile grid LPM.VaR(seq(0, 1, .05), 0, .);
#   * NNS.reg(dim.red.method = "equal", point.only = TRUE) estimates (no smooth);
#   * plain mean aggregation, distance_wrt = 2 * h_step, and a plain rowMeans
#     (na.rm = TRUE) blend across bandwidths.
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

  l <- ncol(x)
  n <- nrow(x)

  results <- list()

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

  original.eval.points.min <- eval.points
  original.eval.points.max <- eval.points
  original.eval.points     <- eval.points

  # v0.5.7 bandwidths.
  h_s <- 1/log(length(x), c(2, 10))
  h_s <- c(h_s, 10 * h_s)
  if (NNS.dep(x[, wrt], y)$Dependence < .5) h_s <- 2 * h_s

  for (index in seq_along(h_s)) {
    h <- h_s[index]

    if (is.vector(eval.points) || dim(eval.points)[2] == 1) {
      eval.points <- unlist(eval.points)

      h_step <- abs(mean(diff(LPM.VaR(seq(.01, 1, h), 0, x[, wrt]))))

      if (is.finite(h_step)) {
        original.eval.points.min <- original.eval.points.min - h_step
        original.eval.points.max <- h_step + original.eval.points.max

        deriv.points <- apply(x, 2, function(z) LPM.VaR(seq(0, 1, .05), 0, z))
        if (dim(deriv.points)[2] != dim(x)[2]) {
          deriv.points <- matrix(deriv.points, ncol = l, byrow = FALSE)
        }
        sampsize <- length(seq(0, 1, .05))

        # base-R replacement for the data.table replicate + set().
        deriv.points <- do.call(rbind, replicate(3 * length(eval.points), deriv.points, simplify = FALSE))
        deriv.points[ , as.integer(wrt)] <- rep(as.vector(rbind(original.eval.points.min,
                                                                eval.points,
                                                                original.eval.points.max)),
                                                 each = sampsize, length.out = nrow(deriv.points))
        colnames(deriv.points) <- colnames(x)

        distance_wrt <- 2 * h_step

        position <- rep(rep(c("l", "m", "u"), each = sampsize), length.out = nrow(deriv.points))
        id <- rep(1:length(eval.points), each = 3 * sampsize, length.out = nrow(deriv.points))

        if (messages) {
          message(paste("Currently evaluating the ", nrow(deriv.points), " required points"),
                  "\r", appendLF = TRUE)
        }

        estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "equal",
                             plot = FALSE, threshold = 0, order = NULL, point.only = TRUE, ncores = 1)$Point.est
        estimates <- as.numeric(estimates)

        # base-R replacement for the data.table by = id plain-mean aggregation.
        ids <- 1:length(eval.points)
        lower   <- vapply(ids, function(g) mean(estimates[position == "l" & id == g]), numeric(1))
        two.f.x <- 2 * vapply(ids, function(g) mean(estimates[position == "m" & id == g]), numeric(1))
        upper   <- vapply(ids, function(g) mean(estimates[position == "u" & id == g]), numeric(1))
        rise <- upper - lower
      } else {
        kk <- length(unlist(eval.points))
        lower <- rep(NA_real_, kk); two.f.x <- rep(NA_real_, kk)
        upper <- rep(NA_real_, kk); rise <- rep(NA_real_, kk); distance_wrt <- NA_real_
      }

    } else {

      n <- dim(eval.points)[1]
      original.eval.points <- eval.points

      h_step <- abs(mean(diff(LPM.VaR(seq(.01, 1, h), 0, x[, wrt]))))

      if (is.finite(h_step)) {
        original.eval.points.min[ , wrt] <- original.eval.points.min[ , wrt] - h_step
        original.eval.points.max[ , wrt] <- h_step + original.eval.points.max[ , wrt]

        deriv.points <- rbind(original.eval.points.min,
                              original.eval.points,
                              original.eval.points.max)

        if (messages) {
          message("Currently generating NNS.reg finite difference estimates...bandwidth ",
                  index, " of ", length(h_s), "\r", appendLF = FALSE)
        }

        estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "equal",
                             plot = FALSE, threshold = 0, order = NULL, point.only = TRUE, ncores = 1)$Point.est
        estimates <- as.numeric(estimates)

        lower   <- head(estimates, n)
        two.f.x <- 2 * estimates[(n + 1):(2 * n)]
        upper   <- tail(estimates, n)
        rise <- upper - lower
        distance_wrt <- 2 * h_step
      } else {
        lower <- rep(NA_real_, n); two.f.x <- rep(NA_real_, n)
        upper <- rep(NA_real_, n); rise <- rep(NA_real_, n); distance_wrt <- NA_real_
      }
    }

    if (mixed) {
      if (is.null(dim(eval.points))) {
        if (length(eval.points) != 2) stop("Mixed Derivatives are only for 2 IV")
      } else {
        if (ncol(eval.points) != 2) stop("Mixed Derivatives are only for 2 IV")
      }

      if (!is.null(dim(eval.points))) {
        h_step_1 <- abs(mean(diff(LPM.VaR(seq(.01, 1, h), 0, x[ , 1]))))
        h_step_2 <- abs(mean(diff(LPM.VaR(seq(.01, 1, h), 0, x[ , 2]))))
        # Corner order per eval point: (+,+), (-,+), (+,-), (-,-); handles
        # multiple evaluation rows without interleaving (matches the Python port).
        ep <- as.matrix(eval.points)
        n_eval <- nrow(ep)
        mixed.deriv.points <- cbind(
          rep(ep[,1], each = 4L) + rep(c(h_step_1, -h_step_1, h_step_1, -h_step_1), times = n_eval),
          rep(ep[,2], each = 4L) + rep(c(h_step_2, h_step_2, -h_step_2, -h_step_2), times = n_eval)
        )
        colnames(mixed.deriv.points) <- colnames(ep)
        mixed.distances <- 4 * h_step_1 * h_step_2
      } else {
        h_step_1 <- h_step; h_step_2 <- h_step
        mixed.deriv.points <- matrix(c(h_step + eval.points,
                                       eval.points[1] - h_step, h_step + eval.points[2],
                                       h_step + eval.points[1], eval.points[2] - h_step,
                                       eval.points - h_step), ncol = 2, byrow = TRUE)
        mixed.distances <- (2 * h_step) * (2 * h_step)
      }

      if (is.finite(h_step_1) && is.finite(h_step_2)) {
        mixed.estimates <- NNS.reg(x, y, point.est = mixed.deriv.points, dim.red.method = "equal",
                                   plot = FALSE, threshold = 0, order = NULL, point.only = TRUE, ncores = 1)$Point.est
        mixed.estimates <- as.numeric(mixed.estimates)
        z <- matrix(mixed.estimates, ncol = 4, byrow = TRUE)
        z <- z[,1] + z[,4] - z[,2] - z[,3]
        mixed_band <- z / mixed.distances
      } else {
        mixed_band <- rep(NA_real_, if (is.null(dim(eval.points))) 1L else nrow(eval.points))
      }

      results[[index]] <- list("First"  = as.numeric(unlist(rise / distance_wrt)),
                               "Second" = as.numeric(unlist((upper - two.f.x + lower) / ((distance_wrt) ^ 2))),
                               "Mixed"  = mixed_band)
    } else {
      results[[index]] <- list("First"  = as.numeric(unlist(rise / distance_wrt)),
                               "Second" = as.numeric(unlist((upper - two.f.x + lower) / ((distance_wrt) ^ 2))))
    }
  }

  if (mixed) {
    final_results <- list("First"  = rowMeans(do.call(cbind, lapply(results, `[[`, 1)), na.rm = TRUE),
                          "Second" = rowMeans(do.call(cbind, lapply(results, `[[`, 2)), na.rm = TRUE),
                          "Mixed"  = rowMeans(do.call(cbind, lapply(results, `[[`, 3)), na.rm = TRUE))
  } else {
    final_results <- list("First"  = rowMeans(do.call(cbind, lapply(results, `[[`, 1)), na.rm = TRUE),
                          "Second" = rowMeans(do.call(cbind, lapply(results, `[[`, 2)), na.rm = TRUE))
  }

  if (messages) message("", "\r", appendLF = TRUE)
  return(final_results)
}

dy.d_ <- Vectorize(dy.d_, vectorize.args = c("wrt"))
