#' Finite-difference weights for a local polynomial stencil
#'
#' Construct polynomial finite-difference weights on distinct local nodes.
#' Invalid or degenerate stencils return NA coefficients.
#'
#' @param nodes Numeric vector of distinct stencil coordinates.
#' @param at Numeric evaluation coordinate.
#' @param degree Integer derivative order; defaults to 1.
#' @return Numeric vector of finite-difference weights, or NA coefficients
#'   when the stencil is invalid.
#' @keywords internal
#' @export

.nns_fd_weights <- function(nodes, at, degree = 1L) {
  if (length(nodes) <= degree || length(unique(nodes)) != length(nodes))
    return(rep(NA_real_, length(nodes)))
  scale <- max(abs(nodes - at))
  if (!is.finite(scale) || scale == 0) return(rep(NA_real_, length(nodes)))
  z <- (nodes - at) / scale
  powers <- outer(seq.int(0, length(nodes) - 1L), z, function(k, value) value^k)
  rhs <- numeric(length(nodes)); rhs[degree + 1L] <- factorial(degree)
  as.numeric(solve(powers, rhs)) / scale^degree
}

.nns_dy_stencil <- function(col, at, probability.width) {
  support <- sort(unique(col))
  lower <- support[1L]; upper <- tail(support, 1L)
  if (!is.finite(at) || at < lower || at > upper)
    stop("Derivative evaluation points must lie inside the observed coordinate support.")
  if (length(support) < 2L) return(list(nodes = at, first = NA_real_, second = NA_real_))
  if (length(support) == 2L) {
    nodes <- support
  } else {
    share <- LPM.ratio(1, at, col)
    lo <- LPM.VaR(max(0, share - probability.width), 1, col)
    hi <- LPM.VaR(min(1, share + probability.width), 1, col)
    # A local half-window, never an entire-support extrapolation step.
    h <- min((hi - lo) / 2, (upper - lower) / 4)
    if (!is.finite(h) || h <= 0) h <- min(diff(support))
    if (at == lower) {
      nodes <- at + c(0, h, 2 * h)
    } else if (at == upper) {
      nodes <- at - c(2 * h, h, 0)
    } else {
      nodes <- c(max(lower, at - h), at, min(upper, at + h))
    }
  }
  list(nodes = nodes, first = .nns_fd_weights(nodes, at, 1L),
       second = .nns_fd_weights(nodes, at, 2L))
}

#' Partial Derivative dy/d_[wrt]
#'
#' Returns the numerical partial derivative of \code{y} with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{NNS.reg} estimates as \code{f(x + h)} and \code{f(x - h)} values.
#'
#' @param x a numeric matrix or data frame.
#' @param y a numeric vector with compatible dimensions to \code{x}.
#' @param wrt integer; Selects the regressor to differentiate with respect to (vectorized).
#' @param eval.points numeric or options: ("obs", "apd", "mean", "median", "last"); Regressor points to be evaluated.
#' \itemize{
#' \item Numeric values must be in matrix or data.frame form to be evaluated for each regressor, otherwise, a vector of points will evaluate only at the \code{wrt} regressor, averaging over the observed rows of the other regressors (at most \code{max(21, 10000 / k)} of them for \code{k} points, taken evenly through the data).  See examples for use cases.
#' \item Set to \code{(eval.points = "obs")} (default) to find the average partial derivative at every observation of the variable with respect to \emph{for specific tuples of given observations.}
#' \item Set to \code{(eval.points = "apd")} to find the average partial derivative at every observation of the variable with respect to \emph{over the entire distribution of other regressors}, represented by their observed rows.
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
#' @details
#' Derivatives are finite-difference summaries of the fitted \link{NNS.reg} surface, averaged over four local probability windows derived from degree-1 partial-moment quantiles. Window sizes depend on the number of observations; they are capped before constructing a local coordinate half-width. All evaluation points and perturbed coordinates must stay inside the observed coordinate ranges. Interior unequal spacings use polynomial finite-difference weights; support endpoints use one-sided three-point formulas. This avoids averaging artificial full-support extrapolation bands into a local derivative.
#'
#' The coordinate-range check does not establish that every perturbed tuple lies in a densely observed part of the joint support. Sparse cells, neighbour changes, fitted-surface bias, and correlated predictors can affect the result, including at the boundary. These are numerical derivatives of an estimated surface, with no universal recovery or coverage guarantee. A two-valued coordinate returns its finite contrast divided by the coordinate difference and an undefined (NA) second derivative. Constant coordinates return NA. Mixed derivatives require two predictors and explicit matrix evaluation points.
#'
#' @note Binary regressors provide a finite contrast, not an identified continuous derivative between the two levels.
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
#' dy.d_(B, y, wrt = 1, eval.points = t(c(.5, median(x_2))))
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' ## only supply 1 value in [eval.points], or a vector of [eval.points]:
#' dy.d_(B, y, wrt = 1, eval.points = .5)
#'
#' dy.d_(B, y, wrt = 1, eval.points = fivenum(B[,1]))
#'
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' ## for every observation of 1st regressor:
#' apd <- dy.d_(B, y, wrt = 1, eval.points = "apd")
#' plot(B[,1], apd[,1]$First)
#'
#' ## Central 95% range across estimated partial derivatives (descriptive)
#' ### Lower CI
#' LPM.VaR(.025, 0, apd[,1]$First)
#'
#' ### Upper CI
#' UPM.VaR(.025, 0, apd[,1]$First)
#' }
#' @export

dy.d_ <- function(x, y, wrt, eval.points = "obs", mixed = FALSE, messages = TRUE) {
  if (is.null(ncol(x)) || ncol(x) < 2L)
    stop("Please supply at least two predictor columns; use dy.dx for one predictor.")
  if (anyNA(cbind(x, y))) stop("You have some missing values, please address.")
  blocks <- lapply(seq_len(ncol(x)), function(j) factor_2_dummy_FR(x[, j]))
  x <- as.matrix(do.call(cbind, blocks))
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  y <- as.numeric(y)
  if (length(y) != nrow(x) || any(!is.finite(c(x, y))))
    stop("x and y must have compatible finite observations.")
  if (length(wrt) != 1L || !is.finite(wrt) || wrt != floor(wrt) || wrt < 1L || wrt > ncol(x))
    stop("wrt must identify one expanded predictor column.")
  if (is.character(eval.points)) {
    eval.points <- switch(tolower(eval.points),
      median = matrix(apply(x, 2, median), nrow = 1),
      mean = matrix(colMeans(x), nrow = 1), last = tail(x, 1),
      apd = x[, wrt], obs = x,
      stop("Unknown eval.points option."))
  }
  marginal <- is.null(dim(eval.points)) || ncol(eval.points) == 1L
  if (marginal) {
    # A one-column matrix or data frame holds values of the wrt regressor.
    points <- as.numeric(unlist(eval.points, use.names = FALSE))
    count <- min(nrow(x), max(21L, floor(10000 / max(1L, length(points)))))
    other <- x[.nns_spread_rows(nrow(x), count), , drop = FALSE]
    contexts <- lapply(points, function(value) {
      z <- other; z[, wrt] <- value; z
    })
  } else {
    ev <- as.matrix(eval.points)
    if (ncol(ev) != ncol(x)) stop("eval.points must have one column per expanded predictor.")
    if (any(!is.finite(ev)) || any(sweep(ev, 2, apply(x, 2, min), "<")) ||
        any(sweep(ev, 2, apply(x, 2, max), ">")))
      stop("Derivative evaluation points must lie inside the observed coordinate support.")
    contexts <- lapply(seq_len(nrow(ev)), function(i) ev[i, , drop = FALSE])
    points <- ev[, wrt]
  }
  if (!length(points)) stop("eval.points must contain at least one point.")
  # Four genuinely local probability windows, based on observations rather
  # than the matrix's element count. The former tenfold windows saturated.
  base <- 1 / log(max(nrow(x), 3), c(2, 10))
  widths <- pmin(.25, c(base / 2, base))
  if (.nns_dep_local(x[, wrt], y, asym = FALSE) < .5) widths <- pmin(.25, 2 * widths)
  queries <- list(); specifications <- list()
  for (b in seq_along(widths)) for (i in seq_along(contexts)) {
    context <- contexts[[i]]
    stencil <- .nns_dy_stencil(x[, wrt], points[i], widths[b])
    nodes <- lapply(stencil$nodes, function(value) {
      z <- context; z[, wrt] <- value; z
    })
    index <- length(queries) + 1L
    queries[[index]] <- do.call(rbind, nodes)
    specifications[[index]] <- list(band = b, point = i, n = nrow(context), stencil = stencil)
  }
  mixed.spec <- list()
  if (mixed && ncol(x) == 2L && !marginal) {
    for (b in seq_along(widths)) for (i in seq_along(contexts)) {
      at <- contexts[[i]][1, ]
      s1 <- .nns_dy_stencil(x[, 1], at[1], widths[b])
      s2 <- .nns_dy_stencil(x[, 2], at[2], widths[b])
      grid <- as.matrix(expand.grid(s1$nodes, s2$nodes))
      colnames(grid) <- colnames(x)
      index <- length(queries) + 1L
      queries[[index]] <- grid
      mixed.spec[[length(mixed.spec) + 1L]] <- list(index = index, band = b, point = i,
        weight = as.numeric(outer(s1$first, s2$first)))
    }
  } else if (mixed) {
    stop("Mixed derivatives require two predictors and explicit matrix evaluation points.")
  }
  if (messages) message("Generating local NNS.reg finite-difference estimates...Regressor ", wrt)
  lengths <- vapply(queries, nrow, integer(1))
  query <- do.call(rbind, queries); colnames(query) <- colnames(x)
  predicted <- as.numeric(NNS.reg(x, y, type = "REG", point.est = query, plot = FALSE, ncores = 1)$Point.est)
  offsets <- c(0L, cumsum(lengths))
  part <- function(index) predicted[seq.int(offsets[index] + 1L, offsets[index + 1L])]
  first <- second <- matrix(NA_real_, length(points), length(widths))
  for (index in seq_along(specifications)) {
    spec <- specifications[[index]]
    values <- colMeans(matrix(part(index), nrow = spec$n))
    values <- values - values[1L]
    first[spec$point, spec$band] <- sum(spec$stencil$first * values)
    second[spec$point, spec$band] <- sum(spec$stencil$second * values)
  }
  average <- function(z) {
    out <- rowMeans(z, na.rm = TRUE); out[is.nan(out)] <- NA_real_; out
  }
  result <- list(First = average(first), Second = average(second))
  if (mixed) {
    values <- matrix(NA_real_, length(points), length(widths))
    for (spec in mixed.spec) {
      f <- part(spec$index)
      values[spec$point, spec$band] <- sum(spec$weight * (f - f[1L]))
    }
    result$Mixed <- average(values)
  }
  result
}

dy.d_ <- Vectorize(dy.d_, vectorize.args = c("wrt"))


# k of n row indices spread through the data without a fixed stride: the rows
# whose golden-ratio (Weyl) sequence values are smallest. A fixed stride can
# line up with periodic or sorted row order and keep only some phases; an
# irrational rotation visits every residue class evenly.
.nns_spread_rows <- function(n, k) {
  if (k >= n) return(seq_len(n))
  sort(order((seq_len(n) * 0.6180339887498949) %% 1)[seq_len(k)])
}
