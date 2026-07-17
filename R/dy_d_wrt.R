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
  
  xstar <- function(mat) as.numeric(rowMeans(as.matrix(mat)))
  
  xstar.train <- xstar(x)
  xstar.train.design <- cbind(Xstar = xstar.train, Xstar2 = xstar.train)
  
  fd.estimates <- function(test_points) {
    xs <- xstar(test_points)
    as.numeric(
      NNS.stack(
        IVs.train = xstar.train.design,
        DV.train = y,
        IVs.test = cbind(Xstar = xs, Xstar2 = xs),
        method = 1,
        status = FALSE,
        order = NULL,
        folds = 1,
        ncores = 1
      )$stack
    )
  }
  
  if (length(wrt) != 1L || !is.numeric(wrt) ||
      wrt < 1L || wrt > l || wrt != as.integer(wrt)) {
    stop("`wrt` must select exactly one column of the expanded predictor matrix.")
  }
  wrt <- as.integer(wrt)
  
  if (l != 2L) mixed <- FALSE
  
  if (is.character(eval.points)) {
    if (length(eval.points) != 1L) {
      stop("Character `eval.points` must contain exactly one option.")
    }
    
    eval.points <- switch(
      tolower(eval.points),
      "median" = t(apply(x, 2, median)),
      "last"   = tail(x, 1),
      "mean"   = t(apply(x, 2, mean)),
      "apd"    = as.vector(x[, wrt, drop = FALSE]),
      "obs"    = x,
      stop("Unknown `eval.points` option.")
    )
  }
  
  
  base.eval.points <- eval.points
  
  norm.matrix <- apply(x, 2, function(z) NNS.rescale(z, 0, 1))
  
  zz <- max(
    NNS.dep(x[, wrt], y, asym = TRUE)$Dependence,
    NNS.copula(cbind(x[, wrt], x[, wrt], y)),
    NNS.copula(cbind(norm.matrix[, wrt], norm.matrix[, wrt], y)),
    na.rm = TRUE
  )
  
  root_n <- floor(sqrt(n))
  if (root_n < 2L) stop("Insufficient observations to construct finite-difference bandwidths.")
  
  # Keep the original five-bandwidth design, but prevent duplicate list indices.
  h_s <- round(exp(seq(log(2), log(root_n), length.out = 5)))
  
  results <- vector(mode = "list", length(h_s))
  
  base_h <- gravity(abs(diff(x[, wrt])))
  if (!is.finite(base_h) || base_h == 0) {
    rng <- abs(max(x[, wrt]) - min(x[, wrt]))
    if (!is.finite(rng) || rng == 0) {
      stop("Regressor `wrt` is constant; derivative is undefined.")
    }
    base_h <- rng / length(x[, wrt])
  }
  
  is_vector_branch <- is.null(dim(base.eval.points)) ||
    (!is.null(ncol(base.eval.points)) && ncol(base.eval.points) == 1L)
  
  if (is_vector_branch) {
    ep <- as.numeric(base.eval.points)
    
    seq_by <- max(0.01, (1 - zz) / 2)
    probs <- seq(0, 1, by = seq_by)
    if (tail(probs, 1) < 1) probs <- c(probs, 1)
    
    deriv.grid <- apply(x, 2, function(z) LPM.VaR(probs, 1, z))
    if (is.null(dim(deriv.grid)) || ncol(deriv.grid) != l) {
      deriv.grid <- matrix(deriv.grid, ncol = l, byrow = FALSE)
    }
    sampsize <- nrow(deriv.grid)
  } else {
    ep_matrix <- as.matrix(base.eval.points)
    if (ncol(ep_matrix) != l) {
      stop("Matrix/data-frame `eval.points` must have one column per expanded predictor.")
    }
  }
  
  for (index in seq_along(h_s)) {
    h_step <- base_h * h_s[index]
    
    if (!is.finite(h_step) || h_step <= 0) {
      stop("A non-positive finite-difference step was generated.")
    }
    
    if (is_vector_branch) {
      # Build each id as lower block, midpoint block, upper block.
      blocks <- lapply(seq_along(ep), function(g) {
        lower_grid <- deriv.grid
        middle_grid <- deriv.grid
        upper_grid <- deriv.grid
        
        lower_grid[, wrt] <- ep[g] - h_step
        middle_grid[, wrt] <- ep[g]
        upper_grid[, wrt] <- ep[g] + h_step
        
        rbind(lower_grid, middle_grid, upper_grid)
      })
      
      deriv.points <- do.call(rbind, blocks)
      colnames(deriv.points) <- colnames(x)
      
      id <- rep(seq_along(ep), each = 3L * sampsize)
      position <- rep(
        rep(c("l", "m", "u"), each = sampsize),
        times = length(ep)
      )
      
      if (messages) {
        message(
          "Currently evaluating the ", nrow(deriv.points),
          " required points ", index, " of ", length(h_s), "\r",
          appendLF = FALSE
        )
      }
      
      estimates <- fd.estimates(deriv.points)
      
      if (length(estimates) != nrow(deriv.points)) {
        stop("NNS.reg returned an unexpected number of point estimates.")
      }
      
      ids <- seq_along(ep)
      lower <- vapply(
        ids,
        function(g) gravity(estimates[position == "l" & id == g]),
        numeric(1)
      )
      f.x <- vapply(
        ids,
        function(g) gravity(estimates[position == "m" & id == g]),
        numeric(1)
      )
      upper <- vapply(
        ids,
        function(g) gravity(estimates[position == "u" & id == g]),
        numeric(1)
      )
      
      mixed_eval_points <- NULL
      
    } else {
      # CRITICAL FIX: create fresh lower and upper points for every bandwidth.
      current.min <- ep_matrix
      current.max <- ep_matrix
      current.min[, wrt] <- ep_matrix[, wrt] - h_step
      current.max[, wrt] <- ep_matrix[, wrt] + h_step
      
      deriv.points <- rbind(current.min, ep_matrix, current.max)
      n_eval <- nrow(ep_matrix)
      
      if (messages) {
        message(
          "Currently generating NNS.reg finite difference estimates...bandwidth ",
          index, " of ", length(h_s), "\r",
          appendLF = FALSE
        )
      }
      
      estimates <- fd.estimates(deriv.points)
      
      if (length(estimates) != 3L * n_eval) {
        stop("NNS.reg returned an unexpected number of point estimates.")
      }
      
      lower <- estimates[seq_len(n_eval)]
      f.x <- estimates[n_eval + seq_len(n_eval)]
      upper <- estimates[2L * n_eval + seq_len(n_eval)]
      
      mixed_eval_points <- ep_matrix
    }
    
    first_deriv <- (upper - lower) / (2 * h_step)
    second_deriv <- (upper - 2 * f.x + lower) / (h_step ^ 2)
    
    if (mixed) {
      if (is_vector_branch) {
        # Mixed derivatives require complete two-dimensional tuples.
        if (length(ep) != 2L) {
          stop("Mixed derivatives require a complete two-predictor evaluation tuple.")
        }
        mixed_eval_points <- matrix(ep, nrow = 1L)
        colnames(mixed_eval_points) <- colnames(x)
      }
      
      if (ncol(mixed_eval_points) != 2L) {
        stop("Mixed derivatives are only available for exactly two predictors.")
      }
      
      h_step_1 <- gravity(abs(diff(x[, 1]))) * h_s[index]
      h_step_2 <- gravity(abs(diff(x[, 2]))) * h_s[index]
      
      if (!is.finite(h_step_1) || h_step_1 == 0) {
        rng_1 <- abs(max(x[, 1]) - min(x[, 1]))
        h_step_1 <- (rng_1 / nrow(x)) * h_s[index]
      }
      if (!is.finite(h_step_2) || h_step_2 == 0) {
        rng_2 <- abs(max(x[, 2]) - min(x[, 2]))
        h_step_2 <- (rng_2 / nrow(x)) * h_s[index]
      }
      
      if (h_step_1 <= 0 || h_step_2 <= 0) {
        stop("Unable to construct valid mixed-derivative bandwidths.")
      }
      
      mixed.deriv.points <- do.call(
        rbind,
        lapply(seq_len(nrow(mixed_eval_points)), function(i) {
          p <- mixed_eval_points[i, ]
          rbind(
            c(p[1] + h_step_1, p[2] + h_step_2),
            c(p[1] - h_step_1, p[2] + h_step_2),
            c(p[1] + h_step_1, p[2] - h_step_2),
            c(p[1] - h_step_1, p[2] - h_step_2)
          )
        })
      )
      colnames(mixed.deriv.points) <- colnames(x)
      
      mixed.estimates <- fd.estimates(mixed.deriv.points)
      
      if (length(mixed.estimates) != 4L * nrow(mixed_eval_points)) {
        stop("NNS.reg returned an unexpected number of mixed-derivative estimates.")
      }
      
      z <- matrix(mixed.estimates, ncol = 4L, byrow = TRUE)
      mixed_deriv <- (z[, 1] - z[, 2] - z[, 3] + z[, 4]) /
        (4 * h_step_1 * h_step_2)
      
      results[[index]] <- list(
        First = first_deriv,
        Second = second_deriv,
        Mixed = mixed_deriv
      )
    } else {
      results[[index]] <- list(
        First = first_deriv,
        Second = second_deriv
      )
    }
  }
  
  weighted_mean <- function(values) {
    m <- do.call(cbind, values)
    if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
    weights <- rev(seq_len(ncol(m)))
    as.numeric(m %*% (weights / sum(weights)))
  }
  
  final_results <- list(
    First = weighted_mean(lapply(results, `[[`, "First")),
    Second = weighted_mean(lapply(results, `[[`, "Second"))
  )
  
  if (mixed) {
    final_results$Mixed <- weighted_mean(lapply(results, `[[`, "Mixed"))
  }
  
  if (messages) message("", "\r", appendLF = TRUE)
  final_results
}

dy.d_ <- Vectorize(dy.d_, vectorize.args = c("wrt"))