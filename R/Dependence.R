#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param asym logical; \code{FALSE} (default) returns the symmetric measures. If \code{TRUE}, returns the asymmetric (directional) measures: the partition is built on \code{x} only and the result measures the dependence of \code{y} on \code{x} (\code{y} as a function of \code{x}); swap the arguments for the other direction.
#' @param p.value logical; \code{FALSE} (default) If \code{TRUE} (bivariate only), compares the results against 100 random permutations of \code{y}: the Dependence p-value is the upper-tail permutation p-value \code{(1 + #(null >= observed)) / (1 + 100)}, the Correlation p-value is the two-sided permutation p-value \code{min(1, 2 * min(1 + #(null <= observed), 1 + #(null >= observed)) / (1 + 100))}, and the 95 percent intervals of the permutation null are returned (and plotted when \code{print.map = TRUE}).
#' @param print.map logical; \code{FALSE} (default) Plots quadrant means, or p-value replicates.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#'
#' @details
#' Dependence is a finite-resolution summary of the recursive partitions of \link{NNS.part}, which split each cell at the gravity of its observations (\code{noise.reduction = "off"}), not at its mean. At each cell, degree-0 (frequency) and degree-1 (area) concordance around the cell means are measured relative to a within-cell permutation reference and rescaled. Their cell-size-weighted averages are evaluated across partition levels. The reported Dependence is the strongest level, floored at zero; Correlation is a signed summary at that selected level and is bounded by Dependence in absolute value. With \code{asym = FALSE}, partitions on x, on y, and on joint quadrants are considered. With \code{asym = TRUE}, only partitions on x are considered, so the direction describes the resolution supplied by the conditioning variable rather than a causal effect.
#'
#' The split depth, stopping rules, finite occupancy, and maximization across levels affect the score. A zero score does not prove independence, a noiseless nonlinear function need not score one, and finite independent samples can produce positive scores. The permutation calculation repeats the same selection rule, providing a reference for the selected statistic rather than removing every finite-sample effect from the score itself.
#'
#' The degree-0 cell reference uses the hypergeometric distribution, with a normal approximation for cell sizes above 400. The degree-1 reference enumerates distinct arrangements only when at most 720 are required (all arrangements for up to six distinct values); larger cells use an approximation. These are properties of the cell calibration, separate from the 100 whole-sample permutations used by \code{p.value = TRUE}.
#'
#' Whole-sample shuffling assumes exchangeability of observations under the null. It is not an automatic test for serially dependent or clustered observations. The returned fields named \code{95\% CIs} are central intervals of the permutation null, not confidence intervals for a population dependence parameter. Inclusive-tail counts with the plus-one correction avoid zero Monte Carlo p-values; ties make the test conservative.
#'
#' @note
#' For asymmetrical \code{(asym = TRUE)} matrices, entry \code{[i, j]} is the directional dependence ([row variable i] ---> [column variable j]), i.e. \code{NNS.dep(x[, i], x[, j], asym = TRUE)}, the dependence of column \code{j} on row variable \code{i}.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.dep(x, y)
#'
#' ## Correlation / Dependence Matrix
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' B <- cbind(x, y, z)
#' NNS.dep(B)
#' }
#' @export

NNS.dep <- function(x,
                     y         = NULL,
                     asym      = FALSE,
                     p.value   = FALSE,
                     print.map = FALSE) {
  
  # ---- helper coercion ------------------------------------------------------
  .coerce_vec <- function(z, nm) {
    if (is.null(z)) return(NULL)
    
    if (any(class(z) %in% c("tbl", "data.table"))) {
      if (!is.null(ncol(z)) && ncol(z) == 1L) {
        z <- as.vector(unlist(z))
      } else {
        stop(sprintf("%s must be a vector or single-column object in the bivariate path.", nm))
      }
    }
    
    if (is.data.frame(z)) {
      if (ncol(z) == 1L) {
        z <- z[[1L]]
      } else {
        stop(sprintf("%s must be a vector or single-column object in the bivariate path.", nm))
      }
    }
    
    as.numeric(z)
  }
  
  # ---- class coercions ------------------------------------------------------
  if (!is.null(y)) {
    x <- .coerce_vec(x, "x")
    y <- .coerce_vec(y, "y")
  } else {
    if (any(class(x) %in% c("tbl", "data.table"))) x <- as.data.frame(x)
    if (is.data.frame(x)) x <- data.matrix(x)
  }
  
  # ---- missing values -------------------------------------------------------
  if (anyNA(x)) stop("x has missing values, please address.")
  if (!is.null(y) && anyNA(y)) stop("y has missing values, please address.")
  if (any(is.infinite(x))) stop("x must be finite.")
  if (!is.null(y) && any(is.infinite(y))) stop("y must be finite.")
  
  # ---- p.value permutation setup --------------------------------------------
  if (p.value) {
    if (is.null(y)) stop("p.value = TRUE requires both x and y.")
    if (length(x) != length(y)) stop("x and y must have the same length.")
    
    y_p <- replicate(100L, sample.int(length(y)))
    x   <- cbind(x, y, matrix(y[y_p], ncol = ncol(y_p), byrow = FALSE))
    y   <- NULL
  }
  
  # ---- matrix / p.value path ------------------------------------------------
  if (is.null(y)) {
    if (p.value) {
      nns.mc <- apply(x, 2L, function(g) NNS.dep(x[, 1L], g, asym = asym))
      cors   <- unlist(lapply(nns.mc, `[[`, "Correlation"))
      deps   <- unlist(lapply(nns.mc, `[[`, "Dependence"))
      
      cor_lower_CI <- LPM.VaR(.025, 0, cors[-c(1L, 2L)])
      cor_upper_CI <- UPM.VaR(.025, 0, cors[-c(1L, 2L)])
      dep_lower_CI <- LPM.VaR(.025, 0, deps[-c(1L, 2L)])
      dep_upper_CI <- UPM.VaR(.025, 0, deps[-c(1L, 2L)])
      
      if (print.map) {
        # par(no.readonly = TRUE) opens a device even when nothing is drawn.
        original.par <- par(no.readonly = TRUE)
        on.exit(par(original.par), add = TRUE)
        par(mfrow = c(1L, 2L))
        hist(cors[-c(1L, 2L)], main = "NNS Correlation", xlab = NULL,
             xlim = c(min(cors), max(cors[-1L])))
        abline(v = cors[2L], col = "red", lwd = 2)
        mtext("Result", side = 3L, col = "red", at = cors[2L])
        abline(v = cor_lower_CI, col = "red", lwd = 2, lty = 3)
        abline(v = cor_upper_CI, col = "red", lwd = 2, lty = 3)
        
        hist(deps[-c(1L, 2L)], main = "NNS Dependence", xlab = NULL,
             xlim = c(min(deps), max(deps[-1L])))
        abline(v = deps[2L], col = "red", lwd = 2)
        mtext("Result", side = 3L, col = "red", at = deps[2L])
        abline(v = dep_lower_CI, col = "red", lwd = 2, lty = 3)
        abline(v = dep_upper_CI, col = "red", lwd = 2, lty = 3)
      }
      
      return(list(
        "Correlation"         = as.numeric(cors[2L]),
        "Correlation p.value" = min(1, 2 * min(1 + sum(cors[-c(1L, 2L)] <= cors[2L]),
                                               1 + sum(cors[-c(1L, 2L)] >= cors[2L])) /
                                       (1 + length(cors[-c(1L, 2L)]))),
        "Correlation 95% CIs" = c(cor_lower_CI, cor_upper_CI),
        "Dependence"          = as.numeric(deps[2L]),
        # Dependence is a magnitude: only the upper tail is evidence against independence.
        "Dependence p.value"  = (1 + sum(deps[-c(1L, 2L)] >= deps[2L])) /
                                (1 + length(deps[-c(1L, 2L)])),
        "Dependence 95% CIs"  = c(dep_lower_CI, dep_upper_CI)
      ))
    }
    
    return(NNS.dep.matrix(x, asym = asym))
  }
  
  # ---- bivariate path -------------------------------------------------------
  if (length(x) != length(y)) stop("x and y must have the same length.")
  
  l   <- length(x)
  obs <- max(8L, as.integer(l / 8L))

  # Three partitions feed the measure: the recursive gravity split on x (y as a
  # function of x), on y, and the joint quadrants of (x, y). asym = TRUE only
  # needs the first. NNS_dep_pair_cpp reads every level of each partition from
  # the prefix-coded cell labels.
  #
  # Native fast path: call NNS_part_cpp in quadrants-only mode and skip the
  # data.table wrappers NNS.part() builds.  Quadrants are identical to the full
  # path, so the result is bit-identical.  Only when not plotting
  # (Voronoi/print.map needs the full partition).  Gated by
  # getOption("NNS.native").
  if (isTRUE(getOption("NNS.native", TRUE)) && !isTRUE(print.map)) {
    ord <- max(ceiling(log(l, 2)), 1L)             # NNS.part default order
    part <- function(a, b, type) {
      as.character(NNS_part_cpp(a, b, type = type, order_in = as.integer(ord),
                                obs_req = as.integer(obs), min_obs_stop = FALSE,
                                noise_reduction = "off", quadrants_only = TRUE)$quadrant)
    }
    qxy <- part(x, y, "XONLY")
    if (isTRUE(asym)) {
      qyx <- qj <- character(0)
    } else {
      qyx <- part(y, x, "XONLY")
      qj <- part(x, y, NULL)
    }
    return(NNS_dep_pair_cpp(x = as.numeric(x), y = as.numeric(y),
                            quad_xy = qxy, quad_yx = qyx, quad_joint = qj,
                            asym = isTRUE(asym)))
  }

  PART_xy <- suppressWarnings(
    NNS.part(x, y, order = NULL, obs.req = obs,
             min.obs.stop = FALSE, type = "XONLY", Voronoi = print.map)
  )
  
  if (nrow(PART_xy$regression.points) == 0L)
    return(list("Correlation" = 0, "Dependence" = 0))

  if (isTRUE(asym)) {
    qyx <- qj <- character(0)
  } else {
    qyx <- as.character(suppressWarnings(
      NNS.part(y, x, order = NULL, obs.req = obs,
               min.obs.stop = FALSE, type = "XONLY", Voronoi = FALSE)
    )$dt$quadrant)
    qj <- as.character(suppressWarnings(
      NNS.part(x, y, order = NULL, obs.req = obs,
               min.obs.stop = FALSE, type = NULL, Voronoi = FALSE)
    )$dt$quadrant)
  }
  
  NNS_dep_pair_cpp(
    x          = as.numeric(x),
    y          = as.numeric(y),
    quad_xy    = as.character(PART_xy$dt$quadrant),
    quad_yx    = qyx,
    quad_joint = qj,
    asym       = isTRUE(asym)
  )
}


# Local (finest-cell) dependence, NNS.dep's measure before the multi-level
# redesign (#95). NNS.reg's default order and dy.d_'s step rule read it as a
# fine-scale signal-to-noise measure (see NNS_dep_local_cpp); NNS.dep itself
# reports the multi-level measure.
.nns_dep_local <- function(x, y, asym = TRUE) {
  x <- as.numeric(x); y <- as.numeric(y)
  l <- length(x)
  obs <- max(8L, as.integer(l / 8L))
  ord <- max(ceiling(log(l, 2)), 1L)
  part <- function(a, b) {
    as.character(NNS_part_cpp(a, b, type = "XONLY", order_in = as.integer(ord),
                              obs_req = as.integer(obs), min_obs_stop = FALSE,
                              noise_reduction = "off", quadrants_only = TRUE)$quadrant)
  }
  qxy <- part(x, y)
  qyx <- if (isTRUE(asym)) character(0) else part(y, x)
  NNS_dep_local_cpp(x, y, qxy, qyx, isTRUE(asym))
}


NNS.dep.matrix <- function(x, asym = FALSE){
  
  n <- ncol(x)
  if(is.null(n)){
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }
  
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
  
  x <- data.matrix(x)
  
  upper_lower <- function(x, y, asym){
    basic_dep <- NNS.dep(x, y, print.map = FALSE, asym = asym)
    if(asym){
      asym_dep <- NNS.dep(y, x, print.map = FALSE, asym = asym)
      return(list("Upper_cor" = basic_dep$Correlation,
                  "Upper_dep" = basic_dep$Dependence,
                  "Lower_cor" = asym_dep$Correlation,
                  "Lower_dep" = asym_dep$Dependence))
    } else {
      return(list("Upper_cor" = basic_dep$Correlation,
                  "Upper_dep" = basic_dep$Dependence,
                  "Lower_cor" = basic_dep$Correlation,
                  "Lower_dep" = basic_dep$Dependence))
    }
  }
  
  raw.both <- lapply(1 : (n-1), function(i) sapply((i + 1) : n, function(b) upper_lower(x[ , i], x[ , b], asym = asym)))
  
  
  raw.both <- unlist(raw.both)
  l <- length(raw.both)
  
  raw.rhos_upper <- raw.both[seq(1, l, 4)]
  raw.deps_upper <- raw.both[seq(2, l, 4)]
  raw.rhos_lower <- raw.both[seq(3, l, 4)]
  raw.deps_lower <- raw.both[seq(4, l, 4)]
  
  rhos <- matrix(0, n, n)
  deps <- matrix(0, n, n)
  
  if(!asym){
    rhos[lower.tri(rhos, diag = FALSE)] <- (unlist(raw.rhos_upper) + unlist(raw.rhos_lower)) / 2
    deps[lower.tri(deps, diag = FALSE)] <- (unlist(raw.deps_upper) + unlist(raw.deps_lower)) / 2
    
    rhos[upper.tri(rhos)] <- t(rhos)[upper.tri(rhos)]
    deps[upper.tri(deps)] <- t(deps)[upper.tri(deps)]
  } else {
    rhos[lower.tri(rhos, diag = FALSE)] <- unlist(raw.rhos_lower)
    deps[lower.tri(deps, diag = FALSE)] <- unlist(raw.deps_lower)
    
    rhos_upper <- matrix(0, n, n)
    deps_upper <- matrix(0, n, n)
    
    rhos[is.na(rhos)] <- 0
    deps[is.na(deps)] <- 0
    
    rhos_upper[lower.tri(rhos_upper, diag=FALSE)] <- unlist(raw.rhos_upper)
    rhos_upper <- t(rhos_upper)
    
    deps_upper[lower.tri(deps_upper, diag=FALSE)] <- unlist(raw.deps_upper)
    deps_upper <- t(deps_upper)
    
    rhos <- rhos + rhos_upper
    deps <- deps + deps_upper
  }
  
  diag(rhos) <- 1
  diag(deps) <- 1
  
  colnames(rhos) <- colnames(x)
  colnames(deps) <- colnames(x)
  rownames(rhos) <- colnames(x)
  rownames(deps) <- colnames(x)
  
  return(.NNS.out(list("Correlation" = rhos,
              "Dependence" = deps)))

}
