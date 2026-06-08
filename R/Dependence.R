#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param asym logical; \code{FALSE} (default) Allows for asymmetrical dependencies.
#' @param p.value logical; \code{FALSE} (default) Generates 100 independent random permutations to test results against and plots 95 percent confidence intervals along with all results.
#' @param print.map logical; \code{FALSE} (default) Plots quadrant means, or p-value replicates.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#'
#' @note
#' For asymmetrical \code{(asym = TRUE)} matrices, directional dependence is returned as ([column variable] ---> [row variable]).
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
      original.par <- par(no.readonly = TRUE)
      on.exit(par(original.par), add = TRUE)
      
      nns.mc <- apply(x, 2L, function(g) NNS.dep(x[, 1L], g))
      cors   <- unlist(lapply(nns.mc, `[[`, "Correlation"))
      deps   <- unlist(lapply(nns.mc, `[[`, "Dependence"))
      
      cor_lower_CI <- LPM.VaR(.025, 0, cors[-c(1L, 2L)])
      cor_upper_CI <- UPM.VaR(.025, 0, cors[-c(1L, 2L)])
      dep_lower_CI <- LPM.VaR(.025, 0, deps[-c(1L, 2L)])
      dep_upper_CI <- UPM.VaR(.025, 0, deps[-c(1L, 2L)])
      
      if (print.map) {
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
        "Correlation p.value" = min(LPM(0, cors[2L], cors[-c(1L, 2L)]),
                                    UPM(0, cors[2L], cors[-c(1L, 2L)])),
        "Correlation 95% CIs" = c(cor_lower_CI, cor_upper_CI),
        "Dependence"          = as.numeric(deps[2L]),
        "Dependence p.value"  = min(LPM(0, deps[2L], deps[-c(1L, 2L)]),
                                    UPM(0, deps[2L], deps[-c(1L, 2L)])),
        "Dependence 95% CIs"  = c(dep_lower_CI, dep_upper_CI)
      ))
    }
    
    return(NNS.dep.matrix(x, asym = asym))
  }
  
  # ---- bivariate path -------------------------------------------------------
  if (length(x) != length(y)) stop("x and y must have the same length.")
  
  l   <- length(x)
  obs <- max(8L, as.integer(l / 8L))
  
  PART_xy <- suppressWarnings(
    NNS.part(x, y, order = NULL, obs.req = obs,
             min.obs.stop = FALSE, type = "XONLY", Voronoi = print.map)
  )
  PART_yx <- suppressWarnings(
    NNS.part(y, x, order = NULL, obs.req = obs,
             min.obs.stop = FALSE, type = "XONLY", Voronoi = FALSE)
  )
  
  if (nrow(PART_xy$regression.points) == 0L)
    return(list("Correlation" = 0, "Dependence" = 0))
  
  NNS_dep_pair_cpp(
    x       = as.numeric(x),
    y       = as.numeric(y),
    quad_xy = as.character(PART_xy$dt$quadrant),
    quad_yx = as.character(PART_yx$dt$quadrant),
    asym    = isTRUE(asym)
  )
}


NNS.dep.matrix <- function(x, order = NULL, degree = NULL, asym = FALSE){
  
  n <- ncol(x)
  if(is.null(n)){
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }
  
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
  
  x <- data.matrix(x)
  
  if(nrow(x) < 20 ) order <- 2
  
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
  
  return(list("Correlation" = rhos,
              "Dependence" = deps))
  
}
