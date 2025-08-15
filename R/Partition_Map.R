#' NNS Partition Map
#'
#' Creates partitions based on partial moment quadrant centroids, iteratively assigning identifications to observations based on those quadrants (unsupervised partitional and hierarchial clustering method).  Basis for correlation, dependence \link{NNS.dep}, regression \link{NNS.reg} routines.
#'
#' @param x a numeric vector.
#' @param y a numeric vector with compatible dimensions to \code{x}.
#' @param Voronoi logical; \code{FALSE} (default) Displays a Voronoi type diagram using partial moment quadrants.
#' @param type \code{NULL} (default) Controls the partitioning basis.  Set to \code{(type = "XONLY")} for X-axis based partitioning.  Defaults to \code{NULL} for both X and Y-axis partitioning.
#' @param order integer; Number of partial moment quadrants to be generated.  \code{(order = "max")} will institute a perfect fit.
#' @param obs.req integer; (8 default) Required observations per cluster where quadrants will not be further partitioned if observations are not greater than the entered value.  Reduces minimum number of necessary observations in a quadrant to 1 when \code{(obs.req = 1)}.
#' @param min.obs.stop logical; \code{TRUE} (default) Stopping condition where quadrants will not be further partitioned if a single cluster contains less than the entered value of \code{obs.req}.
#' @param noise.reduction the method of determining regression points options for the dependent variable \code{y}: ("mean", "median", "mode", "off"); \code{(noise.reduction = "mean")} uses means for partitions.  \code{(noise.reduction = "median")} uses medians instead of means for partitions, while \code{(noise.reduction = "mode")} uses modes instead of means for partitions.  Defaults to \code{(noise.reduction = "off")} where an overall central tendency measure is used, which is the default for the independent variable \code{x}.
#' @return Returns:
#'  \itemize{
#'   \item{\code{"dt"}} a \code{data.table} of \code{x} and \code{y} observations with their partition assignment \code{"quadrant"} in the 3rd column and their prior partition assignment \code{"prior.quadrant"} in the 4th column.
#'   \item{\code{"regression.points"}} the \code{data.table} of regression points for that given \code{(order = ...)}.
#'   \item{\code{"order"}}  the \code{order} of the final partition given \code{"min.obs.stop"} stopping condition.
#'   }
#'
#' @note \code{min.obs.stop = FALSE} will not generate regression points due to unequal partitioning of quadrants from individual cluster observations.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.part(x, y)
#'
#' ## Data.table of observations and partitions
#' NNS.part(x, y, order = 1)$dt
#'
#' ## Regression points
#' NNS.part(x, y, order = 1)$regression.points
#'
#' ## Voronoi style plot
#' NNS.part(x, y, Voronoi = TRUE)
#'
#' ## Examine final counts by quadrant
#' DT <- NNS.part(x, y)$dt
#' DT[ , counts := .N, by = quadrant]
#' DT
#' }
#' @export

NNS.part = function(x, y, Voronoi = FALSE, type = NULL,
                    order = NULL, obs.req = 8, min.obs.stop = TRUE,
                    noise.reduction = "off") {
  noise.reduction <- tolower(noise.reduction)
  if (!noise.reduction %in% c("mean","median","mode","mode_class","off"))
    stop("noise.reduction must be one of 'mean','median','mode','mode_class','off'")
  
  if (any(class(x) %in% c("tbl","data.table"))) x <- unlist(x, use.names = FALSE)
  if (any(class(y) %in% c("tbl","data.table"))) y <- unlist(y, use.names = FALSE)
  x <- as.numeric(x); y <- as.numeric(y)
  if (is.null(obs.req)) obs.req <- 8
  if (!is.null(order) && order == 0) order <- 1
  
  # --- centralized inline reducers (no function-object switching) ---
  nr_x <- function(v) {
    if (noise.reduction == "mean")        mean(v)
    else if (noise.reduction == "median") median(v)
    else if (noise.reduction == "mode")   NNS::NNS.mode(v, discrete = TRUE, multi = FALSE)
    else if (noise.reduction == "mode_class") NNS::NNS.gravity(v)  # gravity for x
    else                                  NNS::NNS.gravity(v)      # "off"
  }
  nr_y <- function(v) {
    if (noise.reduction == "mean")        mean(v)
    else if (noise.reduction == "median") median(v)
    else if (noise.reduction == "mode")   NNS::NNS.mode(v, discrete = TRUE, multi = FALSE)
    else if (noise.reduction == "mode_class") NNS::NNS.mode(v, discrete = TRUE, multi = FALSE)
    else                                  NNS::NNS.gravity(v)      # "off"
  }
  # -----------------------------------------------------------------
  
  if (Voronoi) {
    mc <- match.call(); x.label <- deparse(mc$x); y.label <- deparse(mc$y)
    plot(x, y, col = "steelblue", cex.lab = 1.5, xlab = x.label, ylab = y.label)
  }
  
  PART <- data.table::data.table(x = x, y = y, quadrant = "q", prior.quadrant = "pq")
  PART[, counts := .N, by = quadrant]
  PART[, old.counts := .N, by = prior.quadrant]
  
  n <- length(x)
  default.order <- max(ceiling(log(n, 2)), 1)
  if (is.null(order)) order <- default.order
  OR <- obs.req
  
  drawSegments <- function() {
    if (is.null(type)) {
      PART[split.rows, {
        yh <- nr_y(y); xv <- nr_x(x)
        segments(min(x), yh, max(x), yh, lty = 3)
        segments(xv, min(y), xv, max(y), lty = 3)
      }, by = quadrant]
    } else {
      bounds <- PART[, .(min = min(x), max = max(x)), by = quadrant]
      abline(v = bounds$min, lty = 3); abline(v = bounds$max, lty = 3)
    }
  }
  
  obs_assignment <- function() {
    RP[, prior.quadrant := quadrant]
    PART[split.rows, prior.quadrant := quadrant]
    if (is.null(type)) {
      PART[RP, on = .(quadrant), `:=`(q_new = { lox <- x.x <= i.x; loy <- x.y <= i.y; 1L + lox + loy * 2L })]
    } else {
      PART[RP, on = .(quadrant), `:=`(q_new = { lox <- x.x > i.x; 1L + lox })]
    }
    PART[split.rows, quadrant := paste0(quadrant, q_new)]
  }
  
  i <- 0L
  while (TRUE) {
    if (nrow(PART) > n || i >= order || i >= floor(log(n, 2))) break
    PART[, counts := .N, by = quadrant]
    split.rows <- PART[counts > OR, which = TRUE]
    if (length(split.rows) == 0) break
    
    if (Voronoi) drawSegments()
    
    RP <- PART[split.rows, .( x = nr_x(x), y = nr_y(y) ), by = quadrant]
    
    obs_assignment()
    i <- i + 1L
    if (min.obs.stop) {
      PART[, counts := .N, by = quadrant]
      if (min(PART$counts) <= OR) break
    }
  }
  
  PART[, c("counts","old.counts","q_new") := NULL]
  
  RP <- PART[, .( x = nr_x(x), y = nr_y(y) ), by = prior.quadrant]
  data.table::setnames(RP, "prior.quadrant", "quadrant")
  
  if (is.discrete(x)) {
    RP[, x := ifelse(x %% 1 < 0.5, floor(x), ceiling(x))]
  }
  RP <- data.table::setorder(RP, quadrant)
  
  if (Voronoi) {
    title(main = paste0("NNS Order = ", i), cex.main = 2)
    points(RP$x, RP$y, pch = 15, lwd = 2, col = "red")
  }
  
  list(order = i, dt = PART[], regression.points = RP)
}
