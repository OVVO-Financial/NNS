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

NNS.part <- function(x, y, Voronoi = FALSE, type = NULL,
                     order = NULL, obs.req = 8, min.obs.stop = TRUE,
                     noise.reduction = "off") {
  noise.reduction <- tolower(noise.reduction)
  ok <- c("mean","median","mode","mode_class","off")
  if (!noise.reduction %in% ok)
    stop("noise.reduction must be one of ", paste(shQuote(ok), collapse = ", "))
  
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
  if(any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))
  
  if (is.null(obs.req)) obs.req <- 8L
  if (!is.null(order) && order == 0) order <- 1L
  
  n <- length(x)
  default.order <- max(ceiling(log(n, 2)), 1L)
  if (is.null(order)) order <- default.order
  
  out <- NNS_part_cpp(
    x = x, y = y,
    type = if (is.null(type)) NULL else as.character(type),
    order_in = as.integer(order),
    obs_req = as.integer(obs.req),
    min_obs_stop = isTRUE(min.obs.stop),
    noise_reduction = noise.reduction
  )
  
  PART <- data.table::as.data.table(out$dt)
  RP   <- data.table::as.data.table(out$`regression.points`)
  data.table::setorder(RP, quadrant)
  

  if (is.discrete(x)) RP[, x := ifelse(x %% 1 < 0.5, floor(x), ceiling(x))]
  
  if (isTRUE(Voronoi)) {
    mc <- match.call(); x.label <- deparse(mc$x); y.label <- deparse(mc$y)
    plot(x, y, col = "steelblue", cex.lab = 1.5, xlab = x.label, ylab = y.label)
    
    if (is.null(type)) {
      # draw dashed split segments (per-iteration, per-split group)
      sh <- out$segments_h
      if (NROW(sh)) segments(sh$x0, sh$y, sh$x1, sh$y, lty = 3)
      sv <- out$segments_v
      if (NROW(sv)) segments(sv$x,  sv$y0, sv$x,  sv$y1, lty = 3)
    } else {
      # XONLY: vertical ablines at group bounds each iteration
      vl <- out$vlines
      if (length(vl)) abline(v = vl, lty = 3)
    }
    
    points(RP$x, RP$y, pch = 15, lwd = 2, col = "red")
    title(main = paste0("NNS Order = ", out$order), cex.main = 2)
  }
  
  # Return the same shape as original
  list(order = as.integer(out$order), dt = PART[], regression.points = RP[])
}
