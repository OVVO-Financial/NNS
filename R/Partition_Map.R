#' NNS Partition Map
#'
#' Creates partitions based on partial-moment quadrant centroids. Numeric orders
#' use the compiled recursive partitioner. `order = "max"` returns the maximum
#' observable partition representation without passing an invalid integer to C++.
#'
#' @param x Numeric vector.
#' @param y Numeric vector of the same length as x.
#' @param Voronoi Logical; draw the partition map.
#' @param type NULL or "XONLY".
#' @param order NULL, a positive integer, or "max".
#' @param obs.req Nonnegative integer minimum-observation stopping control.
#' @param min.obs.stop Logical stopping control.
#' @param noise.reduction One of "mean", "median", "mode", "mode_class", or "off".
#'
#' @return A list containing `order`, `dt`, and `regression.points`.
#' @export
#' 
#' 
NNS.part <- function(x, y, Voronoi = FALSE, type = NULL,
                     order = NULL, obs.req = 8, min.obs.stop = TRUE,
                     noise.reduction = "off") {
  
  # Capture the original calls before x and y are coerced or subsetted.
  # Otherwise substitute(x) / substitute(y) deparse the entire numeric vectors
  # and place those values on the plot axes.
  x.label <- paste(deparse(substitute(x)), collapse = " ")
  y.label <- paste(deparse(substitute(y)), collapse = " ")
  
  if (inherits(x, c("tbl", "data.table")) || is.data.frame(x)) {
    if (NCOL(x) != 1L) {
      stop("[x] must be a vector or one-column object.", call. = FALSE)
    }
    x <- x[[1L]]
  }
  if (inherits(y, c("tbl", "data.table")) || is.data.frame(y)) {
    if (NCOL(y) != 1L) {
      stop("[y] must be a vector or one-column object.", call. = FALSE)
    }
    y <- y[[1L]]
  }
  if (is.matrix(x)) {
    if (ncol(x) != 1L) {
      stop("[x] must be a vector or one-column object.", call. = FALSE)
    }
    x <- x[, 1L]
  }
  if (is.matrix(y)) {
    if (ncol(y) != 1L) {
      stop("[y] must be a vector or one-column object.", call. = FALSE)
    }
    y <- y[, 1L]
  }
  
  if (length(x) != length(y)) {
    stop("[x] and [y] must have the same length.", call. = FALSE)
  }
  if (length(x) < 1L) {
    stop("[x] and [y] must not be empty.", call. = FALSE)
  }
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Complete-case handling: drop pairs with NA/NaN in either variable but
  # keep infinities, matching the historical NNS.part contract.
  complete <- !(is.na(x) | is.na(y))
  if (!all(complete)) {
    x <- x[complete]
    y <- y[complete]
  }
  if (length(x) < 1L) {
    stop(
      "[x] and [y] must contain at least one complete (non-missing) pair.",
      call. = FALSE
    )
  }
  
  Voronoi <- .nns_reg_scalar_logical(Voronoi, "Voronoi")
  min.obs.stop <- .nns_reg_scalar_logical(min.obs.stop, "min.obs.stop")
  
  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || is.na(type) ||
        tolower(type) != "xonly") {
      stop("[type] must be NULL or 'XONLY'.", call. = FALSE)
    }
    type <- "XONLY"
  }
  
  if (!is.character(noise.reduction) || length(noise.reduction) != 1L ||
      is.na(noise.reduction)) {
    stop("Invalid [noise.reduction].", call. = FALSE)
  }
  noise.reduction <- tolower(noise.reduction)
  allowed <- c("mean", "median", "mode", "mode_class", "off")
  if (!noise.reduction %in% allowed) {
    stop(
      "[noise.reduction] must be one of ",
      paste(shQuote(allowed), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  if (is.null(obs.req)) obs.req <- 8L
  if (!is.numeric(obs.req) || length(obs.req) != 1L ||
      !is.finite(obs.req) || obs.req < 0 || obs.req != floor(obs.req)) {
    stop("[obs.req] must be a nonnegative integer.", call. = FALSE)
  }
  obs.req <- as.integer(obs.req)
  
  order.max <- is.character(order) && length(order) == 1L &&
    !is.na(order) && tolower(order) == "max"
  
  if (!is.null(order) && !order.max) {
    if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
        order < 1 || order != floor(order)) {
      stop(
        "[order] must be NULL, 'max', or a positive integer.",
        call. = FALSE
      )
    }
    order <- as.integer(order)
  }
  
  # Explicit maximum representation. For two-axis partitioning each observation
  # is retained as its own limiting point. For XONLY, equal x values necessarily
  # share a partition and their y values are reduced coherently.
  if (order.max) {
    if (is.null(type)) {
      quadrant <- paste0("q", seq_along(x))
      prior <- rep("pq", length(x))
      PART <- data.frame(
        x = x,
        y = y,
        quadrant = quadrant,
        prior.quadrant = prior,
        stringsAsFactors = FALSE
      )
      RP <- data.frame(
        x = x,
        y = y,
        quadrant = quadrant,
        stringsAsFactors = FALSE
      )
      final.order <- length(x)
    } else {
      ux <- sort(unique(x))
      reducer <- function(z) {
        if (noise.reduction == "mean") {
          mean(z)
        } else if (noise.reduction == "median") {
          stats::median(z)
        } else if (noise.reduction %in% c("mode", "mode_class")) {
          mode_class(z)
        } else {
          gravity(z)
        }
      }
      y.by.x <- vapply(ux, function(v) reducer(y[x == v]), numeric(1L))
      match.id <- match(x, ux)
      quadrant <- paste0("q", match.id)
      PART <- data.frame(
        x = x,
        y = y,
        quadrant = quadrant,
        prior.quadrant = rep("pq", length(x)),
        stringsAsFactors = FALSE
      )
      RP <- data.frame(
        x = ux,
        y = y.by.x,
        quadrant = paste0("q", seq_along(ux)),
        stringsAsFactors = FALSE
      )
      final.order <- length(ux)
    }
    
    if (Voronoi) {
      graphics::plot(
        x, y,
        col = "steelblue",
        cex.lab = 1.5,
        xlab = x.label,
        ylab = y.label,
        mgp = c(2.5, 0.5, 0)
      )
      graphics::points(RP$x, RP$y, pch = 15, lwd = 2, col = "red")
      graphics::title(
        main = paste0("NNS Order = ", final.order),
        cex.main = 2
      )
    }
    
    return(list(
      order = as.integer(final.order),
      dt = .NNS.df(PART),
      regression.points = .NNS.df(RP)
    ))
  }
  
  n <- length(x)
  if (is.null(order)) order <- max(ceiling(log(n, 2)), 1L)
  
  out <- NNS_part_cpp(
    x = x,
    y = y,
    type = if (is.null(type)) NULL else type,
    order_in = as.integer(order),
    obs_req = obs.req,
    min_obs_stop = min.obs.stop,
    noise_reduction = noise.reduction
  )
  
  PART <- as.data.frame(out$dt, stringsAsFactors = FALSE)
  RP <- as.data.frame(out$regression.points, stringsAsFactors = FALSE)
  RP <- RP[order(RP$quadrant, method = "radix"), , drop = FALSE]
  rownames(RP) <- NULL
  
  if (is.discrete(x)) {
    finite <- is.finite(RP$x)
    RP$x[finite] <- ifelse(
      RP$x[finite] %% 1 < 0.5,
      floor(RP$x[finite]),
      ceiling(RP$x[finite])
    )
  }
  
  if (Voronoi) {
    graphics::plot(
      x, y,
      col = "steelblue",
      cex.lab = 1.5,
      xlab = x.label,
      ylab = y.label,
      mgp = c(2.5, 0.5, 0)
    )
    if (is.null(type)) {
      sh <- out$segments_h
      if (NROW(sh)) {
        graphics::segments(sh$x0, sh$y, sh$x1, sh$y, lty = 3)
      }
      sv <- out$segments_v
      if (NROW(sv)) {
        graphics::segments(sv$x, sv$y0, sv$x, sv$y1, lty = 3)
      }
    } else {
      vl <- out$vlines
      if (length(vl)) graphics::abline(v = vl, lty = 3)
    }
    graphics::points(RP$x, RP$y, pch = 15, lwd = 2, col = "red")
    graphics::title(
      main = paste0("NNS Order = ", out$order),
      cex.main = 2
    )
  }
  
  list(
    order = as.integer(out$order),
    dt = .NNS.df(PART),
    regression.points = .NNS.df(RP)
  )
}
