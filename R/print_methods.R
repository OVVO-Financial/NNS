# Compact, data.table-style printing for NNS table returns.
#
# Many NNS routines return potentially large tables, or lists containing them.
# Previously these were data.tables, whose print method shows only the head/tail
# and dimensions, avoiding flooding the console. To restore that behaviour
# uniformly without re-adding data.table, a thin S3 class is attached to every
# returned data.frame / matrix:
#   * "NNS.data.frame" inherits from data.frame
#   * "NNS.matrix"     inherits from matrix/array
# so every operation (indexing, $, is.data.frame(), is.matrix(), as.data.frame(),
# arithmetic, matrix algebra) is unchanged -- only the auto-print is compacted.
#
# .NNS.out() is the single entry point used at every user-facing return: it
# recurses through lists and tags each data.frame / matrix leaf, so no individual
# element can be missed.


# Tag one tabular object (data.frame or matrix). Idempotent; no-op otherwise.
.NNS.df <- function(x) {
  if (is.data.frame(x)) {
    if (!inherits(x, "NNS.data.frame")) {
      class(x) <- c("NNS.data.frame", class(x))
    }
  } else if (is.matrix(x)) {
    if (!inherits(x, "NNS.matrix")) {
      class(x) <- c("NNS.matrix", "matrix", "array")
    }
  }
  
  x
}


# Recursively tag every data.frame / matrix within a returned object. Only plain
# unclassed lists are descended -- S3 objects (hclust, lm, ts, ...) are treated
# as leaves and left untouched, so their internals are never altered.
#
# Important: use x[i] <- list(value), not x[[i]] <- value, because assigning NULL
# with [[<- deletes the list element and can cause "subscript out of bounds".
.NNS.out <- function(x) {
  if (is.data.frame(x) || is.matrix(x)) return(.NNS.df(x))
  
  if (is.list(x) && is.null(attr(x, "class"))) {
    for (i in seq_along(x)) {
      x[i] <- list(.NNS.out(x[[i]]))
    }
    return(x)
  }
  
  x
}


# Normalize print-row arguments.
.NNS.print_n <- function(n, default = 6L) {
  if (length(n) == 0L) return(default)
  
  n <- suppressWarnings(as.numeric(n[1L]))
  
  if (is.na(n)) return(default)
  if (is.infinite(n)) return(Inf)
  
  max(0L, as.integer(n))
}


# Compact print helper:
#   * if small: print whole object normally
#   * if large: print head, spacer, tail, footer
# The object itself is never truncated or coerced; this is print-layer only.
.NNS.print_compact <- function(y, nr, nc, n = 6L, tailn = n,
                               more, spacer = "---", ...) {
  
  n <- .NNS.print_n(n, default = 6L)
  tailn <- .NNS.print_n(tailn, default = n)
  
  if (is.null(nr) || nr <= 0L ||
      is.infinite(n) || is.infinite(tailn) ||
      nr <= n + tailn) {
    print(y, ...)
    return(invisible(y))
  }
  
  head_lines <- if (n > 0L) {
    capture.output(print(utils::head(y, n), ...))
  } else {
    character()
  }
  
  tail_lines <- if (tailn > 0L) {
    capture.output(print(utils::tail(y, tailn), ...))
  } else {
    character()
  }
  
  # Avoid repeating the column header before the tail block.
  if (length(head_lines) > 0L &&
      length(tail_lines) > 1L &&
      trimws(head_lines[1L]) == trimws(tail_lines[1L])) {
    tail_lines <- tail_lines[-1L]
  }
  
  footer <- if (identical(n, tailn)) {
    sprintf(
      "--- [ %d rows x %d cols ]; showing first and last %d. Use %s for all. ---",
      nr, nc, n, more
    )
  } else {
    sprintf(
      "--- [ %d rows x %d cols ]; showing first %d and last %d. Use %s for all. ---",
      nr, nc, n, tailn, more
    )
  }
  
  cat(paste(c(head_lines, spacer, tail_lines, footer), collapse = "\n"),
      "\n", sep = "")
  
  invisible(y)
}


#' @export
print.NNS.data.frame <- function(x, n = 6L, tailn = n, ...) {
  y <- x
  class(y) <- setdiff(class(y), "NNS.data.frame")   # plain data.frame print
  
  .NNS.print_compact(
    y,
    nr = nrow(x),
    nc = ncol(x),
    n = n,
    tailn = tailn,
    more = "as.data.frame(x)",
    ...
  )
  
  invisible(x)
}


#' @export
print.NNS.matrix <- function(x, n = 6L, tailn = n, ...) {
  y <- x
  class(y) <- NULL                                  # plain matrix print
  
  .NNS.print_compact(
    y,
    nr = nrow(x),
    nc = ncol(x),
    n = n,
    tailn = tailn,
    more = "as.matrix(x)",
    ...
  )
  
  invisible(x)
}


# Grouped gravity consolidation helper: returns gravity(y) within each x group,
# broadcast back in original order (drop-in for ave(y, x, FUN = gravity)).
# Fast path: when every x is unique, gravity() of each singleton group is the
# value itself, so the per-group C-call gravity work is skipped entirely --
# the common case for regression points, making consolidation faster than the
# previous data.table path which still invoked gravity per group.
.grouped_gravity <- function(x, y) {
  if (!anyDuplicated(x)) return(y)
  ave(y, x, FUN = gravity)
}