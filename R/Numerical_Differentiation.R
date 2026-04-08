#' NNS Numerical Differentiation
#'
#' Determines numerical derivative of a given univariate function using projected secant lines on the y-axis. These projected points infer finite steps \code{h}, in the finite step method.
#'
#' @param f an expression or call or a formula with no lhs.
#' @param point numeric; Point to be evaluated for derivative of a given function \code{f}.
#' @param h numeric [0, ...]; Initial step for secant projection. Defaults to \code{(h = abs(point) * 0.1 + 0.01)}.
#' @param tol numeric; Sets the tolerance for the stopping condition of the inferred \code{h}. Defaults to \code{(tol = 1e-10)}.
#' @param max.iter integer; \code{NULL} (default) Maximum number of bisection iterations. \code{NULL} sets the limit to \code{100L}. For noisy functions the bisection may stall before \code{tol} is reached; \code{max.iter} provides a hard upper bound.
#' @param digits numeric; Sets the number of digits specification of the output. Defaults to \code{(digits = 12)}.
#' @param print.trace logical; \code{FALSE} (default) Displays each iteration, lower y-intercept, upper y-intercept and inferred \code{h}.
#' @param plot logical; plots range, secant lines and y-intercept convergence.
#' @return Returns a matrix of values, intercepts, derivatives, inferred step sizes for multiple methods of estimation.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' \dontrun{
#' f <- function(x) sin(x) / x
#' NNS.diff(f, 4.1)
#'
#' ## Noisy function with explicit iteration cap
#' f_noisy <- function(x) sin(x) + rnorm(1, 0, 0.001)
#' NNS.diff(f_noisy, 1.0, max.iter = 100)
#' }
#' @export

NNS.diff <- function(f, point, h = abs(point) * 0.1 + 0.01, tol = 1e-10, max.iter = NULL,
                     digits = 12, print.trace = FALSE, plot = FALSE){
  
  if(!is.function(f)) stop("'f' must be a function.")
  if(!is.numeric(point) || length(point) != 1L || is.na(point) || !is.finite(point)) {
    stop("'point' must be a single finite numeric value.")
  }
  if(!is.numeric(h) || length(h) != 1L || is.na(h) || !is.finite(h) || h <= 0) {
    stop("'h' must be a single finite numeric value > 0.")
  }
  if(!is.numeric(tol) || length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single finite numeric value > 0.")
  }
  if(is.null(max.iter)) {
    max.iter <- 100L
  } else {
    if(!is.numeric(max.iter) || length(max.iter) != 1L || is.na(max.iter) || !is.finite(max.iter) || max.iter < 1) {
      stop("'max.iter' must be a single finite integer >= 1.")
    }
    max.iter <- as.integer(max.iter)
  }
  if(!is.numeric(digits) || length(digits) != 1L || is.na(digits) || !is.finite(digits) || digits < 0) {
    stop("'digits' must be a single finite numeric value >= 0.")
  }
  if(!is.logical(print.trace) || length(print.trace) != 1L || is.na(print.trace)) {
    stop("'print.trace' must be a single TRUE or FALSE.")
  }
  if(!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be a single TRUE or FALSE.")
  }
  
  
  
  Finite.step <- function(f, point, h){
    f.x <- f(point)
    f.x.h.min <- f(point - h)
    f.x.h.pos <- f(point + h)
    
    neg.step <- (f.x - f.x.h.min) / h
    pos.step <- (f.x.h.pos - f.x) / h
    
    c("f(x-h)" = neg.step,
      "f(x+h)" = pos.step,
      "Averaged Finite Step" = mean(c(neg.step, pos.step)))
  }
  
  safe.range <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if(length(x) == 0L) return(c(-1, 1))
    r <- range(x)
    if(r[1] == r[2]) r <- r + c(-1, 1)
    r
  }
  
  Bs <- numeric()
  Bl <- numeric()
  Bu <- numeric()
  
  f.x <- f(point)
  if(!is.numeric(f.x) || length(f.x) != 1L || is.na(f.x) || !is.finite(f.x)) {
    stop("'f(point)' must return a single finite numeric value.")
  }
  
  f.x.h.lower <- f(point - h)
  f.x.h.upper <- f(point + h)
  
  if(any(!is.finite(c(f.x.h.lower, f.x.h.upper)))) {
    stop("'f(point +/- h)' must return finite numeric values.")
  }
  
  left.slope  <- (f.x - f.x.h.lower) / h
  right.slope <- (f.x.h.upper - f.x) / h
  
  B1 <- f.x - left.slope * point
  B2 <- f.x - right.slope * point
  
  low.B <- min(c(B1, B2))
  high.B <- max(c(B1, B2))
  
  lower.B <- low.B
  upper.B <- high.B
  
  # ---------------------------------------------------------------------------
  # FIX 1:
  # If both projected secants share the same intercept, that usually means
  # the local slope is already identified, not that the derivative fails.
  # Return the common secant slope and associated diagnostics.
  # ---------------------------------------------------------------------------
  if(isTRUE(all.equal(lower.B, upper.B, tolerance = .Machine$double.eps^0.5))) {
    
    initial.fs <- c(
      "f(x-h)" = left.slope,
      "f(x+h)" = right.slope,
      "Averaged Finite Step" = mean(c(left.slope, right.slope))
    )
    
    slope <- mean(c(left.slope, right.slope))
    inferred.h <- 0
    i <- 0L
    converged <- TRUE
    termination.code <- 0L
    final.B <- B1
    
    return(round(
      as.matrix(
        c("Value of f(x) at point" = unname(f.x),
          "Final y-intercept (B)" = unname(final.B),
          "DERIVATIVE" = unname(slope),
          "Inferred h" = unname(inferred.h),
          "iterations" = unname(i),
          "converged" = unname(as.integer(converged)),
          "termination.code" = unname(termination.code),
          "Initial h finite step: f(x-h)" = unname(initial.fs["f(x-h)"]),
          "Initial h finite step: f(x+h)" = unname(initial.fs["f(x+h)"]),
          "Initial h averaged finite step" = unname(initial.fs["Averaged Finite Step"]),
          "Inferred h finite step: f(x-h)" = NA_real_,
          "Inferred h finite step: f(x+h)" = NA_real_,
          "Inferred h averaged finite step" = NA_real_,
          "Complex Step Derivative (Inferred h)" = NA_real_)
      ),
      digits
    ))
  }
  
  new.B <- mean(c(lower.B, upper.B))
  i <- 1L
  converged <- FALSE
  termination.code <- 2L
  inferred.h <- NA_real_
  
  while(i >= 1L){
    
    Bl[i] <- lower.B
    Bu[i] <- upper.B
    Bs[i] <- new.B
    
    new.f <- function(x) -f.x + ((f.x - f(point - x)) / x) * point + new.B
    
    inferred.h <- tryCatch(
      uniroot(new.f, c(-2 * h, 2 * h), extendInt = "yes")$root,
      error = function(e) NA_real_
    )
    
    if(print.trace) {
      print(c("Iteration" = as.integer(i),
              "h" = inferred.h,
              "Lower B" = lower.B,
              "Upper B" = upper.B))
    }
    
    if(!is.finite(inferred.h)) {
      termination.code <- 2L
      break
    }
    
    if(abs(inferred.h) < tol) {
      converged <- TRUE
      termination.code <- 0L
      break
    }
    
    if(i >= max.iter) {
      termination.code <- 1L
      break
    }
    
    if(B1 == high.B){
      if(sign(inferred.h) < 0) {
        lower.B <- new.B
      } else {
        upper.B <- new.B
      }
    } else {
      if(sign(inferred.h) < 0) {
        upper.B <- new.B
      } else {
        lower.B <- new.B
      }
    }
    
    new.B <- mean(c(lower.B, upper.B))
    i <- i + 1L
  }
  
  final.B <- mean(c(upper.B, lower.B))
  
  if(is.finite(inferred.h)) inferred.h <- abs(inferred.h)
  
  if(abs(point) < .Machine$double.eps^0.5) {
    slope <- mean(Finite.step(f, point, h)[c("f(x-h)", "f(x+h)")])
  } else {
    slope <- (f.x - final.B) / point
  }
  
  complex.step <- NA_real_
  if(is.finite(inferred.h) && inferred.h != 0) {
    z <- complex(real = point, imaginary = inferred.h)
    f.z <- tryCatch(f(z), error = function(e) NA_complex_)
    if(length(f.z) == 1L && !is.na(f.z)) {
      complex.step <- Im(f.z) / Im(z)
    }
  }
  
  initial.fs <- Finite.step(f, point, h)
  
  inferred.fs <- if(is.finite(inferred.h) && inferred.h != 0) {
    Finite.step(f, point, inferred.h)
  } else {
    c("f(x-h)" = NA_real_,
      "f(x+h)" = NA_real_,
      "Averaged Finite Step" = NA_real_)
  }
  
  if(plot) {
    original.par <- par(no.readonly = TRUE)
    on.exit(par(original.par), add = TRUE)
    
    par(mfrow = c(1, 3))
    
    x.seq.wide <- seq(point - (100 * h), point + (100 * h), length.out = 1000L)
    y.seq.wide <- suppressWarnings(tryCatch(f(x.seq.wide), error = function(e) rep(NA_real_, length(x.seq.wide))))
    ylim1 <- safe.range(c(B1, B2, y.seq.wide))
    
    plot(f,
         xlim = c(min(c(point - (100 * h), point + (100 * h), 0)),
                  max(c(point - (100 * h), point + (100 * h), 0))),
         col = "azure4",
         ylab = "f(x)",
         lwd = 2,
         ylim = ylim1,
         main = "f(x) and initial y-intercept range")
    abline(h = 0, v = 0, col = "grey")
    points(point, f.x, pch = 19, col = "green")
    points(point - h, f.x.h.lower, col = ifelse(B1 == high.B, "steelblue", "red"), pch = 19)
    points(point + h, f.x.h.upper, col = ifelse(B1 == high.B, "red", "steelblue"), pch = 19)
    points(x = rep(0, 2), y = c(B1, B2),
           col = c(ifelse(B1 == high.B, "steelblue", "red"),
                   ifelse(B1 == high.B, "red", "steelblue")),
           pch = 1)
    segments(0, B1, point - h, f.x.h.lower, col = ifelse(B1 == high.B, "steelblue", "red"), lty = 2)
    segments(0, B2, point + h, f.x.h.upper, col = ifelse(B1 == high.B, "red", "steelblue"), lty = 2)
    
    plot(f,
         col = "azure4",
         ylab = "f(x)",
         lwd = 3,
         main = "f(x) narrowed range and secant lines",
         xlim = c(min(c(point - h, point + h, 0)),
                  max(c(point + h, point - h, 0))),
         ylim = safe.range(c(B1, B2, f.x.h.lower, f.x.h.upper)))
    abline(h = 0, v = 0, col = "grey")
    points(point, f.x, pch = 19, col = "red")
    points(point - h, f.x.h.lower, col = ifelse(B1 == high.B, "steelblue", "red"), pch = 19)
    points(point + h, f.x.h.upper, col = ifelse(B1 == high.B, "red", "steelblue"), pch = 19)
    points(point, f.x, pch = 19, col = "green")
    segments(0, B1, point - h, f.x.h.lower, col = ifelse(B1 == high.B, "steelblue", "red"), lty = 2)
    segments(0, B2, point + h, f.x.h.upper, col = ifelse(B1 == high.B, "red", "steelblue"), lty = 2)
    points(x = rep(0, 2), y = c(B1, B2),
           col = c(ifelse(B1 == high.B, "steelblue", "red"),
                   ifelse(B1 == high.B, "red", "steelblue")),
           pch = 1)
    
    plot(Bs,
         ylim = safe.range(c(Bl, Bu)),
         xlab = "Iterations",
         ylab = "y-intercept",
         col = "green",
         pch = 19,
         main = "Iterated range of y-intercept")
    points(Bl, col = "red")
    points(Bu, col = "steelblue")
    legend("topright",
           c("Upper y-intercept", "Lower y-intercept", "Mean y-intercept"),
           col = c("steelblue", "red", "green"),
           pch = c(1, 1, 19),
           bty = "n")
  }
  
  round(
    as.matrix(
      c("Value of f(x) at point" = unname(f.x),
        "Final y-intercept (B)" = unname(final.B),
        "DERIVATIVE" = unname(slope),
        "Inferred h" = unname(inferred.h),
        "iterations" = unname(i),
        "converged" = unname(as.integer(converged)),
        "termination.code" = unname(termination.code),
        "Initial h finite step: f(x-h)" = unname(initial.fs["f(x-h)"]),
        "Initial h finite step: f(x+h)" = unname(initial.fs["f(x+h)"]),
        "Initial h averaged finite step" = unname(initial.fs["Averaged Finite Step"]),
        "Inferred h finite step: f(x-h)" = unname(inferred.fs["f(x-h)"]),
        "Inferred h finite step: f(x+h)" = unname(inferred.fs["f(x+h)"]),
        "Inferred h averaged finite step" = unname(inferred.fs["Averaged Finite Step"]),
        "Complex Step Derivative (Inferred h)" = unname(complex.step))
    ),
    digits
  )
}
