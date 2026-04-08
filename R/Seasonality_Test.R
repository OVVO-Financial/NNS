#' NNS Seasonality Test
#'
#' Seasonality test based on the coefficient of variation for the variable and lagged component series.  A result of 1 signifies no seasonality present.
#'
#' @param variable a numeric vector.
#' @param modulo integer(s); NULL (default) Used to find the nearest multiple(s) in the reported seasonal period.
#' @param mod.only logical; \code{TRUE} (default) Limits the number of seasonal periods returned to the specified \code{modulo}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference.
#' @return Returns a matrix of all periods exhibiting less coefficient of variation than the variable with \code{"all.periods"}; and the single period exhibiting the least coefficient of variation versus the variable with \code{"best.period"}; as well as a vector of \code{"periods"} for easy call into \link{NNS.ARMA.optim}.  If no seasonality is detected, \code{NNS.seas} will return ("No Seasonality Detected").
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## To call strongest period based on coefficient of variation:
#' NNS.seas(x, plot = FALSE)$best.period
#'
#' ## Using modulos for logical seasonal inference:
#' NNS.seas(x, modulo = c(2,3,5,7), plot = FALSE)
#' }
#' @export



NNS.seas <- function(variable,
                     modulo = NULL,
                     mod.only = TRUE,
                     plot = TRUE) {
  # API per NNS manual / Rd (arguments & defaults)  :contentReference[oaicite:3]{index=3}
  # Coerce tbl/data.table to numeric vector (repo reference)  :contentReference[oaicite:4]{index=4}
  if (any(class(variable) %in% c("tbl", "data.table"))) {
    variable <- as.vector(unlist(variable, use.names = FALSE))
  }
  if (!is.numeric(variable)) stop("Variable must be numeric")
  if (anyNA(variable)) stop("You have some missing values, please address.")
  if (any(is.infinite(variable))) stop("Infinite values not allowed")
  
  ans <- NNS_seas_cpp(
    variable           = variable,
    modulo             = if (is.null(modulo)) NULL else as.integer(modulo),
    mod_only           = isTRUE(mod.only)
  )
  
  # Plot (diagnostic)
  if (isTRUE(plot)) {
    M <- ans$all.periods
    
    # nothing to plot
    if (is.null(M) || nrow(M) == 0L) return(ans)
    
    mean_var <- mean(variable)
    if (mean_var != 0) {
      overall_cv <- abs(stats::sd(variable) / mean_var)
    } else {
      # fallback carried back from C++
      overall_cv <- M$`Variable.Coefficient.of.Variation`[1L]
    }
    overall_cv <- as.numeric(overall_cv)
    
    # Predictive strength in [0,1] with guards
    n <- nrow(M)
    if (is.finite(overall_cv) && overall_cv > 0) {
      strength <- 1 - (M[["Coefficient.of.Variation"]] / overall_cv)
      strength <- pmin(pmax(strength, 0), 1)
      steel_pal <- grDevices::colorRampPalette(
        c("steelblue1", "steelblue2", "steelblue3", "steelblue4"))
      palette_ <- steel_pal(100L)
      idx <- pmax(1L, as.integer(round(strength * 99L + 1L)))
      point_colors <- palette_[idx]
    } else {
      point_colors <- rep("steelblue3", n)
    }
    
    # y-limits: nonnegative and wide enough to show the reference line if finite
    ymax <- if (is.finite(overall_cv) && overall_cv > 0) 2 * overall_cv else max(M[["Coefficient.of.Variation"]], 1)
    ylim <- c(0, ymax)
    
    plot(M[["Period"]], M[["Coefficient.of.Variation"]],
         xlab = "Period",
         ylab = "Component Series CV",
         main = "Seasonality Detection via Predictive Power\n(Lower CV = Tighter Distribution = More Predictable)",
         ylim = ylim,
         col = point_colors, pch = 19)
    
    # highlight best period (table is keyed ascending by CV)
    points(M[["Period"]][1L], M[["Coefficient.of.Variation"]][1L],
           pch = 19, col = "red", cex = 1.5)
    
    # Reference CV line and centered label (only when finite)
    if (is.finite(overall_cv)) {
      abline(h = overall_cv, col = "red", lty = 2)
      usr <- graphics::par("usr")
      xmid <- mean(usr[1:2])
      graphics::text(xmid, overall_cv,
                     labels = "Overall Series CV\n(Predictive Power Threshold)",
                     adj = c(0.5, 0.5), col = "red", xpd = NA)
    }
  }
  
  
  # Return results
  ans
}



