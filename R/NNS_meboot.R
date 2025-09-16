#' NNS meboot
#'
#' Adapted maximum entropy bootstrap routine from \code{meboot} \url{https://cran.r-project.org/package=meboot}.
#'
#' @param x vector of data.
#' @param reps numeric; number of replicates to generate.
#' @param rho numeric [-1,1] (vectorized); A \code{rho} must be provided, otherwise a blank list will be returned.
#' @param type options("spearman", "pearson", "NNScor", "NNSdep"); \code{type = "spearman"}(default) dependence metric desired.
#' @param drift logical; \code{drift = TRUE} (default) preserves the drift of the original series.
#' @param target_drift numerical; \code{target_drift = NULL} (default) Specifies the desired drift when \code{drift = TRUE}, i.e. a risk-free rate of return.
#' @param target_drift_scale numerical; instead of calculating a \code{target_drift}, provide a scalar to the existing drift when \code{drift = TRUE}.
#' @param trim numeric [0,1]; The mean trimming proportion, defaults to \code{trim = 0.1}.
#' @param xmin numeric; the lower limit for the left tail.
#' @param xmax numeric; the upper limit for the right tail.
#' @param reachbnd logical; If \code{TRUE} potentially reached bounds (xmin = smallest value - trimmed mean and
#' xmax = largest value + trimmed mean) are given when the random draw happens to be equal to 0 and 1, respectively.
#' @param expand.sd logical; If \code{TRUE} the standard deviation in the ensemble is expanded. See \code{expand.sd} in \code{meboot::meboot}.
#' @param force.clt logical; If \code{TRUE} the ensemble is forced to satisfy the central limit theorem. See \code{force.clt} in \code{meboot::meboot}.
#' @param scl.adjustment logical; If \code{TRUE} scale adjustment is performed to ensure that the population variance of the transformed series equals the variance of the data.
#' @param sym logical; If \code{TRUE} an adjustment is performed to ensure that the ME density is symmetric.
#' @param elaps logical; If \code{TRUE} elapsed time during computations is displayed.
#' @param digits integer; 6 (default) number of digits to round output to.
#' @param colsubj numeric; the column in \code{x} that contains the individual index. It is ignored if the input data \code{x} is not a \code{pdata.frame} object.
#' @param coldata numeric; the column in \code{x} that contains the data of the variable to create the ensemble. It is ignored if the input data \code{x} is not a \code{pdata.frame} object.
#' @param coltimes numeric; an optional argument indicating the column that contains the times at which the observations for each individual are observed. It is ignored if the input data \code{x}
#' is not a \code{pdata.frame} object.
#' @param ... possible argument \code{fiv} to be passed to \code{expand.sd}.
#'
#' @return Returns the following row names in a matrix:
#' \itemize{
#'   \item{x} original data provided as input.
#' \item{replicates} maximum entropy bootstrap replicates.
#' \item{ensemble} average observation over all replicates.
#' \item{xx} sorted order stats (xx[1] is minimum value).
#' \item{z} class intervals limits.
#' \item{dv} deviations of consecutive data values.
#' \item{dvtrim} trimmed mean of dv.
#' \item{xmin} data minimum for ensemble=xx[1]-dvtrim.
#' \item{xmax} data x maximum for ensemble=xx[n]+dvtrim.
#' \item{desintxb} desired interval means.
#' \item{ordxx} ordered x values.
#' \item{kappa} scale adjustment to the variance of ME density.
#' \item{elaps} elapsed time.
#' }
#' 
#' @note Vectorized \code{rho} and \code{drift} parameters will not vectorize both simultaneously.  Also, do not specify \code{target_drift = NULL}.
#'
#' @references
#' \itemize{
#' \item Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank Correlations in Maximum Entropy Bootstrap and Improved Monte Carlo Simulations.  \doi{10.2139/ssrn.3621614}
#'
#' \item Vinod, H.D. (2013), Maximum Entropy Bootstrap Algorithm Enhancements.  \doi{10.2139/ssrn.2285041}
#'
#' \item Vinod, H.D. (2006), Maximum Entropy Ensembles for Time Series Inference in Economics,
#' \emph{Journal of Asian Economics}, \bold{17}(6), pp. 955-978.
#'
#' \item Vinod, H.D. (2004), Ranking mutual funds using unconventional utility theory and stochastic dominance, \emph{Journal of Empirical Finance}, \bold{11}(3), pp. 353-377.
#' }
#'
#' @examples
#' \dontrun{
#' # To generate an orthogonal rank correlated time-series to AirPassengers
#' boots <- NNS.meboot(AirPassengers, reps = 100, rho = 0, xmin = 0)
#'
#' # Verify correlation of replicates ensemble to original
#' cor(boots["ensemble",]$ensemble, AirPassengers, method = "spearman")
#'
#' # Plot all replicates
#' matplot(boots["replicates",]$replicates , type = 'l')
#'
#' # Plot ensemble
#' lines(boots["ensemble",]$ensemble, lwd = 3)
#' 
#' # Plot original
#' lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red")
#' 
#' ### Vectorized drift with a single rho
#' boots <- NNS.meboot(AirPassengers, reps = 10, rho = 0, xmin = 0, target_drift = c(1,7))
#' matplot(do.call(cbind, boots["replicates", ]), type = "l")
#' lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red")
#' 
#' ### Vectorized rho with a single target drift
#' boots <- NNS.meboot(AirPassengers, reps = 10, rho = c(0, .5, 1), xmin = 0, target_drift = 3)
#' matplot(do.call(cbind, boots["replicates", ]), type = "l")
#' lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red")
#' 
#' ### Vectorized rho with a single target drift scale
#' boots <- NNS.meboot(AirPassengers, reps = 10, rho = c(0, .5, 1), xmin = 0, target_drift_scale = 0.5)
#' matplot(do.call(cbind, boots["replicates", ]), type = "l")
#' lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red") 
#' }
#' @export

NNS.meboot <- function(x,
                       reps = 999,
                       rho = NULL,
                       type = "spearman",
                       drift = TRUE,
                       target_drift = NULL,
                       target_drift_scale = NULL,
                       trim = 0.10,
                       xmin = NULL,
                       xmax = NULL,
                       reachbnd = TRUE,
                       expand.sd = TRUE,
                       force.clt = TRUE,
                       scl.adjustment = FALSE, sym = FALSE, elaps = FALSE,
                       digits = 6,
                       colsubj, coldata, coltimes, ...){
  
  if (length(x) == 1) return(list(x = x))
  type <- tolower(type)
  if (any(class(x) %in% c("tbl","data.table"))) x <- as.vector(unlist(x))
  if (anyNA(x)) stop("You have some missing values, please address.")
  
  trim <- list(trim = trim, xmin = xmin, xmax = xmax)
  trimval <- if (is.null(trim$trim)) 0.1 else trim$trim
  n <- length(x)
  
  # --- Fit original linear trend ONCE (time order) and get residuals
  orig_lm        <- fast_lm(1:n, x)
  orig_intercept <- orig_lm$coef[1]
  orig_drift     <- orig_lm$coef[2]
  orig_res       <- orig_lm$residuals
  
  # Choose reconstruction slope (t = 1:n); for drift=FALSE baseline is flat at intercept (t = 0 fitted value)
  if (!is.null(target_drift) || !is.null(target_drift_scale)) drift <- TRUE
  if (drift) {
    if (!is.null(target_drift_scale))      target_drift <- orig_drift * target_drift_scale
    else if (is.null(target_drift))        target_drift <- orig_drift
    recon_slope <- target_drift
  } else {
    recon_slope <- 0
  }
  baseline <- orig_intercept + recon_slope * (1:n)
  
  # ===== MEBOOT CORE ON RESIDUALS ONLY =====
  # Order stats, indices, symmetry, midpoints, tails computed from residuals
  rr     <- orig_res
  xx     <- sort(rr)
  ordxx  <- order(rr)
  ordxx_2 <- rev(ordxx)
  
  if (sym) {
    xxr <- rev(xx)
    xx  <- mean(xx) + 0.5 * (xx - xxr)
  }
  
  z      <- (xx[-1] + xx[-n]) / 2
  dv     <- abs(diff(as.numeric(rr)))
  dvtrim <- mean(dv, trim = trimval)
  
  if (is.list(trim)) {
    xmin <- if (is.null(trim$xmin)) xx[1] - dvtrim else trim$xmin
    xmax <- if (is.null(trim$xmax)) xx[n] + dvtrim else trim$xmax
    if (!is.null(trim$xmin) || !is.null(trim$xmax)) {
      if (isTRUE(force.clt)) { expand.sd <- FALSE; force.clt <- FALSE }
    }
  } else { xmin <- xx[1] - dvtrim; xmax <- xx[n] + dvtrim }
  
  # Theil–Laitinen interval means on residuals
  aux <- colSums(t(cbind(xx[-c(1,2)], xx[-c(1,n)], xx[-c((n-1),n)])) * c(0.25, 0.5, 0.25))
  desintxb <- c(0.75*xx[1] + 0.25*xx[2], aux, 0.25*xx[n-1] + 0.75*xx[n])
  
  # Quantile draws from max-entropy bootstrap IN RESIDUAL SPACE
  res_mat <- matrix(rr, nrow = n, ncol = reps)
  res_mat <- apply(res_mat, 2, NNS.meboot.part, n, z, xmin, xmax, desintxb, reachbnd)  
  qseq <- apply(res_mat, 2, sort)
  res_mat[ordxx, ] <- qseq
  
  # ===== Optional dependence targeting ρ in residual space (per replicate, time-aligned) =====
  if (!is.null(rho)) {
    rho_vec <- if (length(rho) == 1L) rep(rho, reps) else rep_len(rho, reps)
    
    # Ranks of original residuals for aligned vs anti-aligned extremes
    r_o    <- rank(orig_res, ties.method = "average")
    r_anti <- max(r_o) + 1 - r_o
    
    for (i in 1:reps) {
      # start from each residual replicate column
      res_i      <- res_mat[, i]
      res_sorted <- sort(res_i)
      e <- res_sorted[r_o]     # aligned with ranks of orig_res
      m <- res_sorted[r_anti]  # anti-aligned
      
      rho_target <- rho_vec[i]
      obj <- function(ab){
        a <- ab[1]; b <- ab[2]
        comb <- (a*m + b*e) / (a + b)
        if (type %in% c("spearman","pearson")) {
          abs(cor(comb, orig_res, method = type) - rho_target)
        } else if (type == "nnsdep") {
          abs(NNS.dep(comb, orig_res)$Dependence - rho_target)   
        } else {
          abs(NNS.dep(comb, orig_res)$Correlation - rho_target)
        }
      }
      opt <- optim(c(0.5, 0.5), obj, control = list(abstol = 0.01))
      res_mat[, i] <- (opt$par[1]*m + opt$par[2]*e) / sum(abs(opt$par))
    }
  }
  
  # ===== Variance expansion ON RESIDUALS (match sd to original residuals) =====
  res_mat <- NNS.meboot.expand.sd(x = orig_res, ensemble = res_mat, ...) 
  
  # ===== Reconstruct levels: baseline + residuals =====
  ensemble <- sweep(res_mat, 1, baseline, "+")
  
  # Keep legacy “identical(ordxx_2, ordxx)” reshuffle 
  if (identical(ordxx_2, ordxx)) {
    if (reps > 1) ensemble <- t(apply(ensemble, 1, function(z) sample(z, size = reps, replace = TRUE)))
  }
  
  # Optional level scaling toward sd(x) 
  if (isTRUE(expand.sd)) {
    ensemble <- NNS.meboot.expand.sd(x = x, ensemble = ensemble, ...)   
  }
  
  # Optional CLT enforcement
  if (force.clt && reps > 1) ensemble <- force.clt(x = x, ensemble = ensemble)
  
  # Optional ME-density scale adjustment (same as before)
  if (scl.adjustment){
    zz <- c(xmin, z, xmax)
    v  <- diff(zz^2) / 12
    xb <- mean(x)
    s1 <- sum((desintxb - xb)^2)
    uv <- (s1 + sum(v)) / n
    desired.sd  <- sd(x)
    actualME.sd <- sqrt(uv)
    if (actualME.sd <= 0) stop("actualME.sd<=0 Error")
    kappa <- (desired.sd / actualME.sd) - 1
    ensemble <- ensemble + kappa * (ensemble - xb)
  } else kappa <- NULL
  
  # Enforce min / max if provided
  if (!is.null(trim[[2]])) ensemble <- apply(ensemble, 2, function(z) pmax(trim[[2]], z))
  if (!is.null(trim[[3]])) ensemble <- apply(ensemble, 2, function(z) pmin(trim[[3]], z))
  
  # ts attributes
  if (is.ts(x)) {
    ensemble <- ts(ensemble, frequency = frequency(x), start = start(x))
    if (reps > 1) dimnames(ensemble)[[2]] <- paste("Series", 1:reps)
  } else {
    if (reps > 1) dimnames(ensemble)[[2]] <- paste("Replicate", 1:reps)
  }
  
  final <- list(x = x,
                replicates = round(ensemble, digits = digits),
                ensemble = Rfast::rowmeans(ensemble),
                xx = xx, z = z, dv = dv, dvtrim = dvtrim,
                xmin = xmin, xmax = xmax, desintxb = desintxb,
                ordxx = ordxx, kappa = kappa)
  return(final)
}

NNS.meboot <- Vectorize(NNS.meboot,
                        vectorize.args = c("rho", "target_drift", "target_drift_scale"))
