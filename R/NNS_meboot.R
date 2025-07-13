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
                       expand.sd = TRUE, force.clt = TRUE,
                       scl.adjustment = FALSE, sym = FALSE, elaps = FALSE,
                       digits = 6,
                       colsubj, coldata, coltimes,...)
{
  
  if(length(x)==1) return(list(x=x))
  
  type <- tolower(type)
  
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
  
  if(sum(is.na(x)) > 0) stop("You have some missing values, please address.")
  
  trim <- list(trim=trim, xmin=xmin, xmax=xmax)
  
  trimval <- if (is.null(trim$trim)) 0.1 else trim$trim
  
  n <- length(x)
  
  # Sort the original data in increasing order and
  # store the ordering index vector.
  
  xx <- sort(x)
  ordxx <- order(x)
  
  
  if(!is.null(target_drift) || !is.null(target_drift_scale)) drift <- TRUE 

  ordxx_2 <- rev(ordxx) 
  
  
  # symmetry
  
  if (sym)
  {
    xxr <- rev(xx) #reordered values
    xx.sym <- mean(xx) + 0.5*(xx - xxr) #symmetrized order stats
    xx <- xx.sym #replace order stats by symmetrized ones
  }
  
  # Compute intermediate points on the sorted series.
  
  z <- (xx[-1] + xx[-n])/2
  
  # Compute lower limit for left tail ('xmin') and
  # upper limit for right tail ('xmax').
  # This is done by computing the 'trim' (e.g. 10%) trimmed mean
  # of deviations among all consecutive observations ('dv').
  # Thus the tails are uniform distributed.
  
  dv <- abs(diff(as.numeric(x)))
  dvtrim <- mean(dv, trim=trimval)
  
  if (is.list(trim))
  {
    if (is.null(trim$xmin))
    {
      xmin <- xx[1] - dvtrim
    } else
      xmin <- trim$xmin
    
    if (is.null(trim$xmax))
    {
      xmax <- xx[n] + dvtrim
    } else
      xmax <- trim$xmax
    
    if (!is.null(trim$xmin) || !is.null(trim$xmax))
    {
      if (isTRUE(force.clt))
      {
        expand.sd <- FALSE
        force.clt <- FALSE
        # warning("expand.sd and force.clt were set to FALSE in order to ",
        #         "enforce the limits xmin/xmax.")
      }
    }
  } else {
    xmin <- xx[1] - dvtrim
    xmax <- xx[n] + dvtrim
  }
  
  
  # Compute the mean of the maximum entropy density within each
  # interval in such a way that the 'mean preserving constraint'
  # is satisfied. (Denoted as m_t in the reference paper.)
  # The first and last interval means have distinct formulas.
  # See Theil and Laitinen (1980) for details.
  
  aux <- colSums( t(cbind(xx[-c(1,2)], xx[-c(1,n)], xx[-c((n-1),n)]))*c(0.25,0.5,0.25) )
  desintxb <- c(0.75*xx[1]+0.25*xx[2], aux, 0.25*xx[n-1]+0.75*xx[n])
  
  # Generate random numbers from the [0,1] uniform interval and
  # compute sample quantiles at those points.
  
  # Generate random numbers from the [0,1] uniform interval.
  
  ensemble <- matrix(x, nrow=n, ncol=reps)
  ensemble <- apply(ensemble, 2, NNS.meboot.part,
                    n, z, xmin, xmax, desintxb, reachbnd)
  
  # So far the object 'ensemble' contains the quantiles.
  # Now give them time series dependence and heterogeneity.
  
  qseq <- apply(ensemble, 2, sort)
  
  
  # 'qseq' has monotonic series, the correct series is obtained
  # after applying the order according to 'ordxx' defined above.
  
  ensemble[ordxx,] <- qseq
  
  
  matrix2 = matrix(0, nrow=length(x), ncol = reps)
  matrix2[ordxx_2,] <- qseq
  
  # Intial search
  
  e <- c(ensemble)
  m <- c(matrix2)
  l <- length(e)
  
  func <- function(ab, d = drift, ty = type) {
    a <- ab[1]
    b <- ab[2]
    
    # Compute the adjusted ensemble
    combined <- (a * m + b * e) / (a + b)
    
    # Check correlation or dependence structure
    if (ty == "spearman" || ty == "pearson") {
      error <- abs(cor(combined, e, method = ty) - rho)
    } else if (ty == "nnsdep") {
      error <- abs(NNS.dep(combined, e)$Dependence - rho)
    } else {
      error <- abs(NNS.dep(combined, e)$Correlation - rho)
    }
    
    return(error)
  }
  
  res <- optim(c(.5, .5), func, control=list(abstol = .01))
  
  ensemble <- (res$par[1]*matrix2 + res$par[2]*ensemble) / (sum(abs(res$par)))
  
  
  # Drift
  orig_coef <- fast_lm(1:n, x)$coef
  orig_intercept <- orig_coef[1]
  orig_drift <- orig_coef[2]
  
  new_coef <- apply(ensemble, 2, function(i) fast_lm(1:n, i)$coef)
  slopes <- new_coef[2,]
  
  if(drift){
    if(!is.null(target_drift_scale)) target_drift <- orig_drift * target_drift_scale else if(is.null(target_drift)) target_drift <- orig_drift
    new_slopes <- (target_drift - slopes)
    ensemble <- ensemble + t(t(sapply(new_slopes, function(slope) cumsum(rep(slope, n)))))
    
    new_intercepts <- orig_intercept - new_coef[1,]
    ensemble <- sweep(ensemble, 2, new_intercepts, FUN = "+")
  } 
  
  
  
  if(identical(ordxx_2, ordxx)){
    if(reps>1) ensemble <- t(apply(ensemble, 1, function(x) sample(x, size = reps, replace = TRUE)))
  }
  
  
  if(expand.sd) ensemble <- NNS.meboot.expand.sd(x=x, ensemble=ensemble, ...)
  
  if(force.clt && reps > 1) ensemble <- force.clt(x=x, ensemble=ensemble)
  
  
  
  # scale adjustment
  
  if (scl.adjustment){
    zz <- c(xmin,z,xmax) #extended list of z values
    v <- diff(zz^2) / 12
    xb <- mean(x)
    s1 <- sum((desintxb - xb)^2)
    uv <- (s1 + sum(v)) / n
    desired.sd <- sd(x)
    actualME.sd <- sqrt(uv)
    if (actualME.sd <= 0)
      stop("actualME.sd<=0 Error")
    out <- desired.sd / actualME.sd
    kappa <- out - 1
    
    ensemble <- ensemble + kappa * (ensemble - xb)
  } else kappa <- NULL
  
  # Force min / max values
  if(!is.null(trim[[2]])) ensemble <- apply(ensemble, 2, function(z) pmax(trim[[2]], z))
  if(!is.null(trim[[3]])) ensemble <- apply(ensemble, 2, function(z) pmin(trim[[3]], z))
  
  
  
  
  if(is.ts(x)){
    ensemble <- ts(ensemble, frequency=frequency(x), start=start(x))
    if(reps>1) dimnames(ensemble)[[2]] <- paste("Series", 1:reps)
  } else {
    if(reps>1) dimnames(ensemble)[[2]] <- paste("Replicate", 1:reps)
  }
  
  
  
  final <- list(x=x, replicates=round(ensemble, digits = digits), ensemble=Rfast::rowmeans(ensemble), xx=xx, z=z, dv=dv, dvtrim=dvtrim, xmin=xmin,
                xmax=xmax, desintxb=desintxb, ordxx=ordxx, kappa = kappa)
  
  return(final)
}

NNS.meboot <- Vectorize(NNS.meboot, vectorize.args = c("rho", "target_drift", "target_drift_scale"))