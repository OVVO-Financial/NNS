#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#'
#' @param degree numeric; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Set to \code{target = mean(variable)} for classical equivalences, but does not have to be. (Vectorized)
#' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
#' @param excess_ret logical; \code{FALSE} (default)
#' @return LPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' LPM(0, mean(x), x)
#' @export

LPM <- function(degree, target, variable, excess_ret = FALSE) {
  target   <- as.numeric(target)
  variable <- as.numeric(variable)
  
  if (!excess_ret && length(target) > 1) {
    return(.Call("_NNS_LPM_CPv", degree, target, variable))
  }
  
  .Call("_NNS_LPM_RCPP", degree, target, variable, excess_ret)
  
}


#' Upper Partial Moment
#'
#' This function generates a univariate upper partial moment for any degree or target.
#'
#' @param degree numeric; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Set to \code{target = mean(variable)} for classical equivalences, but does not have to be. (Vectorized)#' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
#' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
#' @param excess_ret logical; \code{FALSE} (default)
#' @return UPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' UPM(0, mean(x), x)
#' @export

UPM <- function(degree, target, variable, excess_ret = FALSE) {
  target   <- as.numeric(target)
  variable <- as.numeric(variable)
  
  if (!excess_ret && length(target) > 1) {
    return(.Call("_NNS_UPM_CPv", degree, target, variable))
  }
  
  .Call("_NNS_UPM_RCPP", degree, target, variable, excess_ret)
  
}


#' Co‑Lower Partial Moment nD
#'
#' This function generates an n‑dimensional co‑lower partial moment (n >= 2) for any degree or target.
#'
#' @param data   A numeric matrix with observations in rows and variables in columns.
#' @param target A numeric vector, length equal to ncol(data).
#' @param degree numeric; degree for lower deviations (0 = frequency, 1 = area).
#' @param norm   logical; if \code{TRUE} (default) normalize to the maximum observed value (→ [0,1]), otherwise return the raw moment.
#' @return Numeric; the n‑dimensional co‑lower partial moment.
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(200), ncol = 4)
#' Co.LPM_nD(mat, rep(0, ncol(mat)), degree = 1, norm = FALSE)
#' }
#' @export
Co.LPM_nD <- function(data, target, degree = 0.0, norm = TRUE) {
  data   <- as.matrix(data)
  target <- as.numeric(target)
  degree <- as.numeric(degree)
  norm   <- as.logical(norm)
  
  .Call("_NNS_CoLPM_nD_RCPP", data, target, degree, norm)
}


#' Co‑Upper Partial Moment nD
#'
#' This function generates an n‑dimensional co‑upper partial moment (n >= 2) for any degree or target.
#'
#' @param data   A numeric matrix with observations in rows and variables in columns.
#' @param target A numeric vector, length equal to ncol(data).
#' @param degree numeric; degree for upper deviations (0 = frequency, 1 = area).
#' @param norm   logical; if \code{TRUE} (default) normalize to the maximum observed value (→ [0,1]), otherwise return the raw moment.
#' @return Numeric; the n‑dimensional co‑upper partial moment.
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(200), ncol = 4)
#' Co.UPM_nD(mat, rep(0, ncol(mat)), degree = 1, norm = FALSE)
#' }
#' @export
Co.UPM_nD <- function(data, target, degree = 0.0, norm = TRUE) {
  data   <- as.matrix(data)
  target <- as.numeric(target)
  degree <- as.numeric(degree)
  norm   <- as.logical(norm)
  
  .Call("_NNS_CoUPM_nD_RCPP", data, target, degree, norm)
}


#' Divergent Partial Moment nD
#'
#' This function generates the aggregate n‑dimensional divergent partial moment (n >= 2) for any degree or target.
#'
#' @param data   A numeric matrix with observations in rows and variables in columns.
#' @param target A numeric vector, length equal to ncol(data).
#' @param degree numeric; degree for upper deviations (0 = frequency, 1 = area).
#' @param norm   logical; if \code{TRUE} (default) normalize to the maximum observed value (→ [0,1]), otherwise return the raw moment.
#' @return Numeric; the n‑dimensional co‑upper partial moment.
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(200), ncol = 4)
#' DPM_nD(mat, rep(0, ncol(mat)), degree = 1, norm = FALSE)
#' }
#' @export
DPM_nD <- function(data, target, degree = 0.0, norm = TRUE) {
  data   <- as.matrix(data)
  target <- as.numeric(target)
  degree <- as.numeric(degree)
  norm   <- as.logical(norm)
  
  .Call("_NNS_DPM_nD_RCPP", data, target, degree, norm)
}



#' NNS CDF
#'
#' This function generates an empirical CDF using partial moment ratios \link{LPM.ratio}, and resulting survival, hazard and cumulative hazard functions.
#'
#' @param variable a numeric vector or data.frame of >= 2 variables for joint CDF.
#' @param degree numeric; \code{(degree = 0)} (default) is frequency, \code{(degree = 1)} is area.
#' @param target numeric; \code{NULL} (default) Must lie within support of each variable.
#' @param type options("CDF", "survival", "hazard", "cumulative hazard"); \code{"CDF"} (default) Selects type of function to return for bi-variate analysis.  Multivariate analysis is restricted to \code{"CDF"}.
#' @param plot logical; plots CDF.
#' @return Returns:
#' \itemize{
#'  \item{\code{"Function"}} a data.table containing the observations and resulting CDF of the variable.
#'  \item{\code{"target.value"}} value from the \code{target} argument.
#' }
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"  \doi{10.2139/ssrn.3007373}
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.CDF(x)
#'
#' ## Empirical CDF (degree = 0)
#' NNS.CDF(x)
#'
#' ## Continuous CDF (degree = 1)
#' NNS.CDF(x, 1)
#'
#' ## Joint CDF
#' x <- rnorm(5000) ; y <- rnorm(5000)
#' A <- cbind(x,y)
#'
#' NNS.CDF(A, 0)
#'
#' ## Joint CDF with target
#' NNS.CDF(A, 0, target = rep(0, ncol(A)))
#' }
#' @export


NNS.CDF <- function(variable,
                    degree = 0,
                    target = NULL,
                    type   = "CDF",
                    plot   = TRUE) {
  
  # — Flatten tibbles/data.tables
  if (any(class(variable) %in% c("tbl","data.table")) && ncol(variable)==1) {
    variable <- as.vector(unlist(variable))
  }
  if (any(class(variable) %in% c("tbl","data.table"))) {
    variable <- as.data.frame(variable)
  }
  
  # — Bounds check
  if (!is.null(target)) {
    if (is.null(dim(variable))||ncol(variable)==1) {
      if (target<min(variable)||target>max(variable)) stop("target out of bounds")
    } else {
      if (target[1]<min(variable[,1])||target[1]>max(variable[,1])||
          target[2]<min(variable[,2])||target[2]>max(variable[,2])) stop("target out of bounds")
    }
  }
  
  # — Validate type
  type <- tolower(type)
  if (!type%in%c("cdf","survival","hazard","cumulative hazard")) stop("invalid type")
  
  # — Axis labels
  mc <- match.call(); vc <- mc$variable
  if (is.null(dim(variable))||ncol(variable)==1) {
    vn <- deparse(vc)
  } else if (!is.null(colnames(variable))) {
    xlab <- colnames(variable)[1]; ylab <- colnames(variable)[2]
  } else {
    expr <- deparse(vc); xlab <- paste0(expr,"[,1]"); ylab <- paste0(expr,"[,2]")
  }
  
  # — Univariate branch
  if (is.null(dim(variable))||ncol(variable)==1) {
    x    <- sort(variable)
    pval <- LPM.ratio(degree,x,variable)
    DT   <- data.table::data.table(x, pval)
    colname <- switch(type,
                      cdf="CDF",
                      survival="S(x)",
                      hazard="h(x)",
                      `cumulative hazard`="H(x)")
    data.table::setnames(DT, c("x",colname))
    
    # adjust pval for survival/hazard/cumhaz
    if (type=="survival") DT[[2]] <- 1-DT[[2]]
    if (type=="hazard") {
      n <- length(x); w <- min(10,n-1)
      F <- pval
      proxy <- vapply(seq_along(x),function(i){lo<-max(1,i-w%/%2);hi<-min(n,i+w%/%2);(F[hi]-F[lo])/(x[hi]-x[lo])},numeric(1))
      fit <- NNS.reg(x, pmax(proxy,1e-10), order=NULL, n.best=1, point.est=target, plot=FALSE)
      DT[[2]] <- pmin(pmax(fit$Fitted$y.hat / pmax(1-F,1e-10),0),1e6)
    }
    if (type=="cumulative hazard") DT[[2]] <- pmax(-log(pmax(1-pval,1e-10)),0)
    
    # compute target.value
    if (is.null(target)) {
      Pv <- numeric(0)
    } else {
      Pv <- LPM.ratio(degree,target,variable)
      if (type=="survival") Pv <- 1-Pv
      if (type=="hazard") Pv <- fit$Point.est / pmax(1-pval[which.min(abs(x-target))],1e-10)
      if (type=="cumulative hazard") Pv <- NNS.reg(x,DT[[2]],order=NULL,n.best=1,point.est=target,plot=FALSE)$Point.est
    }
    
    # plotting
    if (plot) {
      plot(DT$x,DT[[2]],type="s",lwd=2,pch=19,col="steelblue",xlab=vn,ylab=colname,main=toupper(type))
      points(DT$x,DT[[2]],pch=19,col="steelblue")
      if(length(Pv)){
        segments(target,0,target,Pv,col="red",lty=2,lwd=2)
        segments(min(x),Pv,target,Pv,col="red",lty=2,lwd=2)
        points(target,Pv,pch=19,col="green")
      }
    }
    
    return(list(Function=DT,target.value=Pv))
  }
  
  # — Multivariate case (d >= 2)
  if (!is.null(dim(variable)) && ncol(variable) >= 2) {
    xlab <- colnames(variable)[1]
    ylab <- if(ncol(variable) >= 2) colnames(variable)[2] else ""
    
    # Compute joint conditional CDF using clpm_nD
    CDF <- apply(variable, 1, function(row) Co.LPM_nD(variable, row, degree = degree))
    
    # Apply transformation based on type
    if (type == "survival") {
      marginal_probs <- apply(variable, 2, function(col) LPM.ratio(degree, col, col))
      CDF <- pmax(0, pmin(1, 1 - rowSums(marginal_probs) + CDF))
    }
    
    if (type == "hazard") {
      f <- NNS.reg(variable, pmax(CDF, 1e-10), order = "max", plot = FALSE)$Fitted$y.hat
      marginals <- apply(variable, 2, function(col) LPM.ratio(degree, col, col))
      CDF <- pmax(f / pmax(1 - rowSums(marginals) + CDF, 1e-10), 0)
    }
    
    if (type == "cumulative hazard") {
      marginals <- apply(variable, 2, function(col) LPM.ratio(degree, col, col))
      CDF <- pmax(-log(pmax(1 - rowSums(marginals) + CDF, 1e-10)), 0)
    }
    
    # Target evaluation
    Pv <- numeric(0)
    if (!is.null(target)) {
      Pv <- Co.LPM_nD(variable, target, degree = degree)
      if (type == "survival") {
        marg_target <- mapply(LPM.ratio, degree, target, as.data.frame(variable))
        Pv <- max(0, min(1, 1 - sum(marg_target) + Pv))
      }
      if (type == "hazard") {
        Pv <- NNS.reg(variable, CDF, order = "max", plot = FALSE, point.est = target)$Point.est /
          pmax(1 - Pv, 1e-10)
      }
      if (type == "cumulative hazard") {
        Pv <- pmax(-log(pmax(1 - Pv, 1e-10)), 0)
      }
    }
    
    if (plot && ncol(variable) == 2) {
      x1 <- variable[, 1]; x2 <- variable[, 2]
      u1 <- LPM.ratio(degree, x1, x1)
      u2 <- LPM.ratio(degree, x2, x2)
      
      rgl::plot3d(u1, u2, CDF,
                  xlab = paste0(xlab, " uniform"), ylab = paste0(ylab, " uniform"), zlab = toupper(type),
                  col  = "steelblue", pch = 19, box = FALSE)
      
      if (length(Pv)) {
        ut1 <- LPM.ratio(degree, target[1], x1)
        ut2 <- LPM.ratio(degree, target[2], x2)
        
        # Target point (green)
        rgl::points3d(ut1, ut2, Pv, col = "green", pch = 19)
        
        # Horizontal segment along x at level Pv
        rgl::segments3d(
          x = c(min(u1), ut1),
          y = c(ut2,     ut2),
          z = c(Pv,      Pv),
          col = "red", lwd = 2, lty = "dashed"
        )
        rgl::text3d(ut1,min(u2), Pv, 
                    text = paste0("x = ", round(target[1], 3)), 
                    col = "red", pos = 2, cex = 0.9)
        
        # Horizontal segment along y at level Pv
        rgl::segments3d(
          x = c(ut1,     ut1),
          y = c(min(u2), ut2),
          z = c(Pv,      Pv),
          col = "red", lwd = 2, lty = "dashed"
        )
        rgl::text3d(min(u1), ut2, Pv, 
                    text = paste0("y = ", round(target[2], 3)), 
                    col = "red", pos = 2, cex = 0.9)
        
        # Final segment to CDF axis (min u1, min u2, Pv)
        rgl::segments3d(
          x = c(ut1,      max(u1)),
          y = c(ut2,      max(u2)),
          z = c(Pv,       Pv),
          col = "red", lwd = 2, lty = "dashed"
        )
        rgl::text3d(max(u1), max(u2), Pv,
                    text = paste0("CDF = ", round(Pv, 4)),
                    col = "red", pos = 2, cex = 0.9)
      }
    }
    
    outDT <- data.table::data.table(variable, CDF = CDF)
    return(list(Function = outDT, target.value = Pv))(list(Function = outDT, target.value = Pv))(list(Function=outDT,target.value=Pv))
  }
}




#' NNS moments
#'
#' This function returns the first 4 moments of the distribution.
#'
#' @param x a numeric vector.
#' @param population logical; \code{TRUE} (default) Performs the population adjustment.  Otherwise returns the sample statistic.
#' @return Returns:
#' \itemize{
#'  \item{\code{"$mean"}} mean of the distribution.
#'  \item{\code{"$variance"}} variance of the distribution.
#'  \item{\code{"$skewness"}} skewness of the distribution.
#'  \item{\code{"$kurtosis"}} excess kurtosis of the distribution.
#' }
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.moments(x)
#' }
#' @export

NNS.moments <- function(x, population = TRUE){
  n <- length(x)
  mean <- UPM(1, 0, x) - LPM(1, 0, x)
  variance <- (UPM(2, mean(x), x) + LPM(2, mean(x), x))
  skew_base <- (UPM(3,mean(x),x) - LPM(3,mean(x),x))
  kurt_base <- (UPM(4,mean(x),x) + LPM(4,mean(x),x))
  
  if(population){
    skewness <- skew_base / variance^(3/2)
    kurtosis <- (kurt_base / variance^2) - 3
  } else {
    skewness <- (n / ((n-1)*(n-2))) * ((n*skew_base) / variance^(3/2))
    kurtosis <- ((n * (n+1)) / ((n-1)*(n-2)*(n-3))) * ((n*kurt_base) / (variance * (n / (n - 1)))^2) - ( (3 * ((n-1)^2)) / ((n-2)*(n-3)))
    variance <- variance * (n / (n - 1))
  }
  
  return(list("mean" = mean,
              "variance" = variance,
              "skewness" = skewness,
              "kurtosis" = kurtosis))
}




#' Partial Moment Matrix
#' @name PM.matrix
#' @title Partial Moment Matrix
#' @description Builds a list containing CUPM, DUPM, DLPM, CLPM and the overall covariance matrix.
#' @param LPM_degree numeric; lower partial moment degree (0 = freq, 1 = area).
#' @param UPM_degree numeric; upper partial moment degree (0 = freq, 1 = area).
#' @param target numeric vector; thresholds for each column (defaults to colMeans).
#' @param variable numeric matrix or data.frame.
#' @param pop_adj logical; TRUE adjusts population vs. sample moments.
#' @param norm logical; default FALSE. If TRUE, each quadrant matrix is cell-wise normalized so their sum is 1 at each (i,j).
#' @return A list: $cupm, $dupm, $dlpm, $clpm, $cov.matrix.
#' @examples
#' set.seed(123)
#' A <- cbind(rnorm(100), rnorm(100), rnorm(100))
#' PM.matrix(1, 1, NULL, A, TRUE)          # uses norm = FALSE by default
#' PM.matrix(1, 1, NULL, A, TRUE, TRUE)    # enable normalization
#' @export
PM.matrix <- function(LPM_degree, UPM_degree, target, variable, pop_adj, norm = FALSE) {
  .Call(`_NNS_PMMatrix_RCPP`, LPM_degree, UPM_degree, target, variable, pop_adj, norm)
}

