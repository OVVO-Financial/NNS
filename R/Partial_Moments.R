#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#'
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Set to \code{target = mean(variable)} for classical equivalences, but does not have to be. (Vectorized)
#' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
#' @param excess_ret Logical; \code{FALSE} (default)
#' @return LPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' LPM(0, mean(x), x)
#' @export
LPM <- function(degree, target, variable, excess_ret = FALSE) {
  .Call(`_NNS_LPM_RCPP`, degree, target, variable, excess_ret)
}

LPM <- Vectorize(LPM, vectorize.args="target", USE.NAMES = FALSE)


#' Upper Partial Moment
#'
#' This function generates a univariate upper partial moment for any degree or target.
#'
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Set to \code{target = mean(variable)} for classical equivalences, but does not have to be. (Vectorized)#' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
#' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
#' @param excess_ret Logical; \code{FALSE} (default)
#' @return UPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' UPM(0, mean(x), x)
#' @export
UPM <- function(degree, target, variable, excess_ret = FALSE) {
  .Call(`_NNS_UPM_RCPP`, degree, target, variable, excess_ret)
}

UPM <- Vectorize(UPM, vectorize.args="target", USE.NAMES = FALSE)


#' NNS CDF
#'
#' This function generates an empirical CDF using partial moment ratios \link{LPM.ratio}, and resulting survival, hazard and cumulative hazard functions.
#'
#' @param variable a numeric vector or data.frame of 2 variables for joint CDF.
#' @param degree integer; \code{(degree = 0)} (default) is frequency, \code{(degree = 1)} is area.
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
#' NNS.CDF(A, 0, target = c(0,0))
#' }
#' @export


NNS.CDF <- function(variable, degree = 0, target = NULL, type = "CDF", plot = TRUE){
  if(any(class(variable) %in% c("tbl", "data.table")) && dim(variable)[2] == 1){ 
    variable <- as.vector(unlist(variable))
  }
  if(any(class(variable) %in% c("tbl", "data.table"))){
    variable <- as.data.frame(variable)
  }
  
  if(!is.null(target)){
    if(is.null(dim(variable)) || dim(variable)[2] == 1){
      if(target < min(variable) || target > max(variable)){
        stop("Please make sure target is within the observed values of variable.")
      }
    } else {
      if(target[1] < min(variable[,1]) || target[1] > max(variable[,1])){
        stop("Please make sure target 1 is within the observed values of variable 1.")
      }
      if(target[2] < min(variable[,2]) || target[2] > max(variable[,2])){
        stop("Please make sure target 2 is within the observed values of variable 2.")
      }
    }
  }
  
  type <- tolower(type)
  if(!(type %in% c("cdf", "survival", "hazard", "cumulative hazard"))){
    stop(paste("Please select a type from: ", "`CDF`, ", "`survival`, ",  "`hazard`, ", "`cumulative hazard`"))
  }
  
  # Univariate Case
  if(is.null(dim(variable)) || dim(variable)[2] == 1){
    overall_target <- sort(variable)
    x <- overall_target
    CDF <- LPM.ratio(degree, overall_target, variable)
    values <- cbind.data.frame(sort(variable), CDF)
    colnames(values) <- c(deparse(substitute(variable)), "CDF")
    if(!is.null(target)){
      P <- LPM.ratio(degree, target, variable)
    } else {
      P <- NULL
    }
    ylabel <- "Probability"
    if(type == "survival"){
      CDF <- 1 - CDF
      P <- 1 - P
      ylabel <- "S(x)"
    } else if(type == "hazard"){
      n <- length(x)
      window <- min(10, n-1)
      f_proxy <- numeric(n)
      for(i in 1:n){
        start <- max(1, i - window %/% 2)
        end <- min(n, i + window %/% 2)
        dx <- x[end] - x[start]
        f_proxy[i] <- (CDF[end] - CDF[start]) / dx
      }
      f_proxy <- pmax(f_proxy, 1e-10)
      reg_fit <- NNS.reg(x, f_proxy, order = NULL, n.best = 1, point.est = if(!is.null(target)) target else NULL, plot = FALSE)
      dens <- pmax(reg_fit$Fitted$y.hat, 1e-10)
      S <- pmax(1e-10, 1 - CDF)
      CDF <- dens / S
      CDF <- pmin(CDF, 1e6)
      CDF <- pmax(0, CDF)
      ylabel <- "h(x)"
      if(!is.null(target)){
        P <- reg_fit$Point.est / S[which.min(abs(x - target))]
        P <- min(P, 1e6)
        P <- max(0, P)
        CDF[is.infinite(CDF)] <- P
      }
    } else if(type == "cumulative hazard"){
      S <- pmax(1e-10, 1 - CDF)
      CDF <- -log(S)
      CDF <- pmax(0, CDF)
      if(!is.null(target)){
        reg_fit <- NNS.reg(x, CDF, order = NULL, n.best = 1, point.est = target, plot = FALSE)
        P <- reg_fit$Point.est
      }
      ylabel <- "H(x)"
    }
    if(plot){
      plot(x, CDF, pch = 19, col = 'steelblue', xlab = deparse(substitute(variable)), ylab = ylabel, main = toupper(type), type = "s", lwd = 2)
      points(x, CDF, pch = 19, col = 'steelblue')
      lines(x, CDF, lty = 2, col = 'steelblue')
      if(!is.null(target)){
        segments(target, 0, target, P, col = "red", lwd = 2, lty = 2)
        segments(min(variable), P, target, P, col = "red", lwd = 2, lty = 2)
        points(target, P, col = "green", pch = 19)
        mtext(text = round(P, 4), col = "red", side = 2, at = P, las = 2)
        mtext(text = round(target, 4), col = "red", side = 1, at = target, las = 1)
      }
    }
    values <- data.table::data.table(cbind.data.frame(x, CDF))
    colnames(values) <- c(deparse(substitute(variable)), ylabel)
    return(
      list(
        "Function" = values,
        "target.value" = P
      )
    )
  } 
  # Bivariate Case
  else {
    overall_target_1 <- variable[,1]
    overall_target_2 <- variable[,2]
    
    sorted_indices <- order(variable[,1], variable[,2])
    sorted_x <- variable[sorted_indices, 1]
    sorted_y <- variable[sorted_indices, 2]
    
    joint_cdf <- (
      Co.LPM(degree, overall_target_1, overall_target_2, overall_target_1, overall_target_2) /
        (
          Co.LPM(degree, overall_target_1, overall_target_2, overall_target_1, overall_target_2) +
            Co.UPM(degree, overall_target_1, overall_target_2, overall_target_1, overall_target_2) +
            D.UPM(degree, degree, overall_target_1, overall_target_2, overall_target_1, overall_target_2) +
            D.LPM(degree, degree, overall_target_1, overall_target_2, overall_target_1, overall_target_2)
        )
    )
    CDF <- joint_cdf[sorted_indices]
    
    marginal_X <- LPM.ratio(degree, sorted_x, overall_target_1)
    marginal_Y <- LPM.ratio(degree, sorted_y, overall_target_2)
    
    ylabel <- "Probability"
    if(type == "survival"){
      CDF <- 1 - marginal_X - marginal_Y + CDF
      CDF <- pmax(0, pmin(1, CDF))
      ylabel <- "S(x, y)"
    } else if(type == "hazard"){
      data_points <- data.frame(x = sorted_x, y = sorted_y)
      dens_proxy <- pmax(joint_cdf[sorted_indices], 1e-10)
      reg_fit <- NNS.reg(data_points, dens_proxy, order = "max", 
                         point.est = if(!is.null(target)) data.frame(x = target[1], y = target[2]) else NULL, 
                         plot = FALSE)
      dens <- pmax(reg_fit$Fitted$y.hat, 1e-10)
      S_xy <- pmax(1e-10, 1 - marginal_X - marginal_Y + CDF)
      CDF <- dens / S_xy
      CDF <- pmax(0, CDF)
      ylabel <- "h(x, y)"
    } else if(type == "cumulative hazard"){
      S_xy <- pmax(1e-10, 1 - marginal_X - marginal_Y + CDF)
      CDF <- -log(S_xy)
      CDF <- pmax(0, CDF)
      ylabel <- "H(x, y)"
    }
    
    if(!is.null(target)){
      P <- (
        Co.LPM(degree, overall_target_1, overall_target_2, target[1], target[2]) /
          (
            Co.LPM(degree, overall_target_1, overall_target_2, target[1], target[2]) +
              Co.UPM(degree, overall_target_1, overall_target_2, target[1], target[2]) +
              D.LPM(degree, degree, overall_target_1, overall_target_2, target[1], target[2]) +
              D.UPM(degree, degree, overall_target_1, overall_target_2, target[1], target[2])
          )
      )
      P_marginal_X <- LPM.ratio(degree, target[1], overall_target_1)
      P_marginal_Y <- LPM.ratio(degree, target[2], overall_target_2)
      if(type == "survival"){
        P <- 1 - P_marginal_X - P_marginal_Y + P
        P <- max(0, min(1, P))
      } else if(type == "hazard"){
        P_dens <- pmax(reg_fit$Point.est, 1e-10)
        S_target <- max(1e-10, 1 - P_marginal_X - P_marginal_Y + P)
        P <- P_dens / S_target
        P <- min(P, 1e6)
        P <- max(0, P)
      } else if(type == "cumulative hazard"){
        S_target <- max(1e-10, 1 - P_marginal_X - P_marginal_Y + P)
        P <- -log(S_target)
        P <- max(0, P)
      }
    } else {
      P <- NULL
    }
    
    if(plot){
      plot3d(
        variable[sorted_indices, 1], variable[sorted_indices, 2], CDF, col = "steelblue",
        xlab = deparse(substitute(variable[,1])), ylab = deparse(substitute(variable[,2])),
        zlab = ylabel, box = FALSE, pch = 19
      )
      if(!is.null(target) && !is.na(P)){
        points3d(target[1], target[2], P, col = "green", pch = 19)
        points3d(target[1], target[2], 0, col = "red", pch = 15, cex = 2)
        lines3d(
          x = c(target[1], max(variable[,1])),
          y = c(target[2], max(variable[,2])),
          z = c(P, P),
          col = "red", lwd = 2, lty = 3
        )
        lines3d(
          x = c(target[1], target[1]),
          y = c(target[2], target[2]),
          z = c(0, P),
          col = "red", lwd = 1, lty = 3
        )
        text3d(
          max(variable[,1]), max(variable[,2]), P, texts = paste0("P = ", round(P, 4)), pos = 4, col = "red"
        )
      }
    }
    
    return(list(
      "Function" = data.table::data.table(cbind(
        data.frame(variable[sorted_indices, ], row.names = NULL), 
        CDF = CDF
      )),
      "target.value" = P
    ))
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
