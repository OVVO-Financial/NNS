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
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"
#' \url{https://www.ssrn.com/abstract=3007373}

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
  
  if(any(class(variable)%in%c("tbl","data.table")) && dim(variable)[2]==1){ 
    variable <- as.vector(unlist(variable))
  }
  if(any(class(variable)%in%c("tbl","data.table"))){
    variable <- as.data.frame(variable)
  }
  
  if(!is.null(target)){
    if(is.null(dim(variable)) || dim(variable)[2]==1){
      if(target<min(variable) || target>max(variable)){
        stop("Please make sure target is within the observed values of variable.")
      }
    } else {
      if(target[1]<min(variable[,1]) || target[1]>max(variable[,1])){
        stop("Please make sure target 1 is within the observed values of variable 1.")
      }
      if(target[2]<min(variable[,2]) || target[2]>max(variable[,2])){
        stop("Please make sure target 2 is within the observed values of variable 2.")
      }
    }
  }
  type <- tolower(type)
  if(!(type%in%c("cdf","survival", "hazard", "cumulative hazard"))){
    stop(paste("Please select a type from: ", "`CDF`, ", "`survival`, ",  "`hazard`, ", "`cumulative hazard`"))
  }
  if(is.null(dim(variable)) || dim(variable)[2] == 1){
    overall_target <- sort(variable)
    x <- overall_target
    if(degree > 0){
      CDF <- LPM.ratio(degree, overall_target, variable)
    } else {
      cdf_fun <- ecdf(x)
      CDF <- cdf_fun(overall_target)
    }
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
    }else if(type == "hazard"){
      CDF <- exp(log(density(x, n = length(x))$y)-log(1-CDF))
      ylabel <- "h(x)"
      P <- NNS.reg(x[-length(x)], CDF[-length(x)], order = "max", point.est = c(x[length(x)], target), plot = FALSE)$Point.est
      CDF[is.infinite(CDF)] <- P[1]
      P <- P[-1]
    }else if(type == "cumulative hazard"){
      CDF <- -log((1 - CDF))
      ylabel <- "H(x)"
      P <- NNS.reg(x[-length(x)], CDF[-length(x)], order = "max", point.est = c(x[length(x)], target), plot = FALSE)$Point.est
      CDF[is.infinite(CDF)] <- P[1]
      P <- P[-1]
    }
    if(plot){
      plot(x, CDF, pch = 19, col = 'steelblue', xlab = deparse(substitute(variable)), ylab = ylabel, main = toupper(type), type = "s", lwd = 2)
      points(x, CDF, pch = 19, col = 'steelblue')
      lines(x, CDF, lty=2, col = 'steelblue')
      if(!is.null(target)){
        segments(target,0,target,P, col = "red", lwd = 2, lty = 2)
        segments(min(variable), P, target, P, col = "red", lwd = 2, lty = 2)
        points(target, P, col = "green", pch = 19)
        mtext(text = round(P,4), col = "red", side = 2, at = P,  las = 2)
        mtext(text = round(target,4), col = "red", side = 1, at = target,  las = 1)
      }
    }
    values <- data.table::data.table(cbind.data.frame(x, CDF))
    colnames(values) <- c(deparse(substitute(variable)), ylabel)
    return(
      list(
        "Function" = values ,
        "target.value" = P
      )
    )
  } else {
    overall_target_1 <- (variable[,1])
    overall_target_2 <- (variable[,2])
    CDF <- (
      Co.LPM(degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) /
      (
        Co.LPM(degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) +
          Co.UPM(degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) +
          D.UPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) +
          D.LPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2)
      )
    )
    if(type == "survival"){
      CDF <- 1 - CDF
    } else if(type == "hazard"){
      CDF <- sort(variable) / (1 - CDF)
    } else if(type == "cumulative hazard"){
      CDF <- -log((1 - CDF))
    }
    if(!is.null(target)){
      P <- (
        Co.LPM(degree, variable[,1], variable[,2], target[1], target[2]) /
        (
          Co.LPM(degree, variable[,1], variable[,2], target[1], target[2]) +
            Co.UPM(degree, variable[,1], variable[,2], target[1], target[2]) +
            D.LPM(degree,degree, variable[,1], variable[,2], target[1], target[2]) +
            D.UPM(degree,degree, variable[,1], variable[,2], target[1], target[2])
        )
      )
    } else {
      P <- NULL
    }
    if(plot){
      plot3d(
        variable[,1], variable[,2], CDF, col = "steelblue",
        xlab = deparse(substitute(variable[,1])), ylab = deparse(substitute(variable[,2])),
        zlab = "Probability", box = FALSE, pch = 19
      )
      if(!is.null(target)){
        points3d(target[1], target[2], P, col = "green", pch = 19)
        points3d(target[1], target[2], 0, col = "red", pch = 15, cex = 2)
        lines3d(
          x= c(target[1], max(variable[,1])),
          y= c(target[2], max(variable[,2])),
          z= c(P, P),
          col = "red", lwd = 2, lty=3
        )
        lines3d(
          x= c(target[1], target[1]),
          y= c(target[2], target[2]),
          z= c(0, P),
          col = "red", lwd = 1, lty=3
        )
        text3d(
          max(variable[,1]), max(variable[,2]), P, texts = paste0("P = ", round(P,4)), pos = 4, col = "red"
        )
      }
      
    }
    
  }
  
  return(list("CDF" = data.table::data.table(cbind((variable), CDF = CDF)),
              "P" = P))
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
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
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
    variance <- variance * (n / (n - 1))
  } else {
    skewness <- (n / ((n-1)*(n-2))) * ((n*skew_base) / variance^(3/2))
    kurtosis <- ((n * (n+1)) / ((n-1)*(n-2)*(n-3))) * ((n*kurt_base) / (variance * (n / (n - 1)))^2) - ( (3 * ((n-1)^2)) / ((n-2)*(n-3)))
  }

  return(list("mean" = mean,
              "variance" = variance,
              "skewness" = skewness,
              "kurtosis" = kurtosis))
}
