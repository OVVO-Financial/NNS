#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#'
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return LPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' LPM(0, mean(x), x)
#' @export

LPM <-  function(degree, target, variable){
    if(degree == 0) return(mean(variable <= target))

    sum((target - (variable[variable <= target])) ^ degree) / length(variable)
}
LPM <- Vectorize(LPM, vectorize.args = 'target')


#' Upper Partial Moment
#'
#' This function generates a univariate upper partial moment for any degree or target.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return UPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' UPM(0, mean(x), x)
#' @export


UPM <-  function(degree, target, variable){
  if(degree == 0) return(mean(variable > target))

  sum(((variable[variable > target]) - target) ^ degree) / length(variable)
}
UPM <- Vectorize(UPM, vectorize.args = 'target')

#' Co-Upper Partial Moment
#' (Upper Right Quadrant 1)
#'
#' This function generates a co-upper partial moment between two equal length variables for any degree or target.
#' @param degree.x integer; Degree for variable X.  \code{(degree.x = 0)} is frequency, \code{(degree.x = 1)} is area.
#' @param degree.y integer; Degree for variable Y.  \code{(degree.y = 0)} is frequency, \code{(degree.y = 1)} is area.
#' @param x a numeric vector.
#' @param y a numeric vector of equal length to \code{x}.
#' @param target.x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)
#' @param target.y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be. (Vectorized)
#' @return Co-UPM of two variables
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' Co.UPM(0, 0, x, y, mean(x), mean(y))
#' @export


Co.UPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  z <- cbind(x, y)
  z <- t(t(z) - c(target.x, target.y))
  z[z<=0] <- NA
  z <- z[complete.cases(z), , drop = FALSE]
  return(z[,1]^degree.x %*% z[,2]^degree.y / length(x))
  }
Co.UPM <- Vectorize(Co.UPM, vectorize.args = c('target.x', 'target.y'))

#' Co-Lower Partial Moment
#' (Lower Left Quadrant 4)
#'
#' This function generates a co-lower partial moment for between two equal length variables for any degree or target.
#' @param degree.x integer; Degree for variable X.  \code{(degree.x = 0)} is frequency, \code{(degree.x = 1)} is area.
#' @param degree.y integer; Degree for variable Y.  \code{(degree.y = 0)} is frequency, \code{(degree.y = 1)} is area.
#' @param x a numeric vector.
#' @param y a numeric vector of equal length to \code{x}.
#' @param target.x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)
#' @param target.y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be. (Vectorized)
#' @return Co-LPM of two variables
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' Co.LPM(0, 0, x, y, mean(x), mean(y))
#' @export

Co.LPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  z <- cbind(x,y)
  z <- t(c(target.x, target.y) - t(z))
  z[z<=0] <- NA
  z <- z[complete.cases(z), , drop = FALSE]
  return(z[,1]^degree.x %*% z[,2]^degree.y / length(x))
}
Co.LPM <- Vectorize(Co.LPM, vectorize.args = c('target.x', 'target.y'))

#' Divergent-Lower Partial Moment
#' (Lower Right Quadrant 3)
#'
#' This function generates a divergent lower partial moment between two equal length variables for any degree or target.
#' @param degree.x integer; Degree for variable X.  \code{(degree.x = 0)} is frequency, \code{(degree.x = 1)} is area.
#' @param degree.y integer; Degree for variable Y.  \code{(degree.y = 0)} is frequency, \code{(degree.y = 1)} is area.
#' @param x a numeric vector.
#' @param y a numeric vector of equal length to \code{x}.
#' @param target.x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)
#' @param target.y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be. (Vectorized)
#' @return Divergent LPM of two variables
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' D.LPM(0, 0, x, y, mean(x), mean(y))
#' @export

D.LPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  z <- cbind(x,y)
  z[,1] <- z[,1] - target.x
  z[,2] <- target.y - z[,2]
  z[z<=0] <- NA
  z <- z[complete.cases(z), , drop = FALSE]
  return(z[,1]^degree.x %*% z[,2]^degree.y / length(x))
}
D.LPM <- Vectorize(D.LPM, vectorize.args = c('target.x', 'target.y'))

#' Divergent-Upper Partial Moment
#' (Upper Left Quadrant 2)
#'
#' This function generates a divergent upper partial moment between two equal length variables for any degree or target.
#' @param degree.x integer; Degree for variable X.  \code{(degree.x = 0)} is frequency, \code{(degree.x = 1)} is area.
#' @param degree.y integer; Degree for variable Y.  \code{(degree.y = 0)} is frequency, \code{(degree.y = 1)} is area.
#' @param x a numeric vector.
#' @param y a numeric vector of equal length to \code{x}.
#' @param target.x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)
#' @param target.y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be. (Vectorized)
#' @return Divergent UPM of two variables
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' D.UPM(0, 0, x, y, mean(x), mean(y))
#' @export

D.UPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  z <- cbind(x,y)
  z[,1] <- target.x - z[,1]
  z[,2] <- z[,2] - target.y
  z[z<=0] <- NA
  z <- z[complete.cases(z), , drop = FALSE]
  return(z[,1]^degree.x %*% z[,2]^degree.y / length(x))
 }
D.UPM <- Vectorize(D.UPM, vectorize.args = c('target.x', 'target.y'))


#' Partial Moment Matrix
#'
#'
#' This function generates a co-partial moment matrix for the specified co-partial moment.
#' @param LPM.degree integer; Degree for \code{variable} below \code{target} deviations.  \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param UPM.degree integer; Degree for \code{variable} above \code{target} deviations.  \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)  \code{(target = "mean")} (default) will set the target as the mean of every variable.
#' @param variable a numeric matrix or data.frame.
#' @param pop.adj logical; \code{FALSE} (default) Adjusts the sample co-partial moment matrices for population statistics.
#' @return Matrix of partial moment quadrant values (CUPM, DUPM, DLPM, CLPM), and overall covariance matrix.  Uncalled quadrants will return a matrix of zeros.
#' @note For divergent asymmetical \code{"D.LPM" and "D.UPM"} matrices, matrix is \code{D.LPM(column,row,...)}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @references Viole, F. (2017) "Bayes' Theorem From Partial Moments"
#' \url{https://www.ssrn.com/abstract=3457377}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' A <- cbind(x,y,z)
#' PM.matrix(LPM.degree = 1, UPM.degree = 1, target = "mean", variable = A)
#'
#' ## Use of vectorized numeric targets (target_x, target_y, target_z)
#' PM.matrix(LPM.degree = 1, UPM.degree = 1, target = c(0, 0.15, .25), variable = A)
#'
#' ## Calling Individual Partial Moment Quadrants
#' cov.mtx <- PM.matrix(LPM.degree = 1, UPM.degree = 1, target = "mean", variable = A)
#' cov.mtx$cupm
#'
#' ## Full covariance matrix
#' cov.mtx$cov.matrix
#' @export


PM.matrix <- function(LPM.degree, UPM.degree, target = "mean", variable, pop.adj=FALSE){

  n <- ncol(variable)
  if(is.null(n)){stop("supply a matrix-like 'variable'")}

    clpms <- list()
    cupms <- list()
    dlpms <- list()
    dupms <- list()

    for(i in 1 : n){
        if(is.numeric(target)){
            clpms[[i]] <- sapply(1 : n, function(b) Co.LPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = LPM.degree, target.x = target[i], target.y = target[b]))

            cupms[[i]] <- sapply(1 : n, function(b) Co.UPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = UPM.degree, target.x = target[i], target.y = target[b]))

            dlpms[[i]] <- sapply(1 : n, function(b) D.LPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = LPM.degree, target.x = target[i], target.y = target[b]))

            dupms[[i]] <- sapply(1 : n, function(b) D.UPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = UPM.degree, target.x = target[i], target.y = target[b]))

        } else {
            clpms[[i]] <- sapply(1 : n, function(b) Co.LPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = LPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))

            cupms[[i]] <- sapply(1 : n, function(b) Co.UPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = UPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))

            dlpms[[i]] <- sapply(1 : n, function(b) D.LPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = LPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))

            dupms[[i]] <- sapply(1 : n, function(b) D.UPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = UPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))
        }
    }

    clpm.matrix <- matrix(unlist(clpms), n, n)
    colnames(clpm.matrix) <- colnames(variable)
    rownames(clpm.matrix) <- colnames(variable)

    cupm.matrix <- matrix(unlist(cupms), n, n)
    colnames(cupm.matrix) <- colnames(variable)
    rownames(cupm.matrix) <- colnames(variable)

    dlpm.matrix <- matrix(unlist(dlpms), n, n)
    diag(dlpm.matrix) <- 0
    colnames(dlpm.matrix) <- colnames(variable)
    rownames(dlpm.matrix) <- colnames(variable)

    dupm.matrix <- matrix(unlist(dupms), n, n)
    diag(dupm.matrix) <- 0
    colnames(dupm.matrix) <- colnames(variable)
    rownames(dupm.matrix) <- colnames(variable)


  if(pop.adj){
    adjustment <- length(variable[ , 1]) / (length(variable[ , 1]) - 1)
    clpm.matrix <- clpm.matrix*adjustment
    cupm.matrix <- cupm.matrix*adjustment
    dlpm.matrix <- dlpm.matrix*adjustment
    dupm.matrix <- dupm.matrix*adjustment
  }

  cov.matrix <- cupm.matrix + clpm.matrix - dupm.matrix - dlpm.matrix

  return(list(cupm = cupm.matrix,
              dupm = dupm.matrix,
              dlpm = dlpm.matrix,
              clpm = clpm.matrix,
              cov.matrix = cov.matrix
              ))
}


#' Lower Partial Moment RATIO
#'
#' This function generates a standardized univariate lower partial moment for any degree or target.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return Standardized LPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @references Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"
#' \url{https://www.ssrn.com/abstract=3007373}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' LPM.ratio(0, mean(x), x)
#'
#' \dontrun{
#' ## Empirical CDF (degree = 0)
#' lpm_cdf <- LPM.ratio(0, sort(x), x)
#' plot(sort(x), lpm_cdf)
#'
#' ## Continuous CDF (degree = 1)
#' lpm_cdf_1 <- LPM.ratio(1, sort(x), x)
#' plot(sort(x), lpm_cdf_1)
#'
#' ## Joint CDF
#' x <- rnorm(5000) ; y <- rnorm(5000)
#' plot3d(x, y, Co.LPM(0, 0, sort(x), sort(y), x, y), col = "blue", xlab = "X", ylab = "Y",
#' zlab = "Probability", box = FALSE)
#' }
#' @export

LPM.ratio <- function(degree, target, variable){
  lpm <- LPM(degree, target, variable)

  if(degree>0){
      area <- lpm + UPM(degree, target, variable)
  } else {
      area <- 1
  }



  return(lpm / area)
}



#' Upper Partial Moment RATIO
#'
#' This function generates a standardized univariate upper partial moment for any degree or target.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return Standardized UPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' UPM.ratio(0, mean(x), x)
#'
#' ## Joint Upper CDF
#' \dontrun{
#' x <- rnorm(5000) ; y <- rnorm(5000)
#' plot3d(x, y, Co.UPM(0, 0, sort(x), sort(y), x, y), col = "blue", xlab = "X", ylab = "Y",
#' zlab = "Probability", box = FALSE)
#' }
#' @export


UPM.ratio <- function(degree, target, variable){
  upm <- UPM(degree, target, variable)

  if(degree>0){
    area <- LPM(degree, target, variable) + upm
  } else {
    area <- 1
  }

  return(upm / area)
}



#' NNS PDF
#'
#' This function generates an empirical PDF using \link{dy.dx} on \link{NNS.CDF}.
#'
#' @param variable a numeric vector.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} (default) is area.
#' @param target a numeric range of values [a,b] where a < b.  \code{NULL} (default) uses the \code{variable} min and max observations respectively.
#' @param bins integer; \code{NULL} Selects number of bins.  Bin width defaults to \code{density(x)$bw}.
#' @param plot logical; plots PDF.
#' @return Returns a data.table containing the intervals used and resulting PDF of the variable.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.PDF(x)
#'
#' ## Custom target range
#' NNS.PDF(x, target = c(-5, 5))
#' @export


NNS.PDF <- function(variable, degree = 1, target = NULL, bins = NULL , plot = TRUE){
  if(is.null(target)){target <- sort(variable)}

# d/dx approximation
  if(is.null(bins)){
      bins <- density(variable)$bw
      tgt <- seq(min(target), max(target), bins)
  } else {
      d.dx <- (abs(max(target)) + abs(min(target))) / bins
      tgt <- seq(min(target), max(target), d.dx)
  }



  CDF <- NNS.CDF(variable, plot = FALSE, degree = degree)$Function
  PDF <- pmax(dy.dx(unlist(CDF[,1]), unlist(CDF[,2]), eval.point = tgt, deriv.method = "FD")$First, 0)

  if(plot){plot(tgt, PDF, col = 'steelblue', type = 'l', lwd = 3, xlab = "X", ylab = "Density")}

  return(data.table::data.table(cbind("Intervals" = tgt, PDF)))
}


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
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.CDF(x)
#'
#' \dontrun{
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

  if(!is.null(target)){
     if(is.null(dim(variable)) || dim(variable)[2]==1){
          if(target<min(variable) || target>max(variable))   stop("Please make sure target is within the observed values of variable.")
     } else {
        if(target[1]<min(variable[,1]) || target[1]>max(variable[,1])) stop("Please make sure target 1 is within the observed values of variable 1.")
        if(target[2]<min(variable[,2]) || target[2]>max(variable[,2])) stop("Please make sure target 2 is within the observed values of variable 2.")
    }
  }

  type <- tolower(type)

  if(!(type%in%c("cdf","survival", "hazard", "cumulative hazard"))) stop(paste("Please select a type from: ", "`CDF`, ", "`survival`, ",  "`hazard`, ", "`cumulative hazard`"))

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
    }


    if(type == "hazard"){
      CDF <- exp(log(density(x, n = length(x))$y)-log(1-CDF))

      ylabel <- "h(x)"
      P <- NNS.reg(x[-length(x)], CDF[-length(x)], order = "max", point.est = c(x[length(x)], target), plot = FALSE)$Point.est
      CDF[is.infinite(CDF)] <- P[1]
      P <- P[-1]
    }


    if(type == "cumulative hazard"){
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


    return(list("Function" = values ,
                "target.value" = P))

  } else {
    overall_target_1 <- (variable[,1])
    overall_target_2 <- (variable[,2])

    CDF <- Co.LPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) /
                (
                 Co.LPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) +
                 Co.UPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) +
                 D.UPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2) +
                 D.LPM(degree,degree, sort(variable[,1]), sort(variable[,2]), overall_target_1, overall_target_2)
                )

    if(type == "survival"){
        CDF <- 1 - CDF
    }

    if(type == "hazard"){
        CDF <- sort(variable) / (1 - CDF)
    }

    if(type == "cumulative hazard"){
        CDF <- -log((1 - CDF))
    }


    if(!is.null(target)){
      P <- Co.LPM(degree,degree, variable[,1], variable[,2], target[1], target[2]) /
                (
                  Co.LPM(degree,degree, variable[,1], variable[,2], target[1], target[2]) +
                  Co.UPM(degree,degree, variable[,1], variable[,2], target[1], target[2]) +
                  D.LPM(degree,degree, variable[,1], variable[,2], target[1], target[2]) +
                  D.UPM(degree,degree, variable[,1], variable[,2], target[1], target[2])
                )


    } else {
      P <- NULL
    }

    if(plot){
      plot3d(variable[,1], variable[,2], CDF, col = "steelblue",
             xlab = deparse(substitute(variable[,1])), ylab = deparse(substitute(variable[,2])),
             zlab = "Probability", box = FALSE, pch = 19)

      if(!is.null(target)){
          points3d(target[1], target[2], P, col = "green", pch = 19)
          points3d(target[1], target[2], 0, col = "red", pch = 15, cex = 2)
          lines3d(x= c(target[1], max(variable[,1])),
                  y= c(target[2], max(variable[,2])),
                  z= c(P, P),
                  col = "red", lwd = 2, lty=3)
          lines3d(x= c(target[1], target[1]),
                  y= c(target[2], target[2]),
                  z= c(0, P),
                  col = "red", lwd = 1, lty=3)
          text3d(max(variable[,1]), max(variable[,2]), P, texts = paste0("P = ", round(P,4)), pos = 4, col = "red")
      }

    }

  }

  return(list("CDF" = data.table::data.table(cbind((variable), CDF = CDF)),
              "P" = P))


}


