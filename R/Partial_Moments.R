#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return LPM of variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' LPM(0, mean(x), x)
#' @export

LPM <-  function(degree, target, variable){
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
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' UPM(0, mean(x), x)
#' @export


UPM <-  function(degree, target, variable){
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
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' Co.UPM(0,0,x,y,mean(x),mean(y))
#' @export


Co.UPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  if(degree.x == 0){x[x == target.x] <- target.x - 1}
  if(degree.y == 0){y[y == target.y] <- target.y - 1}
  z <- cbind(x,y); z <- z[complete.cases(z),]
  x <- z[,1]
  y <- z[,2]
  x <- x - target.x
  y <- y - target.y
  x[x <= 0] <- 0
  y[y <= 0] <- 0
  x[x > 0] <- x[x > 0] ^ degree.x
  y[y > 0] <- y[y > 0] ^ degree.y
  return(x %*% y / length(x))
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
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' Co.LPM(0, 0, x, y, mean(x), mean(y))
#' @export

Co.LPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  if(degree.x == 0){x[x == target.x] <- target.x - 1}
  if(degree.y == 0){y[y == target.y] <- target.y - 1}
  z <- cbind(x,y); z <- z[complete.cases(z),]
  x <- z[,1]
  y <- z[,2]
  x <- target.x - x
  y <- target.y - y
  x[x <= 0] <- 0
  y[y <= 0] <- 0
  x[x > 0] <- x[x > 0] ^ degree.x
  y[y > 0] <- y[y > 0] ^ degree.y
  return(x %*% y / length(x))
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
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' D.LPM(0, 0, x, y, mean(x), mean(y))
#' @export

D.LPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  if(degree.x == 0){x[x == target.x] <- target.x - 1}
  if(degree.y == 0){y[y == target.y] <- target.y - 1}
  z <- cbind(x,y); z <- z[complete.cases(z),]
  x <- z[,1]
  y <- z[,2]
  x <- x - target.x
  y <- target.y - y
  x[x <= 0] <- 0
  y[y <= 0] <- 0
  x[x > 0] <- x[x > 0] ^ degree.x
  y[y > 0] <- y[y > 0] ^ degree.y
  return(x %*% y / length(x))
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
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' D.UPM(0, 0, x, y, mean(x), mean(y))
#' @export

D.UPM <- function(degree.x, degree.y, x, y, target.x = mean(x), target.y = mean(y)){
  if(degree.x == 0){x[x == target.x] <- target.x - 1}
  if(degree.y == 0){y[y == target.y] <- target.y - 1}
  z <- cbind(x,y); z <- z[complete.cases(z),]
  x <- z[,1]
  y <- z[,2]
  x <- target.x - x
  y <- y - target.y
  x[x <= 0] <- 0
  y[y <= 0] <- 0
  x[x > 0] <- x[x > 0] ^ degree.x
  y[y > 0] <- y[y > 0] ^ degree.y
  return(x %*% y / length(x))
 }
D.UPM <- Vectorize(D.UPM, vectorize.args = c('target.x', 'target.y'))


#' Partial Moment Matrix
#'
#'
#' This function generates a co-partial moment matrix for the specified co-partial moment.
#' @param LPM.degree integer; Degree for \code{variable} below \code{target} deviations.  \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param UPM.degree integer; Degree for \code{variable} above \code{target} deviations.  \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)  \code{(target = "mean")} will set the target as the mean of every variable.
#' @param variable a numeric matrix or data.frame.
#' @param pop.adj logical; \code{FALSE} (default) Adjusts the sample co-partial moment matrices for population statistics.
#' @return Matrix of partial moment quadrant values.  Uncalled quadrants will return a matrix of zeros.
#' @note For divergent asymmetical \code{"D.LPM" and "D.UPM"} matrices, matrix is \code{D.LPM(column,row,...)}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' A <- cbind(x,y,z)
#' PM.matrix(LPM.degree = 1, UPM.degree = 1, target = "mean", variable = A)
#'
#' ## Calling Individual Partial Moment Quadrants
#' cov.mtx <- PM.matrix(LPM.degree = 1, UPM.degree = 1, target = "mean", variable = A)
#' cov.mtx$cupm
#'
#' ## Full covariance matrix
#' cov.mtx$cov.matrix
#' @export


PM.matrix <- function(LPM.degree, UPM.degree, target, variable, pop.adj=FALSE){

  n <- ncol(variable)
  if(is.null(n)){stop("supply a matrix-like 'variable'")}

    clpms <- list()

    for(i in 1 : n){
        if(is.numeric(target)){
            clpms[[i]] <- sapply(1 : n, function(b) Co.LPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = LPM.degree, target.x = target, target.y = target))
        } else {
            clpms[[i]] <- sapply(1 : n, function(b) Co.LPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = LPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))
        }
    }

    clpm.matrix <- matrix(unlist(clpms), n, n)
    colnames(clpm.matrix) <- colnames(variable)
    rownames(clpm.matrix) <- colnames(variable)


    cupms <- list()

    for(i in 1 : n){
        if(is.numeric(target)){
            cupms[[i]] <- sapply(1 : n, function(b) Co.UPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = UPM.degree, target.x = target, target.y = target))
        } else {
            cupms[[i]] <- sapply(1 : n, function(b) Co.UPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = UPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))
        }
    }

    cupm.matrix <- matrix(unlist(cupms), n, n)
    colnames(cupm.matrix) <- colnames(variable)
    rownames(cupm.matrix) <- colnames(variable)


    dlpms <- list()

    for(i in 1 : n){
        if(is.numeric(target)){
            dlpms[[i]] <- sapply(1 : n, function(b) D.LPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = LPM.degree, target.x = target, target.y = target))
        } else {
            dlpms[[i]] <- sapply(1 : n, function(b) D.LPM(x = variable[ , i], y = variable[ , b], degree.x = UPM.degree, degree.y = LPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))
        }
    }

    dlpm.matrix <- matrix(unlist(dlpms), n, n)
    diag(dlpm.matrix) <- 0
    colnames(dlpm.matrix) <- colnames(variable)
    rownames(dlpm.matrix) <- colnames(variable)


    dupms <- list()

    for(i in 1 : n){
        if(is.numeric(target)){
            dupms[[i]] <- sapply(1 : n, function(b) D.UPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = UPM.degree, target.x = target, target.y = target))
        } else {
            dupms[[i]] <- sapply(1 : n, function(b) D.UPM(x = variable[ , i], y = variable[ , b], degree.x = LPM.degree, degree.y = UPM.degree, target.x = mean(variable[ , i]), target.y = mean(variable[ , b])))
        }
    }

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

  return(list(cov.matrix = cov.matrix,
              cupm = cupm.matrix,
              dupm = dupm.matrix,
              dlpm = dlpm.matrix,
              clpm = clpm.matrix
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
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' LPM.ratio(0, mean(x), x)
#'
#' ## Continuous CDF
#' LPM.ratio(1, sort(x), x)
#'
#'
#' ## Joint CDF
#' \dontrun{
#' x <- rnorm(5000) ; y <- rnorm(5000)
#' plot3d(x, y, Co.LPM(0, 0, sort(x), sort(y), x, y), col = "blue", xlab = "X", ylab = "Y",
#' zlab = "Probability", box = FALSE)
#' }
#' @export

LPM.ratio <- function(degree, target, variable){
  lpm <- LPM(degree, target, variable)
  upm <- UPM(degree, target, variable)

  lpm / (lpm + upm)
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
#' \url{http://amzn.com/1490523995}
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
  lpm <- LPM(degree, target, variable)
  upm <- UPM(degree, target, variable)

  upm / (lpm + upm)
}


#' NNS PDF
#'
#' This function generates an empirical PDF using continuous CDFs from \link{LPM.ratio}.
#'
#' @param variable a numeric vector.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} (default) is area.
#' @param target a numeric range of values [a,b] where a < b.  \code{NULL} (default) uses the \code{variable} observations.
#' @param bins numeric; \code{NULL} (default) Selects number of observations as default bins.
#' @param plot logical; plots PDF.
#' @return Returns a data.table containing the intervals used and resulting PDF of the variable.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.PDF(x)
#'
#' ## Custom target range
#' NNS.PDF(x, target = c(-5, 5))
#' @export


NNS.PDF <- function(variable, degree = 1, target = NULL, bins = NULL, plot = TRUE){

  if(is.null(target)){target <- sort(variable)}

# d/dx approximation
  if(is.null(bins)){bins <- length(variable)}

  d.dx <- (abs(max(target)) + abs(min(target))) / bins
  tgt <- seq(min(target), max(target), d.dx)
  PDF <- abs((diff(LPM.ratio(degree, tgt, variable),2)))

  Intervals <- (sort(tgt)+(d.dx/2))[1:length(PDF)]

  if(plot){plot(Intervals, PDF, col = 'steelblue', type = 'l', lwd = 3, xlab = "X")}

  return(data.table(cbind("Intervals" = Intervals, PDF)))
}
