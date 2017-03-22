#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return LPM of variable
#' @keywords partial moments, mean, variance, CDF
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' LPM(0,mean(x),x)
#' @export

LPM<- function(degree,target,variable){
sum((target - (variable[variable <= target]))^degree)/length(variable)}
LPM<- Vectorize(LPM,vectorize.args = 'target')



#' Upper Partial Moment
#'
#' This function generates a univariate upper partial moment for any degree or target.
#' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
#' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
#' @param variable a numeric vector.
#' @return UPM of variable
#' @keywords partial moments, mean, variance, upper CDF
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' UPM(0,mean(x),x)
#' @export


UPM<- function(degree,target,variable){
  sum(((variable[variable > target]) - target)^degree)/length(variable)}
UPM<- Vectorize(UPM,vectorize.args = 'target')

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
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' Co.UPM(0,0,x,y,mean(x),mean(y))
#' @export


Co.UPM<- function(degree.x,degree.y,x,y,target.x,target.y){
  x=x-target.x;y=y-target.y
  x[x<=0]<- 0;y[y<=0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
  }
Co.UPM<- Vectorize(Co.UPM,vectorize.args = c('target.x','target.y'))

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
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' Co.LPM(0,0,x,y,mean(x),mean(y))
#' @export

Co.LPM<- function(degree.x,degree.y,x,y,target.x,target.y){
  x=target.x-x;y=target.y-y
  x[x<0]<- 0;y[y<0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
  }
Co.LPM<- Vectorize(Co.LPM,vectorize.args = c('target.x','target.y'))

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
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' D.LPM(0,0,x,y,mean(x),mean(y))
#' @export

D.LPM<- function(degree.x,degree.y,x,y,target.x,target.y){
  x=x-target.x;y=target.y-y
  x[x<=0]<- 0;y[y<0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
  }
D.LPM<- Vectorize(D.LPM,vectorize.args = c('target.x','target.y'))

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
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' D.UPM(0,0,x,y,mean(x),mean(y))
#' @export

D.UPM<- function(degree.x,degree.y,x,y,target.x,target.y){
  x=target.x-x;y=y-target.y
  x[x<0]<- 0;y[y<=0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
 }
D.UPM<- Vectorize(D.UPM,vectorize.args = c('target.x','target.y'))
