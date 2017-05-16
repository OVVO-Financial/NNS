#' NNS Correlation
#'
#' Returns the nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param order integer; Controls the level of quadrant partitioning.  Defualts to \code{(order=NULL)}.  Errors can generally be rectified by setting \code{(order=1)}.
#' @param degree integer; \code{(degree = 0)} is frequency based correlations, while \code{(degree = 1)} is for area based correlations.  Defaults to \code{(degree = 0)} for smaller number of observations.
#' @return Returns nonlinear correlation coefficient between two variables, or nonlinear correlation matrix for matrix input.
#' @keywords nonlinear correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' ## Pairwise Correlation
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.cor(x,y)
#'
#' ## Correlation Matrix
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' B<-cbind(x,y,z)
#' NNS.cor(B)
#'
#' @export

NNS.cor = function(x, y=NULL, order = NULL,
                   degree= NULL){

  if(is.null(y)){
  NNS.dep(x,order=order,degree=degree)$Correlation}
  else{
    NNS.dep(x,y,order=order,degree=degree)$Correlation}


}
