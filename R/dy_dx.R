#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivate of y wrt x for a point of interest.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param order integer; Controls the number of partial moment quadrant means.  Defaults to \code{(order=NULL)} which generates a more accurate derivative for well specified cases.
#' @param stn numeric [0,1]; Signal to noise parameter, sets the threshold of \code{NNS.dep} which reduces \code{"order"} when \code{(order=NULL)}.  Defaults to 0.99 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param eval.point numeric; \code{x} point to be evaluated.  Defaults to \code{(eval.point=median(x))}.  Set to \code{(eval.points="overall")} to find an overall partial derivative estimate.
#' @param deriv.order numeric options: (1,2); 1 (default) For second derivative estimate of \code{f(x)}, set \code{(deriv.order=2)}.
#' @param h numeric [0,...]; Percentage step used for finite step method.  Defaults to \code{h=.01} representing a 1 percent step from the value of the independent variable.
#' @param noise.reduction the method of determing regression points options: ("mean","median","mode","off"); In low signal to noise situations, \code{(noise.reduction="median")} uses medians instead of means for partitions, while \code{(noise.reduction="mode")} uses modes instead of means for partitions.  \code{(noise.reduction="off")}  allows for maximum possible fit in \link{NNS.reg}. Default setting is \code{(noise.reduction="mean")}.
#' @param deriv.method method of derivative estimation {"NNS","FS"}; Determines the partial derivative from the coefficient of the \link{NNS.reg} output when \code{(deriv.method="NNS")} or generates a partial derivative using the finite step method \code{(deriv.method="FS")} (Defualt).
#' @return Returns the value of the partial derivative estimate for the given order.
#' @keywords partial derivative
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' x<-seq(0,2*pi,pi/100); y<-sin(x)
#' dy.dx(x,y,eval.point=1.75)
#' @export

dy.dx <- function(x,y,order=NULL,stn=0.99,eval.point=median(x),deriv.order=1,h=.01,noise.reduction='mean',deriv.method="FS"){

  if(eval.point=='overall'){

  ranges=NNS.reg(x,y,order=order,noise.reduction=noise.reduction,plot=F)$derivative
  range.weights=numeric()
  range.weights=data.table(x,'interval'=findInterval(x,ranges[,X.Lower.Range]))
  range.weights=range.weights[,.N,by='interval']
  range.weights=range.weights$N/length(x)

  return(sum(ranges[,Coefficient]*range.weights))

  } else {

  original.eval.point.min=eval.point
  original.eval.point.max=eval.point

  eval.point.min = (1-h)*original.eval.point.min
  eval.point.max = (1+h)*original.eval.point.max


  run=eval.point.max-eval.point.min

  if(deriv.order==1){

  if(deriv.method=="FS"){
  estimates=NNS.reg(x,y,plot = FALSE,order=order,stn = stn,noise.reduction = noise.reduction,point.est = c(eval.point.min,eval.point.max))$Point.est

    rise=estimates[2]-estimates[1]

    return(rise/run) } else {

     reg.output <- NNS.reg(x,y,plot = FALSE,return.values = TRUE,order=order,stn = stn,noise.reduction = noise.reduction)

     output<- reg.output$derivative
      if(length(output[,Coefficient])==1){return(output[,Coefficient])}
      if((output[,X.Upper.Range][which(eval.point<output[,X.Upper.Range])-1][1])<eval.point){
      return(output[,Coefficient][which(eval.point<output[,X.Upper.Range])][1])}
      else{
        return(mean(c(output[,Coefficient][which(eval.point<output[,X.Upper.Range])][1],output[,X.Lower.Range][which(eval.point<output[,X.Upper.Range])-1][1])))
      }
      }


  }

  else{
    ## Second derivative form:
    # f(x+h) - 2(f(x)) +f(x-h) / h^2

    deriv.points=matrix(c((1+h)*eval.point,eval.point,(1-h)*eval.point),ncol = length(eval.point),byrow = TRUE)
    second.deriv.estimates= NNS.reg(x,y,plot = FALSE,return.values = TRUE,order=order,point.est = deriv.points,stn = stn,noise.reduction = noise.reduction)$Point.est
    f.x_h = second.deriv.estimates[1]

    two_f.x = 2*second.deriv.estimates[2]

    f.x__h = second.deriv.estimates[3]

    run = ((1+h)*eval.point) - eval.point
  return((f.x_h - two_f.x + f.x__h)/(run^2))

  }

    }

}
