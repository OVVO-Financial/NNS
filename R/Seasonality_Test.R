#' NNS Seasonality Test
#'
#' Seasonality test based on the coefficient of variance for the variable and lagged component series.  A result of 1 signifies no seasonality present.
#'
#' @param variable a numeric vector.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference.
#' @return Returns a matrix of all periods exhibiting less coefficient of variance than the variable with \code{"all.periods"}; and the single period exhibiting the least coefficient of variance versus the variable with \code{"best.period"}.  If no seasonality is detected, \code{NNS.seas} will return (1).
#' @keywords seasonality
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#'
#' ## To call strongest period based on coefficient of variance:
#' NNS.seas(x)$best.period
#' @export


NNS.seas <- function(variable,plot=TRUE){

  output <- vector("numeric", length(variable)/4)
  instances <- vector("numeric", length(variable)/4)

  for (i in 1:(length(variable)/4)){
      if (abs(sd(variable[seq(length(variable),1,-i)])/mean(variable[seq(length(variable),1,-i)])) < abs(sd(variable)/mean(variable))){

          instances[i] <- i

          output[i]<- (abs(sd(variable[seq(length(variable),1,-i)])/mean(variable[seq(length(variable),1,-i)])))
      } else {
            instances[i] <- 0
            output[i]<- 0
        }
  }


  if(length(instances[instances>0])>0){

      n<- rep(abs(sd(variable)/mean(variable)),length(instances[instances>0]))

      M<- data.table("Period"=instances[instances>0],"Coefficient.of.Variance"=output[output>0],"Variable.Coefficient.of.Variance"=n,key = "Coefficient.of.Variance")


    }

  if(plot==T){
    if(sum(instances[instances>0])>0){
  plot(instances[instances>0],output[output>0],
         xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",
         ylim = c(0,2*abs(sd(variable)/mean(variable))))

    points(M[1,Period],M[1,Coefficient.of.Variance],pch=19,col='red')

    abline(h=abs(sd(variable)/mean(variable)), col="red",lty=5)
    text(mean(instances[instances>0]),abs(sd(variable)/mean(variable)),pos=3,
         "Variable Coefficient of Variance",col='red')}

           else {
            plot(1,1,pch=19,col='blue', xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",         ylim = c(0,2*abs(sd(variable)/mean(variable))))

             text(1,abs(sd(variable)/mean(variable)),pos=3,"NO SEASONALITY DETECTED",col='red')
           }

  }

  if(sum(instances[instances>0])==0) {return('best.period'=t(1))}

if(length(instances[instances>0])>0){

  return(list("all.periods"=M,"best.period"=M[1,Period]))
}
    else {
     list("all.periods"=NULL, "best.period"=1)}


}
