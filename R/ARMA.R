#' NNS ARMA
#'
#' Autoregressive model incorporating nonlinear regressions of component series.
#'
#' @param variable a numeric vector.
#' @param h integer; 1 (default) Number of periods to forecast.
#' @param training.set \code{NULL} (defualt) numeric; Sets the number of variable observations \code{(variable[1:training.set])} to monitor performance of forecast over in-sample range.
#' @param seasonal.factor logical or integer; \code{TRUE} (default) Automatically selects the best seasonal lag from the seasonality test.  To use weighted average of all seasonal lags set to \code{(seasonal.factor=FALSE)}.  Otherwise, directly input known frequency integer lag to use, i.e. \code{(seasonal.factor=12)} for monthly data.
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to \code{(negative.values=TRUE)}.
#' @param method options:("lin","nonlin","both"); \code{"nonlin"} (default)  To select the regression type of the component series, select \code{(method="both")} where both linear and nonlinear estimates are generated.  To use a nonlineaer regression, set to \code{(method="nonlin")}; to use a linear regression set to \code{(method="lin")}.
#' @param dynamic logical; \code{FALSE} (default) To update the seasonal factor with each forecast point, set to \code{(dynamic=TRUE)}.  The default is \code{(dynamic=FALSE)} to retain the original seasonal factor from the inputted variable for all ensuing \code{h}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the \code{variable} level reference in upper panel.  Lower panel returns original data and forecast.
#' @param seasonal.plot logical; \code{TRUE} (default) Adds the seasonality plot above the forecast.  Will be set to \code{FALSE} if no seasonality is detected or \code{seasonal.factor} is set to an integer value.
#' @param intervals logical; \code{FALSE} (default) Plots the surrounding forecasts around the final estimate when \code{(intervals=TRUE)} and \code{(seasonal.factor=FALSE)}.  There are no other forecasts to plot when a single \code{seasonal.factor} is selected.
#' @return Returns a vector of forecasts of length \code{(h)}.
#' @note \code{(seasonal.factor=FALSE)} can be a very comutationally expensive exercise due to the number of seasonal periods detected.
#' @keywords Autoregressive model
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data and 12 period lag
#' NNS.ARMA(AirPassengers,h=45,training.set=100,seasonal.factor=12,method='nonlin')
#' @export



# Autoregressive Model
NNS.ARMA <- function(variable,h=1,training.set = NULL, seasonal.factor = TRUE ,negative.values = FALSE, method = "nonlin", dynamic = FALSE,plot=TRUE,seasonal.plot=TRUE,intervals=FALSE){

  if(intervals==TRUE && is.numeric(seasonal.factor)){stop('Hmmm...Seems you have "intervals" and "seasonal.factor" selected.  Please set "intervals=F" or "seasonal.factor=F"')}

  if(intervals==TRUE && seasonal.factor==TRUE){stop('Hmmm...Seems you have "intervals" and "seasonal.factor" selected.  Please set "intervals=F" or "seasonal.factor=F"')}


  variable=as.numeric(variable)
  OV = variable

  if(!is.null(training.set)){
    variable = variable[1:training.set]
    FV = variable[1:training.set]
  } else {
    training.set = length(variable)
    variable = variable
    FV = variable
  }

  Estimates = numeric()

  # Weight and lag function for seasonality...
  ARMA.seas.weighting=function(){
    if(is.null(ncol(M))){
      return(list(lag=M[1],Weights=1))}

    if(ncol(M)==1){
      return(list(lag=1,Weights=1))}

    if(ncol(M)>1){
        if(seasonal.factor==TRUE){
            lag = seas.matrix$best.period
            Weights=1
            return(list(lag=lag,Weights=Weights))
        }

      Observation.sum = sum(1/sqrt(M[,1]))
      Observation.weighting = (1/sqrt(M[,1]))

      Lag.sum = sum(M[,3]-M[,2])
      Lag.weighting = (M[,3]-M[,2])


      Weights = (Lag.weighting*Observation.weighting) / sum(Lag.weighting*Observation.weighting)

      # Determine lag from seasonality test
      if(seasonal.factor==FALSE){lag<- M[,1]}
      if(is.numeric(seasonal.factor)){
        lag<- seasonal.factor
        Weights=1}
      return(list(lag=lag,Weights=Weights))
    }}

  #Vectors generator for 1:lag
  generate.vectors=function(lag){
    Component.series = list()
    Component.index = list()

    for (i in 1:length(lag)){
      CS=rev(variable[seq(aa+1,1,-lag[i])])
      CS=CS[!is.na(CS)]
      Component.series[[paste('Series.',i,sep="")]] <- CS
      Component.index[[paste('Index.',i,sep="")]] <- (1:length(CS))
    }
    return(list(Component.index=Component.index,Component.series=Component.series))
  }



  if(dynamic == TRUE){
    for (j in 0:(h-1)){
      if(is.numeric(seasonal.factor)){
        M<-t(seasonal.factor)
        lag=seasonal.factor
        Weights=1
        seasonal.plot=F
      } else {

        seas.matrix = NNS.seas(variable,plot=F)
        if(!is.list(seas.matrix)){
          M<- t(1)} else {
            M<- seas.matrix$all.periods}
        ASW=ARMA.seas.weighting()
        lag=ASW$lag
        Weights=ASW$Weights
      }

      a=length(FV)
      aa=length(variable)

      # Generate vectors for 1:lag
      GV=generate.vectors(lag)
      Component.index=GV$Component.index
      Component.series=GV$Component.series
      # Regression on Component Series
      Regression.Estimates = numeric()
      Coefficients = numeric()
      Estimate.band= list()


      if(method=='nonlin' | method=='both'){

          Regression.Estimates=sapply(1:length(lag),function(i) NNS.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),return.values = TRUE,order = NULL,plot = FALSE,stn=0)$Point.est)

        if(negative.values==FALSE){
          Regression.Estimates=pmax(0,Regression.Estimates)
        }

        Nonlin.estimates=sum(Regression.Estimates*Weights)
      }#Linear == F

      if(method=='lin' | method=='both'){

          Regression.Estimates=sapply(1:length(lag),function(i) coef(lm(Component.series[[i]]~Component.index[[i]]))[1]+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1)))
        if(negative.values==FALSE){
          Regression.Estimates=pmax(0,Regression.Estimates)
        }

        Lin.estimates=sum(Regression.Estimates*Weights)
      }#Linear==T

      if(intervals==TRUE){
        if(method=='both'){
          Estimate.band[[j+1]]=c(Nonlin.estimates,Lin.estimates)}
        if(method=='nonlin'){
          Estimate.band[[j+1]]=Nonlin.estimates}
        if(method=='lin'){
          Estimate.band[[j+1]]=Lin.estimates}
      }


      if(method=='both'){
        Estimates[j+1]=mean(c(Lin.estimates,Nonlin.estimates))}
        else{
            Estimates[j+1] = sum(Regression.Estimates*Weights)
        }

      variable = c(variable,(Estimates[j+1]))
      FV=variable

    } #j loop
  }  #dynamic {}

  else {
    if(is.numeric(seasonal.factor)){
      M<-t(seasonal.factor)
      lag=seasonal.factor
      Weights=1
      seasonal.plot=F
    } else {

    seas.matrix = NNS.seas(variable,plot=F)
    if(!is.list(seas.matrix)){
      M<- t(1)} else {
        M<- seas.matrix$all.periods}
    ASW=ARMA.seas.weighting()
    lag=ASW$lag
    Weights=ASW$Weights
    }

    Estimate.band= list()

    for (j in 0:(h-1)){
      a=length(FV)
      aa=length(variable)

      # Generate vectors for 1:lag
      GV=generate.vectors(lag)
      Component.index=GV$Component.index
      Component.series=GV$Component.series

      # Regression on Component Series
      Regression.Estimates = numeric()
      Coefficients = numeric()


      if(method=='nonlin' | method=='both'){

          Regression.Estimates=sapply(1:length(lag),function(i) NNS.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),return.values = TRUE,order = NULL,plot = FALSE,stn=0)$Point.est)

        if(negative.values==FALSE){
          Regression.Estimates=pmax(0,Regression.Estimates)
        }

        NL.Regression.Estimates=Regression.Estimates
        Nonlin.estimates=sum(Regression.Estimates*Weights)
      }#Linear == F

      if(method=='lin' | method=='both'){

          Regression.Estimates=sapply(1:length(lag),function(i) coef(lm(Component.series[[i]]~Component.index[[i]]))[1]+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1)))

        if(negative.values==FALSE){
          Regression.Estimates=pmax(0,Regression.Estimates)
        }

        L.Regression.Estimates=Regression.Estimates
        Lin.estimates=sum(Regression.Estimates*Weights)

      }#Linear==T

      if(intervals==TRUE){
        if(method=='both'){
          Estimate.band[[j+1]]=c(NL.Regression.Estimates,L.Regression.Estimates)}
        if(method=='nonlin'){
          Estimate.band[[j+1]]=NL.Regression.Estimates}
        if(method=='lin'){
          Estimate.band[[j+1]]=L.Regression.Estimates}
      }

      if(method=='both'){
        Estimates[j+1]=mean(c(Lin.estimates,Nonlin.estimates))}
        else{
            Estimates[j+1] = sum(Regression.Estimates*Weights)
        }

      variable = c(variable,(Estimates[j+1]))
      FV=variable

    } # j loop non-dynamic
    }  # ELSE from (if dynamic)

  #### PLOTTING
  if(plot==T){
    if(seasonal.plot==T){
    if(ncol(M)>1){
      par(mfrow=c(2,1))
      plot(M[,Period],M[,Coefficient.of.Variance],
           xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",         ylim = c(0,2*M[1,Variable.Coefficient.of.Variance]))

      points(M[1,Period],M[1,Coefficient.of.Variance],pch=19,col='red')

      abline(h=M[1,Variable.Coefficient.of.Variance], col="red",lty=5)
      text(median(M[,Period]),M[1,Variable.Coefficient.of.Variance],pos=3,"Variable Coefficient of Variance",col='red')
    } else {
      par(mfrow=c(2,1))
      plot(1,1,pch=19,col='blue', xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",         ylim = c(0,2*abs(sd(FV)/mean(FV))))

      text(1,abs(sd(FV)/mean(FV)),pos=3,"NO SEASONALITY DETECTED",col='red')
    }
    }
    label=names(variable)
    plot(OV, type = 'l',lwd=2,main = "Forecast",col='steelblue',
          xlim=c(1,max((training.set+h),length(OV))),
          ylab=label, ylim=c(min(Estimates, OV),max( OV,Estimates)))

    if(intervals==T){
      for(i in 1:h){
        ys=unlist(Estimate.band[[i]])
        points(rep(training.set+i,length(ys)),ys,col=rgb(1, 0, 0, 0.0125),pch=15)
      }
      lines((training.set+1):(training.set+h),Estimates,type = 'l',lwd=2,lty=1,col='red')

      segments(training.set,FV[training.set],training.set+1,Estimates[1],lwd=2,lty=1,col='red')
      legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,1),col=c('steelblue','red'),lwd=2)

    } else{

      if(training.set[1]<length(OV)){
        lines((training.set+1):(training.set+h),Estimates,type = 'l',lwd=2,lty=3,col='red')

        segments(training.set,FV[training.set],training.set+1,Estimates[1],lwd=2,lty=3,col='red')
        legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,2),col=c('steelblue','red'),lwd=2)

      } else {
        lines((training.set+1):(training.set+h),Estimates,type = 'l',lwd=2,lty=1,col='red')

        segments(training.set,FV[training.set],training.set+1,Estimates[1],lwd=2,lty=1,col='red')
        legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,1),col=c('steelblue','red'),lwd=2)
      }


    }
    points(training.set,variable[training.set],col="green",pch=18)
    points(training.set+h,sum(Regression.Estimates*Weights),col="green",pch=18)

    par(mfrow=c(1,1))
  }
  return(Estimates)

}
