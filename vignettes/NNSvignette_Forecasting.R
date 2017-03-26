## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(NNS)
require(knitr)
require(rgl)

## ----linear--------------------------------------------------------------
nns=NNS.ARMA(AirPassengers,h=44,training.set = 100,method='lin',plot = TRUE,seasonal.plot = FALSE)
sqrt(mean((nns-tail(AirPassengers,44))^2))

## ----nonlinear-----------------------------------------------------------
nns=NNS.ARMA(AirPassengers,h=44,training.set = 100,method='nonlin',plot=TRUE,seasonal.plot = FALSE)
sqrt(mean((nns-tail(AirPassengers,44))^2))

## ----seasonal test-------------------------------------------------------
seas=t(sapply(1:25,function(i) c(i,sqrt(mean((NNS.ARMA(AirPassengers,h=44,training.set = 100,method='lin',seasonal.factor=i,plot=FALSE)-tail(AirPassengers,44))^2)))))
colnames(seas)=c("Period","RMSE")
seas

## ----best fit------------------------------------------------------------
a=seas[which.min(seas[,2]),1]

## ----best nonlinear------------------------------------------------------
nns=NNS.ARMA(AirPassengers,h=44,training.set = 100,method='nonlin',seasonal.factor = a,plot = TRUE,seasonal.plot = FALSE)
sqrt(mean((nns-tail(AirPassengers,44))^2))

## ----best both-----------------------------------------------------------
nns=NNS.ARMA(AirPassengers,h=44,training.set = 100,method='both',seasonal.factor = a,plot=TRUE,seasonal.plot=FALSE)
sqrt(mean((nns-tail(AirPassengers,44))^2))

## ----extension,results='hide'--------------------------------------------
NNS.ARMA(AirPassengers,h=50,seasonal.factor = a,method = 'both',plot = TRUE,seasonal.plot = FALSE)

