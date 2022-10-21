## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(mc.cores=2)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

## ----linear,fig.width=5,fig.height=3,fig.align = "center", warning=FALSE------
nns = NNS.ARMA(AirPassengers, 
               h = 44, 
               training.set = 100, 
               method = "lin", 
               plot = TRUE, 
               seasonal.factor = 12, 
               seasonal.plot = FALSE)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"------------------
nns = NNS.ARMA(AirPassengers, 
               h = 44, 
               training.set = 100, 
               method = "nonlin", 
               plot = FALSE, 
               seasonal.factor = 12, 
               seasonal.plot = FALSE)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----seasonal test, eval=TRUE-------------------------------------------------
seas = t(sapply(1 : 25, function(i) c(i, sqrt( mean( (NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "lin", seasonal.factor = i, plot=FALSE) - tail(AirPassengers, 44)) ^ 2) ) ) ) )

colnames(seas) = c("Period", "RMSE")
seas

## ----best fit, eval=TRUE------------------------------------------------------
a = seas[which.min(seas[ , 2]), 1]

## ----best nonlinear,fig.width=5,fig.height=3,fig.align = "center", eval=TRUE----
nns = NNS.ARMA(AirPassengers, 
               h = 44, 
               training.set = 100, 
               method = "nonlin", 
               seasonal.factor = a, 
               plot = TRUE, seasonal.plot = FALSE)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----modulo, eval=TRUE--------------------------------------------------------
NNS.seas(AirPassengers, modulo = 12, plot = FALSE)

## ----best optim---------------------------------------------------------------
nns.optimal = NNS.ARMA.optim(AirPassengers, 
                             training.set = 100, 
                             seasonal.factor = seq(12, 24, 6),
                             obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
                             objective = "min",
                             conf.intervals = .95)

nns.optimal

## ----best optim2, eval=TRUE---------------------------------------------------
sqrt(mean((nns.optimal$results - tail(AirPassengers, 44)) ^ 2))

## ----best optim3, eval=TRUE---------------------------------------------------
sqrt(mean((nns+nns.optimal$bias.shift - tail(AirPassengers, 44)) ^ 2))

## ----neg, eval=TRUE-----------------------------------------------------------
nns <- pmax(0, nns+nns.optimal$bias.shift)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=TRUE----
NNS.ARMA(AirPassengers, 
         h = 50,
         conf.intervals = .95,
         seasonal.factor = nns.optimal$periods, 
         method  = nns.optimal$method, 
         weights = nns.optimal$weights, 
         plot = TRUE, seasonal.plot = FALSE) + nns.optimal$bias.shift

