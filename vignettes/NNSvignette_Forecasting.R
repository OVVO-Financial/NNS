## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(NNS)
require(knitr)
require(rgl)

## ----linear,fig.width=5,fig.height=3,fig.align = "center"----------------
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "lin", plot = TRUE, seasonal.factor = 12, seasonal.plot = FALSE)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"-------------
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "nonlin", plot = TRUE, seasonal.factor = 12, seasonal.plot = FALSE)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----seasonal test,fig.width=5,fig.height=3,fig.align = "center"---------
seas = t(sapply(1 : 25, function(i) c(i, sqrt(mean((NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "lin", seasonal.factor = i, plot=FALSE) - tail(AirPassengers, 44)) ^ 2)))))
colnames(seas) = c("Period", "RMSE")
seas

## ----best fit------------------------------------------------------------
a = seas[which.min(seas[ , 2]), 1]

## ----best nonlinear,fig.width=5,fig.height=3,fig.align = "center"--------
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "nonlin",seasonal.factor = a, plot = TRUE, seasonal.plot = FALSE)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----best both,fig.width=5,fig.height=3,fig.align = "center"-------------
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "both", seasonal.factor = a, plot = TRUE, seasonal.plot = FALSE)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----best optim,fig.width=5,fig.height=3,fig.align = "center"------------
nns.optimal = NNS.ARMA.optim(AirPassengers, training.set = 100, seasonal.factor = seq(12, 48, 12), method = "comb")

nns.optimal

## ----best optim2,fig.width=5,fig.height=3,fig.align = "center"-----------
nns = NNS.ARMA(AirPassengers, training.set = 100, h = 44, seasonal.factor = nns.optimal$periods, method = nns.optimal$method, plot = TRUE, seasonal.plot = FALSE)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center"----
NNS.ARMA(AirPassengers, h = 50, seasonal.factor = nns.optimal$periods, method  = nns.optimal$method, plot = TRUE, seasonal.plot = FALSE)

