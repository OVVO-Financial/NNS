## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)


## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
require(knitr)
require(rgl)


## ----linear,fig.width=5,fig.height=3,fig.align = "center", warning=FALSE------
nns_lin = NNS.ARMA(AirPassengers,
               h = 44,
               training.set = 100,
               method = "lin",
               plot = TRUE,
               seasonal.factor = 12,
               seasonal.plot = FALSE)

sqrt(mean((nns_lin - tail(AirPassengers, 44)) ^ 2))


## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center", eval = FALSE----
## nns_nonlin = NNS.ARMA(AirPassengers,
##                h = 44,
##                training.set = 100,
##                method = "nonlin",
##                plot = FALSE,
##                seasonal.factor = 12,
##                seasonal.plot = FALSE)
##
## sqrt(mean((nns_nonlin - tail(AirPassengers, 44)) ^ 2))


## ----seasonal-test, eval=TRUE-------------------------------------------------
training <- head(AirPassengers, 100)
inner_n <- 80L
candidates <- 1:16
seas <- data.frame(
  Period = candidates,
  validation_RMSE = vapply(candidates, function(period) {
    predicted <- NNS.ARMA(head(training, inner_n), h = length(training) - inner_n,
                          method = "lin", seasonal.factor = period,
                          plot = FALSE, seasonal.plot = FALSE)
    sqrt(mean((predicted - tail(training, length(training) - inner_n))^2))
  }, numeric(1))
)
a <- seas$Period[which.min(seas$validation_RMSE)]
seas


## ----selected-period-test, fig.width=5,fig.height=3,fig.align="center"--------
selected <- NNS.ARMA(training, h = 44, method = "lin", seasonal.factor = a,
                     plot = TRUE, seasonal.plot = FALSE)
sqrt(mean((selected - tail(AirPassengers, 44))^2))


## ----modulo, eval=TRUE--------------------------------------------------------
NNS.seas(training, modulo = 12, plot = FALSE)


## ----best optim, eval=FALSE---------------------------------------------------
## nns.optimal = NNS.ARMA.optim(training,
##                              training.set = 80, h = 44,
##                              seasonal.factor = seq(12, 60, 6),
##                              obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
##                              objective = "min",
##                              pred.int = .95, plot = TRUE)
##
## nns.optimal


## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
## NNS.ARMA.optim(AirPassengers,
##                 seasonal.factor = seq(12, 60, 6),
##                 obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
##                 objective = "min",
##                 pred.int = .95, h = 50, plot = TRUE)

