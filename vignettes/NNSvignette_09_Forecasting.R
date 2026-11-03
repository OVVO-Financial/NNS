## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, cache = FALSE)
library(NNS)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)


## ----linear,fig.width=5,fig.height=3,fig.align = "center", warning=FALSE------
nns_lin = NNS.ARMA(AirPassengers,
               h = 44,
               training.set = 100,
               method = "lin",
               plot = TRUE,
               seasonal.factor = 12,
               seasonal.plot = FALSE)

sqrt(mean((nns_lin - tail(AirPassengers, 44)) ^ 2))


## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"------------------
nns_nonlin = NNS.ARMA(AirPassengers,
               h = 44,
               training.set = 100,
               method = "nonlin",
               plot = TRUE,
               seasonal.factor = 12,
               seasonal.plot = FALSE)

sqrt(mean((nns_nonlin - tail(AirPassengers, 44)) ^ 2))


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


## ----optimizer, fig.width=6,fig.height=4,fig.align="center"-------------------
optimizer_periods <- c(12L, 18L, 24L)
nns.optimal <- NNS.ARMA.optim(training,
                              training.set = inner_n, h = 44,
                              seasonal.factor = optimizer_periods,
                              obj.fn = expression(sqrt(mean((predicted - actual)^2))),
                              objective = "min",
                              ncores = 1, print.trace = FALSE,
                              pred.int = .95, plot = TRUE)

nns.optimal


## ----optimizer-test-score-----------------------------------------------------
c(internal_validation_RMSE = nns.optimal$obj.fn,
  raw_reserved_test_RMSE = sqrt(mean((nns.optimal$raw.results - tail(AirPassengers, 44))^2)),
  corrected_reserved_test_RMSE = sqrt(mean((nns.optimal$results - tail(AirPassengers, 44))^2)))


## ----extension,fig.width=6,fig.height=4,fig.align="center"--------------------
extended <- NNS.ARMA.optim(AirPassengers,
                            seasonal.factor = seq(12, 36, 6),
                            obj.fn = expression(sqrt(mean((predicted - actual)^2))),
                            objective = "min",
                            ncores = 1, print.trace = FALSE,
                            pred.int = .95, h = 50, plot = TRUE)

head(data.frame(horizon = seq_along(extended$results),
                forecast = extended$results,
                lower = extended$lower.pred.int,
                upper = extended$upper.pred.int))

