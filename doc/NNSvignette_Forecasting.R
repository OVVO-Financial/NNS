## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
library(data.table)
data.table::setDTthreads(2L)
options(mc.cores = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 2)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
library(data.table)
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
#  nns_nonlin = NNS.ARMA(AirPassengers,
#                 h = 44,
#                 training.set = 100,
#                 method = "nonlin",
#                 plot = FALSE,
#                 seasonal.factor = 12,
#                 seasonal.plot = FALSE)
#  
#  sqrt(mean((nns_nonlin - tail(AirPassengers, 44)) ^ 2))

## ----nonlinearres, eval = FALSE-----------------------------------------------
#  [1] 20.55102

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

## ----best optim, eval=FALSE---------------------------------------------------
#  nns.optimal = NNS.ARMA.optim(AirPassengers,
#                               training.set = 100,
#                               seasonal.factor = seq(12, 60, 6),
#                               obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                               objective = "min",
#                               pred.int = .95, plot = TRUE)
#  
#  nns.optimal

## ----optimres, eval=FALSE-----------------------------------------------------
#  [1] "CURRNET METHOD: lin"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT lin OBJECTIVE FUNCTION = 35.3996540135277"
#  [1] "BEST method = 'lin', seasonal.factor = c( 12 )"
#  [1] "BEST lin OBJECTIVE FUNCTION = 35.3996540135277"
#  [1] "CURRNET METHOD: nonlin"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'nonlin' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 20.5510211522245"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 20.5510211522245"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 19.4534618627141"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 19.4534618627141"
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 19.45346
#  
#  $method
#  [1] "both"
#  
#  $shrink
#  [1] FALSE
#  
#  $nns.regress
#  [1] FALSE
#  
#  $bias.shift
#  [1] 8.983377
#  
#  $errors
#   [1] -14.240300 -20.899691 -17.699956 -31.623457 -22.080467 -15.972663
#   [7] -12.633377  -4.289462  -2.371119   2.417484  16.657475  24.096405
#  [13]   8.926263   6.535789   1.826191 -27.314316  14.649727   2.812027
#  [19]  10.087412  26.318135   8.130503   5.994088   5.810046   8.968724
#  [25] -15.343831  11.107010  -9.030058 -37.870074  -6.868421 -13.948830
#  [31] -13.833241 -10.076019 -19.089073 -16.278753  25.441499 -23.904395
#  [37] -35.211740 -11.322375 -38.211436 -42.494907 -15.487474 -37.670592
#  [43] -14.477746  -6.587231
#  
#  $results
#   [1] 349.7431 410.0837 456.2834 444.3599 390.9029 340.0107 301.3500 340.6939
#   [9] 346.6123 329.4009 387.6409 381.0798 376.5337 442.1913 493.9687 479.2886
#  [17] 421.0377 366.8738 324.8321 367.6235 371.9435 351.7882 414.2779 408.6346
#  [25] 405.9523 477.0303 533.3005 516.6684 452.9214 394.5738 349.1661 394.8069
#  [33] 397.6765 374.3703 441.9106 436.5026 434.6776 511.4234 571.8929 553.3783
#  [41] 484.1285 421.4944 373.0609 421.2843
#  
#  $lower.pred.int
#   [1] 302.5739 362.9145 409.1142 397.1907 343.7337 292.8415 254.1808 293.5247
#   [9] 299.4430 282.2317 340.4716 333.9106 329.3645 395.0221 446.7995 432.1194
#  [17] 373.8685 319.7046 277.6629 320.4543 324.7743 304.6190 367.1087 361.4654
#  [25] 358.7831 429.8611 486.1313 469.4992 405.7521 347.4046 301.9968 347.6377
#  [33] 350.5072 327.2011 394.7414 389.3334 387.5084 464.2542 524.7237 506.2091
#  [41] 436.9593 374.3252 325.8917 374.1151
#  
#  $upper.pred.int
#   [1] 384.0671 444.4077 490.6074 478.6839 425.2269 374.3347 335.6740 375.0179
#   [9] 380.9363 363.7249 421.9648 415.4038 410.8577 476.5153 528.2927 513.6126
#  [17] 455.3617 401.1978 359.1561 401.9475 406.2675 386.1122 448.6019 442.9586
#  [25] 440.2763 511.3543 567.6245 550.9924 487.2453 428.8978 383.4900 429.1309
#  [33] 432.0004 408.6943 476.2346 470.8266 469.0016 545.7474 606.2169 587.7023
#  [41] 518.4525 455.8184 407.3849 455.6083

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

