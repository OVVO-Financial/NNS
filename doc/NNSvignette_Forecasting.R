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
#  [1] 19.02029

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
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 20.1959877511828"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 20.1959877511828"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 19.5082249052739"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 19.5082249052739"
#  
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 19.02029
#  
#  $method
#  [1] "nonlin"
#  
#  $shrink
#  [1] FALSE
#  
#  $nns.regress
#  [1] FALSE
#  
#  $bias.shift
#  [1] 12.55657
#  
#  $errors
#   [1] -12.0495905 -19.5023885 -18.2981119 -30.4665605 -21.9967015 -16.3628298 -12.6732257  -6.5326720  -2.6001984   2.4174837  16.6574755  24.0964052  12.0029210   7.8864972
#  [15]  -0.7526824 -26.4198893  13.6743157   1.1898601   9.1072756  21.4715719   6.7813958   5.0184404   4.0969001   7.2867135 -13.0622228  11.2980457 -13.7030809 -38.2538043
#  [29]  -9.2828407 -16.8050053 -15.9037050 -16.7472123 -20.9196379 -17.6078957  22.9215801 -26.3360285 -32.9747032 -11.3052396 -43.9834829 -43.1658811 -18.4833744 -41.0960163
#  [43] -17.0511702 -14.4745600
#  
#  $results
#   [1] 366.8145 433.8018 474.9961 465.6235 408.3846 350.9738 313.6744 349.5269 360.3006 343.5860 404.5382 395.5272 403.1100 481.4129 524.3726 514.0501 449.2981 383.4720 343.8664
#  [20] 381.3211 393.2971 374.1379 440.0170 430.1078 439.1702 529.0693 574.7576 562.7156 490.9552 416.8751 374.4754 414.6693 426.1284 404.2352 476.2664 465.0623 475.5357 576.6961
#  [39] 624.4825 611.2248 532.1255 449.6860 404.8114 447.0003
#  
#  $lower.pred.int
#   [1] 311.2473 378.2346 419.4288 410.0562 352.8174 295.4066 258.1072 293.9597 304.7334 288.0188 348.9710 339.9599 347.5428 425.8457 468.8054 458.4829 393.7308 327.9048 288.2991
#  [20] 325.7539 337.7299 318.5707 384.4498 374.5406 383.6030 473.5021 519.1903 507.1484 435.3880 361.3079 318.9081 359.1021 370.5612 348.6680 420.6992 409.4951 419.9685 521.1289
#  [39] 568.9153 555.6576 476.5583 394.1188 349.2442 391.4331
#  
#  $upper.pred.int
#   [1] 402.1839 469.1712 510.3655 500.9929 443.7540 386.3432 349.0438 384.8963 395.6700 378.9554 439.9076 430.8966 438.4794 516.7823 559.7420 549.4195 484.6675 418.8414 379.2358
#  [20] 416.6905 428.6665 409.5073 475.3864 465.4772 474.5396 564.4387 610.1270 598.0850 526.3246 452.2445 409.8448 450.0387 461.4978 439.6046 511.6358 500.4317 510.9051 612.0655
#  [39] 659.8519 646.5942 567.4949 485.0554 440.1808 482.3697

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

