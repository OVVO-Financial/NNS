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
#  [1] 19.21812

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
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 19.2181153361782"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 19.2181153361782"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 19.9790337412655"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 19.9790337412655"
#  
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 19.21812
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
#  [1] 10.3416
#  
#  $errors
#   [1] -12.0495905 -19.5023885 -18.2981119 -30.4665605 -21.9967015 -16.3628298 -12.6732257  -4.2894621  -2.6001984   2.4174837  16.6574755  24.0964052  12.0029210   7.8864972
#  [15]  -0.7526824 -26.4198893  13.6743157   1.1898601   9.1072756  24.6494525   6.3543872   4.5198310   3.7511736   6.9241735 -13.4927319  10.9518474 -12.5758246 -38.6502806
#  [29]  -9.6293956 -16.2385122 -15.9817320 -12.1192381 -21.4941585 -18.2787520  22.4564209 -26.8238096 -33.5539336 -11.7710337 -42.4668107 -43.6993219 -18.9496482 -40.3338256
#  [43] -17.1561519  -8.7338598
#  
#  $results
#   [1] 364.5996 431.5868 472.7811 463.4085 406.1696 348.7588 311.4594 351.7984 358.0857 341.3710 402.3232 393.3122 400.8951 479.1979 522.1577 511.8352 447.0831 381.2570 341.6514
#  [20] 385.4619 390.2282 370.9257 437.1106 427.1677 436.0942 526.1620 574.7971 559.7077 488.0471 415.7932 372.1043 421.7103 422.7644 400.6785 473.1211 461.8718 472.1623 573.5495
#  [39] 625.3008 607.9430 528.9780 448.9954 402.3864 456.2667
#  
#  $lower.pred.int
#   [1] 311.9511 378.9384 420.1327 410.7600 353.5212 296.1104 258.8110 299.1500 305.4372 288.7226 349.6748 340.6638 348.2466 426.5495 469.5092 459.1867 394.4347 328.6086 289.0030
#  [20] 332.8135 337.5797 318.2773 384.4622 374.5193 383.4458 473.5135 522.1487 507.0593 435.3987 363.1447 319.4559 369.0618 370.1160 348.0301 420.4727 409.2233 419.5139 520.9011
#  [39] 572.6524 555.2945 476.3296 396.3469 349.7380 403.6183
#  
#  $upper.pred.int
#   [1] 398.9146 465.9018 507.0961 497.7235 440.4846 383.0738 345.7744 386.1134 392.4007 375.6860 436.6382 427.6272 435.2101 513.5129 556.4727 546.1502 481.3981 415.5720 375.9664
#  [20] 419.7769 424.5432 405.2407 471.4256 461.4827 470.4092 560.4770 609.1121 594.0227 522.3621 450.1082 406.4193 456.0253 457.0794 434.9936 507.4361 496.1868 506.4773 607.8645
#  [39] 659.6159 642.2580 563.2930 483.3104 436.7015 490.5818
#  

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

