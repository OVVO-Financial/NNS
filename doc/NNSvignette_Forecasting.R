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
require(meboot)

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
#  [1] 19.1362

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
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 19.1362358134362"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 19.1362358134362"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 26.9572377070982"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 26.9572377070982"
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 19.13624
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
#  [1] 0
#  
#  $errors
#   [1] -14.4760573 -21.8681345 -19.6868048 -32.7495054 -23.5509967 -17.4013923 -13.7046909  -6.2483756
#   [9]  -2.2442482   2.3982320  17.3238474  24.2539307   4.6698746  -2.6774464  -8.7261360 -35.8160320
#  [17]   6.1433263  -2.9982440   4.5157617  18.8213406   2.7327818   0.5108170  -0.5554123   3.2115438
#  [25] -23.9853255  -5.7414136 -27.2036192 -53.3243747 -21.6360845 -23.6855725 -23.5649419 -22.1711063
#  [33] -29.0886859 -26.4416647  13.1916867 -34.5942818 -47.9230569 -35.0864901 -63.1557545 -64.2537725
#  [41] -35.8920252 -50.9463999 -27.9879160 -22.7959250
#  
#  $results
#   [1] 338.4777 405.5864 448.7348 437.5736 381.7921 325.4127 288.1276 326.6116 337.5286 320.0636 382.3870
#  [12] 372.3583 373.5542 452.6515 498.8287 485.5634 422.6217 358.1323 318.7048 360.0762 371.9976 351.8174
#  [23] 420.6722 409.0720 409.6204 500.8377 550.0593 534.7215 464.7289 391.8851 349.8560 394.5948 405.5181
#  [34] 382.3369 458.1976 444.8364 444.8414 547.7296 599.8279 582.5432 505.5431 424.5574 380.1445 427.8977
#  
#  $lower.pred.int
#   [1] 286.9866 354.0953 397.2437 386.0826 330.3010 273.9217 236.6365 275.1206 286.0375 268.5725 330.8960
#  [12] 320.8672 322.0632 401.1604 447.3377 434.0723 371.1306 306.6413 267.2137 308.5851 320.5065 300.3264
#  [23] 369.1811 357.5809 358.1294 449.3467 498.5683 483.2305 413.2378 340.3940 298.3649 343.1037 354.0271
#  [34] 330.8458 406.7066 393.3453 393.3503 496.2386 548.3368 531.0521 454.0521 373.0664 328.6535 376.4067
#  
#  $upper.pred.int
#   [1] 368.1141 435.2228 478.3711 467.2100 411.4285 355.0491 317.7639 356.2480 367.1650 349.6999 412.0234
#  [12] 401.9947 403.1906 482.2878 528.4651 515.1997 452.2581 387.7687 348.3411 389.7126 401.6339 381.4538
#  [23] 450.3086 438.7084 439.2568 530.4741 579.6957 564.3579 494.3652 421.5215 379.4924 424.2312 435.1545
#  [34] 411.9733 487.8340 474.4727 474.4777 577.3660 629.4642 612.1795 535.1795 454.1938 409.7809 457.5341

## ----optimsubs, echo = FALSE--------------------------------------------------
nns.optimal = list()
nns.optimal$periods = 12
nns.optimal$weights = NULL
nns.optimal$method = "nonlin"
nns.optimal$shrink = FALSE
nns.optimal$results = c(354.2580, 421.2452, 462.4395, 453.0669, 395.8280, 338.4172, 301.1178, 338.6083, 347.7440, 330.7530, 393.0655, 383.2619, 390.9250, 468.8563, 511.8161, 501.4936, 436.7415, 370.9154, 331.3098, 371.0849, 380.7716, 361.0259, 430.2580, 418.6685, 427.7316, 516.8815, 561.5732, 550.3086, 478.0325, 403.7194, 361.7944, 403.9807, 413.6136, 390.9586, 467.3674, 453.9804, 464.4469, 564.6356, 611.0813, 598.8694, 519.0765, 436.3233, 392.0875, 436.6022)

nns.optimal$bias.shift = 0

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

