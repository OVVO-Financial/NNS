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
#  [1] 20.19599

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
#  [1] 19.50822
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
#  [1] 11.34026
#  
#  $errors
#   [1] -12.0495905 -19.5023885 -18.2981119 -30.4665605 -21.9967015 -16.3628298 -12.6732257  -4.2894621  -2.6001984
#  [10]   2.4174837  16.6574755  24.0964052  12.0029210   7.8864972  -0.7526824 -26.4198893  13.6743157   1.1898601
#  [19]   9.1072756  24.6494525   7.8148305   5.9940877   5.8100458   8.9687243 -11.4805859  12.7282091 -12.4809879
#  [28] -36.8363281  -8.2269378 -16.1171482 -15.1770286 -12.3754742 -19.5291985 -16.2952067  25.2265398 -24.0729594
#  [37] -30.8466826  -9.3810198 -42.3392122 -41.2587312 -17.0627050 -40.1705358 -16.0734602  -9.0786139
#  
#  $results
#   [1] 354.2907 413.8379 458.0421 447.8737 393.3436 341.9774 303.6670 343.0508 348.7401 331.7577 389.9977 383.4367
#  [13] 380.9157 445.1133 493.9374 481.8833 422.2942 367.7482 326.1741 368.2520 374.0967 354.1451 416.6348 410.9915
#  [25] 410.2864 480.0025 531.6140 519.3302 454.3068 395.5812 350.6175 395.6273 399.2349 376.1436 443.3473 438.1093
#  [37] 438.3493 513.6632 569.8738 555.3985 485.0235 422.3001 374.1886 422.0488
#  
#  $lower.pred.int
#   [1] 301.7733 361.3205 405.5248 395.3563 340.8262 289.4601 251.1497 290.5334 296.2227 279.2404 337.4804 330.9193
#  [13] 328.3984 392.5960 441.4200 429.3659 369.7768 315.2308 273.6567 315.7346 321.5793 301.6277 364.1174 358.4741
#  [25] 357.7690 427.4851 479.0966 466.8128 401.7894 343.0638 298.1001 343.1099 346.7175 323.6262 390.8299 385.5920
#  [37] 385.8320 461.1459 517.3564 502.8811 432.5061 369.7827 321.6712 369.5314
#  
#  $upper.pred.int
#   [1] 390.2389 449.7861 493.9904 483.8219 429.2918 377.9257 339.6153 378.9990 384.6883 367.7060 425.9460 419.3849
#  [13] 416.8640 481.0616 529.8856 517.8315 458.2425 403.6964 362.1223 404.2002 410.0450 390.0933 452.5830 446.9397
#  [25] 446.2347 515.9507 567.5622 555.2784 490.2550 431.5294 386.5657 431.5755 435.1831 412.0918 479.2956 474.0576
#  [37] 474.2976 549.6115 605.8220 591.3467 520.9717 458.2483 410.1368 457.9970

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

