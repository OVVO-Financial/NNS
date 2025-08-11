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
#  [1] 18.15208

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
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 18.1435264878535"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 18.1435264878535"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 20.8560044654062"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 20.8560044654062"
#  
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 18.15208
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
#  [1] -8.576982
#  
#  $errors
#   [1]  -5.6787879  -5.2833333  -4.1616162 -17.7909091 -10.3838384  -8.8636364  -7.4526316   3.9393939   7.4882812  12.3750000  29.1132812  34.3281250  20.2205492
#  [14]  27.6022786  20.8336687  -8.4665838  30.8449534  12.9914773  17.5563939  38.3826941  19.2903993  17.4644272  19.3331767  19.8155057  -3.4480488  35.5619032
#  [27]  13.5978472 -16.1723154  12.1689345  -0.7539891  -5.0831451   5.9867956  -3.9068174  -0.7986170  42.1995863 -10.1324609 -19.3155737  19.8071364  -8.0478172
#  [40] -14.4690520   9.3426681 -20.0538349  -2.4281117  14.8998761
#  
#  $results
#   [1] 340.7442 408.1397 452.2614 440.6321 385.0392 329.5594 288.9704 331.3624 338.9113 321.7980 382.5363 373.7511 374.6436 454.0253 503.2567 487.9564 426.2680
#  [18] 363.4145 318.9794 366.8057 370.7134 350.8874 416.7562 407.2385 407.9750 498.9849 553.0209 534.2507 466.5920 397.6690 348.3399 402.4098 404.5162 381.6244
#  [35] 452.6226 442.2906 444.1074 546.2302 605.3752 582.9540 508.7657 432.3692 378.9949 438.3229
#  
#  $lower.pred.int
#   [1] 293.9961 361.3916 405.5133 393.8840 338.2911 282.8113 242.2223 284.6143 292.1632 275.0499 335.7882 327.0030 327.8955 407.2772 456.5086 441.2083 379.5199
#  [18] 316.6664 272.2313 320.0576 323.9653 304.1393 370.0081 360.4904 361.2269 452.2368 506.2727 487.5026 419.8438 350.9209 301.5918 355.6617 357.7681 334.8763
#  [35] 405.8745 395.5424 397.3593 499.4820 558.6271 536.2058 462.0176 385.6211 332.2468 391.5748
#  
#  $upper.pred.int
#   [1] 387.4923 454.8878 499.0095 487.3802 431.7873 376.3075 335.7185 378.1105 385.6594 368.5461 429.2844 420.4993 421.3917 500.7734 550.0048 534.7046 473.0161
#  [18] 410.1626 365.7275 413.5538 417.4615 397.6356 463.5043 453.9866 454.7231 545.7330 599.7690 580.9988 513.3401 444.4171 395.0880 449.1579 451.2643 428.3725
#  [35] 499.3707 489.0387 490.8556 592.9783 652.1233 629.7021 555.5138 479.1173 425.7430 485.0710
#  

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

