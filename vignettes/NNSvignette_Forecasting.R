## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
library(data.table)
data.table::setDTthreads(1L)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)

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
# nns_nonlin = NNS.ARMA(AirPassengers,
#                h = 44,
#                training.set = 100,
#                method = "nonlin",
#                plot = FALSE,
#                seasonal.factor = 12,
#                seasonal.plot = FALSE)
# 
# sqrt(mean((nns_nonlin - tail(AirPassengers, 44)) ^ 2))

## ----nonlinearres, eval = FALSE-----------------------------------------------
# [1] 18.1809

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
# nns.optimal = NNS.ARMA.optim(AirPassengers,
#                              training.set = 100,
#                              seasonal.factor = seq(12, 60, 6),
#                              obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                              objective = "min",
#                              pred.int = .95, plot = TRUE)
# 
# nns.optimal

## ----optimres, eval=FALSE-----------------------------------------------------
# [1] "CURRNET METHOD: lin"
# [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
# [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 12 ) ...)"
# [1] "CURRENT lin OBJECTIVE FUNCTION = 35.3996540135277"
# [1] "BEST method = 'lin', seasonal.factor = c( 12 )"
# [1] "BEST lin OBJECTIVE FUNCTION = 35.3996540135277"
# [1] "CURRNET METHOD: nonlin"
# [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
# [1] "NNS.ARMA(... method =  'nonlin' , seasonal.factor =  c( 12 ) ...)"
# [1] "CURRENT nonlin OBJECTIVE FUNCTION = 18.1809033101955"
# [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
# [1] "BEST nonlin OBJECTIVE FUNCTION = 18.1809033101955"
# [1] "CURRNET METHOD: both"
# [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
# [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
# [1] "CURRENT both OBJECTIVE FUNCTION = 22.7363330823967"
# [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
# [1] "BEST both OBJECTIVE FUNCTION = 22.7363330823967"
# >
# > nns.optimal
# $periods
# [1] 12
# 
# $weights
# NULL
# 
# $obj.fn
# [1] 18.1809
# 
# $method
# [1] "nonlin"
# 
# $shrink
# [1] FALSE
# 
# $nns.regress
# [1] FALSE
# 
# $bias.shift
# [1] 0
# 
# $errors
#  [1]  -6.0626221 -10.8434613 -10.7646998 -22.7134790 -15.3519569 -12.9673866  -9.1626428   3.9393939   7.4882812  12.3750000  29.1132812  34.3281250  19.7002739
# [14]  20.0656989  11.8833952 -15.1389735  24.1108241   7.4289721  15.2385271  38.3826941  19.2903993  17.4644272  19.3331767  19.8155057  -4.0856291  26.3260739
# [27]   2.6153110 -24.3491085   3.9057436  -8.8271346  -7.9236143   5.9867956  -3.9068174  -0.7986170  42.1995863 -10.1324609 -20.0852820   8.6573328 -21.3067790
# [40] -24.3403514  -0.6332912 -29.8418247  -5.8572216  14.8998761
# 
# $results
#  [1] 348.9374 411.1565 454.2353 444.2865 388.6480 334.0326 295.8374 339.9394 347.4883 330.3750 391.1133 382.3281 382.7003 455.0657 502.8834 489.8610 428.1108
# [18] 366.4290 325.2385 375.3827 379.2904 359.4644 425.3332 415.8155 415.9144 498.3261 550.6153 534.6509 466.9057 398.1729 354.0764 410.9868 413.0932 390.2014
# [35] 461.1996 450.8675 451.9147 543.6573 600.6932 581.6596 507.3667 431.1582 384.1428 446.8999
# 
# $lower.pred.int
#  [1] 310.8588 373.0779 416.1567 406.2079 350.5694 295.9540 257.7588 301.8608 309.4097 292.2964 353.0347 344.2495 344.6217 416.9871 464.8048 451.7824 390.0322
# [18] 328.3504 287.1599 337.3041 341.2118 321.3858 387.2546 377.7369 377.8358 460.2475 512.5367 496.5723 428.8271 360.0943 315.9978 372.9082 375.0146 352.1228
# [35] 423.1210 412.7889 413.8361 505.5787 562.6146 543.5810 469.2881 393.0796 346.0642 408.8213
# 
# $upper.pred.int
#  [1] 387.0160 449.2351 492.3139 482.3651 426.7266 372.1112 333.9160 378.0180 385.5669 368.4536 429.1919 420.4067 420.7789 493.1443 540.9620 527.9396 466.1894
# [18] 404.5076 363.3171 413.4613 417.3690 397.5430 463.4118 453.8941 453.9930 536.4047 588.6939 572.7295 504.9843 436.2515 392.1550 449.0654 451.1718 428.2800
# [35] 499.2782 488.9461 489.9933 581.7359 638.7718 619.7382 545.4453 469.2368 422.2214 484.9785
# 

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
# NNS.ARMA.optim(AirPassengers,
#                 seasonal.factor = seq(12, 60, 6),
#                 obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                 objective = "min",
#                 pred.int = .95, h = 50, plot = TRUE)

