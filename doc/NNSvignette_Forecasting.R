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
#  [1] 18.14353

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
#  [1] 18.14353
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
#   [1] -14.5179654 -21.7666667 -19.0986652 -32.3954545 -23.1026335 -16.5032468 -14.4584586  -5.0481602  -2.7280816
#  [10]   2.0902778  16.2233073  23.7751736   8.3364651   4.7594727  -0.7438799 -28.8999586  12.7260481   1.7278815
#  [19]   7.2305779  24.7806328   6.0325355   3.9167637   3.0094798   6.3994115 -18.1291213   6.7874238 -14.0071585
#  [28] -41.9492985 -10.9071507 -15.8409856 -18.4207758 -12.2973363 -21.2007896 -18.3892305  21.5920295 -27.2149430
#  [37] -38.1985891 -16.0122636 -42.9634219 -46.7733009 -19.6083870 -39.1855794 -19.3920235  -8.2085387
#  
#  $results
#   [1] 349.3212 416.7167 460.8384 449.2091 393.6162 338.1364 297.5474 339.9394 347.4883 330.3750 391.1133 382.3281
#  [13] 383.2205 462.6023 511.8337 496.5334 434.8450 371.9915 327.5564 375.3827 379.2428 359.3780 425.2856 415.7766
#  [25] 416.4799 507.4915 561.5928 542.7681 475.1500 406.2466 356.8846 411.0125 413.0095 390.1160 461.0507 450.7479
#  [37] 452.5314 554.7255 613.9660 591.4534 517.3189 440.9503 387.5731 446.9758
#  
#  $lower.pred.int
#   [1] 325.7098 393.1052 437.2269 425.5977 370.0047 314.5249 273.9359 316.3280 323.8768 306.7636 367.5018 358.7167
#  [13] 359.6091 438.9908 488.2222 472.9220 411.2335 348.3800 303.9450 351.7713 355.6314 335.7665 401.6742 392.1652
#  [25] 392.8684 483.8801 537.9814 519.1566 451.5385 382.6352 333.2732 387.4010 389.3981 366.5045 437.4393 427.1365
#  [37] 428.9200 531.1140 590.3546 567.8420 493.7075 417.3388 363.9617 423.3643
#  
#  $upper.pred.int
#   [1] 372.9326 440.3281 484.4498 472.8205 417.2276 361.7478 321.1588 363.5508 371.0997 353.9864 414.7247 405.9396
#  [13] 406.8320 486.2137 535.4451 520.1449 458.4564 395.6029 351.1678 398.9941 402.8543 382.9894 448.8971 439.3880
#  [25] 440.0913 531.1030 585.2043 566.3795 498.7614 429.8580 380.4961 434.6239 436.6210 413.7274 484.6622 474.3593
#  [37] 476.1428 578.3369 637.5775 615.0648 540.9304 464.5617 411.1845 470.5872
#  

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

