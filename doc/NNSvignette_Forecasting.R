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
#  [1] 18.17333

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
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 18.17333
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
#   [1] -14.5179654 -21.7666667 -19.0986652 -32.3954545 -23.1026335 -16.5032468 -14.4584586  -5.0481602  -3.2533610
#  [10]   2.0902778  16.2233073  23.7751736   8.3364651   3.7503262  -0.7438799 -29.8067346  12.7260481   1.7278815
#  [19]   7.2305779  24.7806328   5.3729087   4.3178513   4.4136936   6.3994115 -16.9198724   5.5201743 -14.3855904
#  [28] -43.0879949 -11.1688248 -15.8409856 -18.0183676 -12.6409715 -22.1182760 -17.9006597  23.3339776 -27.7393952
#  [37] -36.6984978 -17.7822175 -44.0355166 -48.2781513 -20.4198570 -39.2330543 -18.8928301  -8.6348233
#  
#  $results
#   [1] 349.3212 416.7167 460.8384 449.2091 393.6162 338.1364 297.5474 339.9394 346.4377 330.3750 391.1133 382.3281 383.2205
#  [14] 460.5840 511.8337 494.7199 434.8450 371.9915 327.5564 375.3827 377.9236 360.1801 428.0941 415.7766 418.8984 504.9570
#  [27] 560.8360 540.4907 474.6266 406.2466 357.6895 410.3252 411.1746 391.0931 464.5346 449.6990 455.5316 551.1856 611.8218
#  [40] 588.4437 515.6960 440.8553 388.5715 446.1232
#  
#  $lower.pred.int
#   [1] 325.5791 392.9746 437.0963 425.4670 369.8741 314.3943 273.8053 316.1973 322.6956 306.6329 367.3712 358.5860 359.4785
#  [14] 436.8419 488.0916 470.9778 411.1029 348.2494 303.8143 351.6406 354.1815 336.4381 404.3520 392.0345 395.1563 481.2149
#  [27] 537.0939 516.7486 450.8846 382.5045 333.9474 386.5831 387.4325 367.3510 440.7925 425.9569 431.7895 527.4435 588.0797
#  [40] 564.7016 491.9539 417.1132 364.8294 422.3811
#  
#  $upper.pred.int
#   [1] 373.0633 440.4588 484.5805 472.9512 417.3582 361.8784 321.2895 363.6815 370.1798 354.1171 414.8554 406.0702 406.9626
#  [14] 484.3261 535.5758 518.4619 458.5870 395.7336 351.2985 399.1248 401.6657 383.9222 451.8361 439.5187 442.6404 528.6991
#  [27] 584.5780 564.2328 498.3687 429.9887 381.4315 434.0673 434.9166 414.8352 488.2767 473.4411 479.2737 574.9276 635.5639
#  [40] 612.1858 539.4381 464.5974 412.3136 469.8653
#  

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

