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
#  [1] 18.77889

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
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 18.7788946994645"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 18.7788946994645"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 26.9404229568573"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 26.9404229568573"
#  
#  
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 18.77889
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
#   [1] -14.47605733 -21.86813452 -19.68680481 -32.74950538 -23.55099666 -17.40139233 -13.70469091  -6.24837561  -2.48949314   2.13237206  16.92258707  23.90974592   4.53920624
#  [14]  -2.80578226  -8.66505386 -35.95512814   6.15047733  -3.05392579   4.47214630  18.85551669   2.32637841   0.08983381  -1.04028830   2.89582233 -24.22859693  -6.13584252
#  [27] -27.38320978 -53.70280267 -21.86253259 -23.90679863 -23.77369245 -22.36041401 -29.22056474 -26.31507238  12.93035387 -34.59586321 -47.79031322 -34.85559482 -62.85716413
#  [40] -64.07451099 -35.67072562 -50.91910417 -27.91240080 -22.72623962
#  
#  $results
#   [1] 338.5400 405.6487 448.7971 437.6360 381.8544 325.4751 288.1899 326.6740 337.1004 319.5942 381.6468 371.7323 373.1716 452.1286 498.6929 485.0008 422.3997 357.9051
#  [19] 318.5071 360.0295 371.1142 350.8923 419.5388 408.2967 408.8936 499.6842 549.5283 533.5434 464.0692 391.2918 349.2869 394.1562 404.3549 381.6612 456.5031 443.9546
#  [37] 443.9982 546.2672 599.0085 581.0612 504.6510 423.7409 379.3403 427.2016
#  
#  $lower.pred.int
#   [1] 287.2344 354.3431 397.4915 386.3304 330.5488 274.1695 236.8843 275.3684 285.7949 268.2886 330.3413 320.4267 321.8661 400.8231 447.3873 433.6952 371.0941 306.5995
#  [19] 267.2016 308.7239 319.8086 299.5867 368.2332 356.9912 357.5880 448.3786 498.2227 482.2379 412.7636 339.9862 297.9813 342.8507 353.0493 330.3556 405.1975 392.6490
#  [37] 392.6926 494.9616 547.7030 529.7556 453.3454 372.4353 328.0347 375.8961
#  
#  $upper.pred.int
#   [1] 368.1156 435.2243 478.3727 467.2115 411.4300 355.0506 317.7655 356.2495 366.6760 349.1697 411.2224 401.3078 402.7472 481.7042 528.2685 514.5763 451.9753 387.4807
#  [19] 348.0827 389.6050 400.6897 380.4679 449.1144 437.8723 438.4692 529.2597 579.1039 563.1190 493.6448 420.8674 378.8625 423.7318 433.9305 411.2367 486.0786 473.5301
#  [37] 473.5738 575.8427 628.5841 610.6368 534.2265 453.3165 408.9159 456.7772

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

