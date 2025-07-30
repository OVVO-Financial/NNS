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
#  [1] 18.26941

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
#  [1] 18.26941
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
#  [1] 11.40284
#  
#  $errors
#   [1] -14.517965 -21.766667 -18.691342 -32.395455 -22.867532 -16.503247 -13.189719  -5.048160  -2.728082   2.090278
#  [11]  16.223307  23.775174   8.336465   4.759473  -0.191766 -28.899959  13.044720   1.727881   8.950315  24.780633
#  [21]   6.032536   3.916764   3.009480   6.399411 -18.129121   6.787424 -13.327920 -41.949298 -10.515104 -15.840986
#  [31] -16.307810 -12.297336 -21.200790 -18.389231  21.592029 -27.214943 -38.198589 -16.012264 -42.144402 -46.773301
#  [41] -19.135661 -39.185579 -16.844322  -8.208539
#  
#  $results
#   [1] 360.7241 428.1195 473.0559 460.6119 405.4892 349.5392 311.4877 351.3422 358.8911 341.7778 402.5161 393.7310
#  [13] 394.6234 474.0051 524.3407 507.9363 446.8851 383.3943 342.3987 386.7855 390.6457 370.7808 436.6885 427.1794
#  [25] 427.8827 518.8944 574.3541 554.1709 487.3369 417.6494 372.5134 422.4153 424.4124 401.5188 472.4536 462.1507
#  [37] 463.9342 566.1283 627.0069 602.8562 529.6672 452.3531 404.0713 458.3786
#  
#  $lower.pred.int
#   [1] 325.7098 393.1052 438.0416 425.5977 370.4749 314.5249 276.4734 316.3280 323.8768 306.7636 367.5018 358.7167
#  [13] 359.6091 438.9908 489.3265 472.9220 411.8709 348.3800 307.3844 351.7713 355.6314 335.7665 401.6742 392.1652
#  [25] 392.8684 483.8801 539.3399 519.1566 452.3226 382.6352 337.4991 387.4010 389.3981 366.5045 437.4393 427.1365
#  [37] 428.9200 531.1140 591.9926 567.8420 494.6530 417.3388 369.0571 423.3643
#  
#  $upper.pred.int
#   [1] 395.7383 463.1338 508.0701 495.6262 440.5035 384.5535 346.5020 386.3565 393.9054 376.7921 437.5304 428.7452
#  [13] 429.6377 509.0194 559.3550 542.9505 481.8994 418.4086 377.4130 421.7998 425.6600 405.7951 471.7027 462.1937
#  [25] 462.8970 553.9086 609.3684 589.1852 522.3512 452.6637 407.5277 457.4296 459.4266 436.5331 507.4678 497.1650
#  [37] 498.9485 601.1426 662.0212 637.8705 564.6815 487.3674 439.0856 493.3929
#  

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=FALSE----
#  NNS.ARMA.optim(AirPassengers,
#                  seasonal.factor = seq(12, 60, 6),
#                  obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                  objective = "min",
#                  pred.int = .95, h = 50, plot = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

