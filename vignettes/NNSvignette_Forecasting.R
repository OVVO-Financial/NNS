## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(mc.cores=2)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

## ----linear,fig.width=5,fig.height=3,fig.align = "center", warning=FALSE------
nns = NNS.ARMA(AirPassengers, 
               h = 44, 
               training.set = 100, 
               method = "lin", 
               plot = TRUE, 
               seasonal.factor = 12, 
               seasonal.plot = FALSE)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center", eval = FALSE----
#  nns = NNS.ARMA(AirPassengers,
#                 h = 44,
#                 training.set = 100,
#                 method = "nonlin",
#                 plot = FALSE,
#                 seasonal.factor = 12,
#                 seasonal.plot = FALSE)
#  
#  sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----nonlinearres, eval = FALSE-----------------------------------------------
#  [1] 19.49762

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
#                               seasonal.factor = seq(12, 24, 6),
#                               obj.fn = expression( sqrt(mean((predicted - actual)^2)) ),
#                               objective = "min",
#                               conf.intervals = .95)
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
#  [1] "CURRENT nonlin OBJECTIVE FUNCTION = 19.4976178189546"
#  [1] "BEST method = 'nonlin' PATH MEMBER = c( 12 )"
#  [1] "BEST nonlin OBJECTIVE FUNCTION = 19.4976178189546"
#  [1] "CURRNET METHOD: both"
#  [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
#  [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 12 ) ...)"
#  [1] "CURRENT both OBJECTIVE FUNCTION = 26.6112299452096"
#  [1] "BEST method = 'both' PATH MEMBER = c( 12 )"
#  [1] "BEST both OBJECTIVE FUNCTION = 26.6112299452096"
#  $periods
#  [1] 12
#  
#  $weights
#  NULL
#  
#  $obj.fn
#  [1] 19.49762
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
#   [1] -12.0495905 -19.5023885 -18.2981119 -30.4665605 -21.9967015 -16.3628298
#   [7] -12.6732257  -5.7137170  -2.6001984   2.2792659  17.1994048  24.2420635
#  [13]   6.6919485  -1.2269250  -8.4029057 -34.4569779   6.9539623  -2.5920976
#  [19]   4.8338436  18.5863427   1.8098569  -0.3087157  -1.1892791   2.5325891
#  [25] -22.4687006  -4.9819699 -27.7262972 -52.7041072 -21.5667488 -23.9122298
#  [31] -23.6982624 -23.0856682 -29.9142644 -27.1628466  12.6507957 -35.1714729
#  [37] -46.1877025 -34.0820674 -63.4664903 -63.3893474 -35.6270575 -51.0256013
#  [43] -27.9853043 -23.5848310
#  
#  $results
#   [1] 354.2580 421.2452 462.4395 453.0669 395.8280 338.4172 301.1178 338.6083
#   [9] 347.7440 330.7530 393.0655 383.2619 390.9250 468.8563 511.8161 501.4936
#  [17] 436.7415 370.9154 331.3098 371.0849 380.7716 361.0259 430.2580 418.6685
#  [25] 427.7316 516.8815 561.5732 550.3086 478.0325 403.7194 361.7944 403.9807
#  [33] 413.6136 390.9586 467.3674 453.9804 464.4469 564.6356 611.0813 598.8694
#  [41] 519.0765 436.3233 392.0875 436.6022
#  
#  $lower.conf.int
#   [1] 291.6700 358.6573 399.8515 390.4789 333.2401 275.8292 238.5299 276.0203
#   [9] 285.1561 268.1650 330.4775 320.6740 328.3370 406.2684 449.2281 438.9056
#  [17] 374.1535 308.3275 268.7218 308.4970 318.1836 298.4380 367.6701 356.0805
#  [25] 365.1437 454.2935 498.9853 487.7206 415.4446 341.1314 299.2065 341.3928
#  [33] 351.0256 328.3706 404.7794 391.3924 401.8589 502.0477 548.4933 536.2814
#  [41] 456.4885 373.7353 329.4996 374.0142
#  
#  $upper.conf.int
#   [1] 416.8459 483.8332 525.0274 515.6548 458.4160 401.0052 363.7058 401.1962
#   [9] 410.3320 393.3409 455.6534 445.8499 453.5129 531.4443 574.4040 564.0815
#  [17] 499.3294 433.5034 393.8977 433.6729 443.3595 423.6139 492.8460 481.2564
#  [25] 490.3196 579.4694 624.1612 612.8965 540.6205 466.3074 424.3824 466.5687
#  [33] 476.2015 453.5465 529.9553 516.5683 527.0349 627.2236 673.6692 661.4573
#  [41] 581.6644 498.9112 454.6755 499.1901

## ----optimsubs, echo = FALSE--------------------------------------------------
nns.optimal = list()
nns.optimal$periods = 12
nns.optimal$weights = NULL
nns.optimal$method = "nonlin"
nns.optimal$shrink = FALSE
nns.optimal$results = c(354.2580, 421.2452, 462.4395, 453.0669, 395.8280, 338.4172, 301.1178, 338.6083, 347.7440, 330.7530, 393.0655, 383.2619, 390.9250, 468.8563, 511.8161, 501.4936, 436.7415, 370.9154, 331.3098, 371.0849, 380.7716, 361.0259, 430.2580, 418.6685, 427.7316, 516.8815, 561.5732, 550.3086, 478.0325, 403.7194, 361.7944, 403.9807, 413.6136, 390.9586, 467.3674, 453.9804, 464.4469, 564.6356, 611.0813, 598.8694, 519.0765, 436.3233, 392.0875, 436.6022)

nns.optimal$bias.shift = 0

## ----best optim2, eval=TRUE---------------------------------------------------
sqrt(mean((nns.optimal$results - tail(AirPassengers, 44)) ^ 2))

## ----best optim3, eval=TRUE---------------------------------------------------
sqrt(mean((nns+nns.optimal$bias.shift - tail(AirPassengers, 44)) ^ 2))

## ----neg, eval=TRUE-----------------------------------------------------------
nns <- pmax(0, nns + nns.optimal$bias.shift)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center", eval=TRUE----
NNS.ARMA(AirPassengers, 
         h = 50,
         conf.intervals = .95,
         seasonal.factor = nns.optimal$periods, 
         method  = nns.optimal$method, 
         weights = nns.optimal$weights, 
         shrink = nns.optimal$shrink,
         plot = TRUE, seasonal.plot = FALSE) + nns.optimal$bias.shift

