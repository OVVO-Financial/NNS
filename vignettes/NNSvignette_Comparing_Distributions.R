## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
library(data.table)
data.table::setDTthreads(2L)
options(mc.cores = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 2)

## ----setup2,message=FALSE,warning = FALSE-------------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

## ----cars, fig.width=7, fig.align='center'------------------------------------
mpg_auto_trans = mtcars[mtcars$am==1, "mpg"]
mpg_man_trans = mtcars[mtcars$am==0, "mpg"]

NNS.ANOVA(control = mpg_man_trans, treatment = mpg_auto_trans, robust = TRUE)

## ----cars2, warning=FALSE-----------------------------------------------------
wilcox.test(mpg ~ am, data=mtcars) 

## ----equalmeans, echo=TRUE, fig.width=7, fig.align='center'-------------------
set.seed(123)
x = rnorm(1000, mean = 0, sd = 1)
y = rnorm(1000, mean = 0, sd = 2)

NNS.ANOVA(control = x, treatment = y,
          means.only = TRUE, robust = TRUE, plot = TRUE)

t.test(x,y)

## ----unequalmeans, echo=TRUE, fig.width=7, fig.align='center'-----------------
set.seed(123)
x = rnorm(1000, mean = 0, sd = 1)
y = rnorm(1000, mean = 1, sd = 1)

NNS.ANOVA(control = x, treatment = y,
          means.only = TRUE, robust = TRUE, plot = TRUE)

t.test(x,y)

## ----stochdom, fig.width=7, fig.align='center'--------------------------------
NNS.FSD(x, y)

## ----stochdomset, eval=FALSE--------------------------------------------------
#  set.seed(123)
#  x1 = rnorm(1000)
#  x2 = x1 + 1
#  x3 = rnorm(1000)
#  x4 = x3 + 1
#  x5 = rnorm(1000)
#  x6 = x5 + 1
#  x7 = rnorm(1000)
#  x8 = x7 + 1
#  
#  NNS.SD.efficient.set(cbind(x1, x2, x3, x4, x5, x6, x7, x8), degree = 1, status = FALSE)
#  [1] "x4" "x2" "x8" "x6"

