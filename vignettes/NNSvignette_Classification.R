## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

## ----rhs, rows.print=18-------------------------------------------------------
NNS.reg(iris[,1:4], iris[,5], residual.plot = FALSE, ncores = 1)$rhs.partitions

## ----NNSBOOST,fig.align = "center", fig.height = 8,fig.width=6.5--------------
set.seed(1234)
test.set = sample(150,10)
 
a = NNS.boost(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              epochs = 10, learner.trials = 10, 
              status = FALSE, balance = TRUE,
              type = "CLASS", folds = 1)

a$results

a$feature.weights

mean( a$results == as.numeric(iris[test.set, 5]) )

## ----NNSstack,fig.align = "center", fig.height = 8,fig.width=6.5,message=FALSE----
b = NNS.stack(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              type = "CLASS", balance = TRUE,
              ncores = 1, folds = 1)

b

mean( b$stack == as.numeric(iris[test.set, 5]) )

