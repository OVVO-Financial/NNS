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

## ----rhs, rows.print=18-------------------------------------------------------
NNS.reg(iris[,1:4], iris[,5], residual.plot = FALSE, ncores = 1)$rhs.partitions

## ----NNSBOOST,fig.align = "center", fig.height = 8,fig.width=6.5, eval=FALSE----
#  test.set = 141:150
#  
#  a = NNS.boost(IVs.train = iris[-test.set, 1:4],
#                DV.train = iris[-test.set, 5],
#                IVs.test = iris[test.set, 1:4],
#                epochs = 10, learner.trials = 10,
#                status = FALSE, balance = TRUE,
#                type = "CLASS", folds = 1)
#  
#  a$results
#  [1] 3 3 3 3 3 3 3 3 3 3
#  
#  a$feature.weights
#  Petal.Width Petal.Length Sepal.Length
#     0.4285714    0.4285714    0.1428571
#  
#  mean( a$results == as.numeric(iris[test.set, 5]) )
#  [1] 1

## ----NNSstack,fig.align = "center", fig.height = 8,fig.width=6.5, message=FALSE, eval= FALSE----
#  b = NNS.stack(IVs.train = iris[-test.set, 1:4],
#                DV.train = iris[-test.set, 5],
#                IVs.test = iris[test.set, 1:4],
#                type = "CLASS", balance = TRUE,
#                ncores = 1, folds = 1)
#  
#  b

## ----stackeval, eval = FALSE--------------------------------------------------
#  $OBJfn.reg
#  [1] 1
#  
#  $NNS.reg.n.best
#  [1] 1
#  
#  $probability.threshold
#  [1] 0.43875
#  
#  $OBJfn.dim.red
#  [1] 0.9798658
#  
#  $NNS.dim.red.threshold
#  [1] 0.93
#  
#  $reg
#   [1] 3 3 3 3 3 3 3 3 3 3
#  
#  $reg.pred.int
#  NULL
#  
#  $dim.red
#   [1] 3 3 3 3 3 3 3 3 3 3
#  
#  $dim.red.pred.int
#  NULL
#  
#  $stack
#   [1] 3 3 3 3 3 3 3 3 3 3
#  
#  $pred.int
#  NULL

## ----stackevalres, eval = FALSE-----------------------------------------------
#  mean( b$stack == as.numeric(iris[test.set, 5]) )

## ----stackreseval, eval = FALSE-----------------------------------------------
#  [1] 1

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

