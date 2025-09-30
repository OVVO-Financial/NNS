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

## ----rhs, rows.print=18-------------------------------------------------------
NNS.reg(iris[,1:4], iris[,5], residual.plot = FALSE, ncores = 1)$rhs.partitions

## ----NNSBOOST,fig.align = "center", fig.height = 8,fig.width=6.5, eval=FALSE----
# test.set = 141:150
# 
# a = NNS.boost(IVs.train = iris[-test.set, 1:4],
#               DV.train = iris[-test.set, 5],
#               IVs.test = iris[test.set, 1:4],
#               epochs = 10, learner.trials = 10,
#               status = FALSE, balance = TRUE,
#               type = "CLASS", folds = 5)
# 
# a
# $results
#  [1] 3 3 3 3 3 3 3 3 3 3
# 
# $pred.int
# NULL
# 
# $feature.weights
#  Petal.Width Petal.Length Sepal.Length
#    0.4285714    0.4285714    0.1428571
# 
# $feature.frequency
#  Petal.Width Petal.Length Sepal.Length
#            3            3            1
# 
# mean( a$results == as.numeric(iris[test.set, 5]) )
# [1] 1

## ----NNSstack,fig.align = "center", fig.height = 8,fig.width=6.5, message=FALSE, eval= FALSE----
# b = NNS.stack(IVs.train = iris[-test.set, 1:4],
#               DV.train = iris[-test.set, 5],
#               IVs.test = iris[test.set, 1:4],
#               type = "CLASS", balance = TRUE,
#               ncores = 1, folds = 5)
# 
# b

## ----stackeval, eval = FALSE--------------------------------------------------
# $OBJfn.reg
# [1] 0.955787
# 
# $NNS.reg.n.best
# [1] 1
# 
# $probability.threshold
# [1] 0.6429167
# 
# $OBJfn.dim.red
# [1] 0.955787
# 
# $NNS.dim.red.threshold
# [1] 0.925
# 
# $reg
#  [1] 3 3 3 3 3 3 3 3 3 3
# 
# $reg.pred.int
# NULL
# 
# $dim.red
#  [1] 3 3 3 3 3 3 3 3 3 3
# 
# $dim.red.pred.int
# NULL
# 
# $stack
#  [1] 3 3 3 3 3 3 3 3 3 3
# 
# $pred.int
# NULL

## ----stackevalres, eval = FALSE-----------------------------------------------
# mean( b$stack == as.numeric(iris[test.set, 5]) )

## ----stackreseval, eval = FALSE-----------------------------------------------
# [1] 1

