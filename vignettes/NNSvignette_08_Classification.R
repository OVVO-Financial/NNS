## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)


## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
require(knitr)
require(rgl)


## ----rhs, rows.print=18-------------------------------------------------------
NNS.reg(iris[,1:4], iris[,5], residual.plot = FALSE, ncores = 1)$rhs.partitions


## ----NNSBOOST,fig.align = "center", fig.height = 8,fig.width=6.5, eval=FALSE----
## set.seed(20261103)
## test.set <- sort(unlist(lapply(split(seq_len(nrow(iris)), iris$Species),
##                                function(rows) sample(rows, 4))))
##
## a = NNS.boost(IVs.train = iris[-test.set, 1:4],
##               DV.train = iris[-test.set, 5],
##               IVs.test = iris[test.set, 1:4],
##               epochs = 10, learner.trials = 10,
##               status = FALSE, balance = TRUE,
##               type = "CLASS")
##
## a
## mean(a$results == as.numeric(iris[test.set, 5]))


## ----NNSstack,fig.align = "center", fig.height = 8,fig.width=6.5, message=FALSE, eval= FALSE----
## b = NNS.stack(IVs.train = iris[-test.set, 1:4],
##               DV.train = iris[-test.set, 5],
##               IVs.test = iris[test.set, 1:4],
##               type = "CLASS", balance = TRUE,
##               ncores = 1, folds = 5)
##
## b


## ----stackeval, eval = FALSE--------------------------------------------------
## mean(b$stack == as.numeric(iris[test.set, 5]))

