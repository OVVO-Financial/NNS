## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
library(data.table)
data.table::setDTthreads(1L)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)

## ----setup2, message=FALSE, warning=FALSE-------------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)

## ----linear-------------------------------------------------------------------
x = seq(-5, 5, .05); y = x ^ 3

for(i in 1 : 4){NNS.part(x, y, order = i, Voronoi = TRUE, obs.req = 0)}

## ----x part,results='hide'----------------------------------------------------
for(i in 1 : 4){NNS.part(x, y, order = i, type = "XONLY", Voronoi = TRUE)}

## ----res2, echo=FALSE---------------------------------------------------------
NNS.part(x,y,order = 4, type = "XONLY")

## ----depreg},results='hide'---------------------------------------------------
for(i in 1 : 3){NNS.part(x, y, order = i, obs.req = 0, Voronoi = TRUE, type = "XONLY") ; NNS.reg(x, y, order = i, ncores = 1)}

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"------------------
NNS.reg(x, y, ncores = 1)

## ----nonlinear multi,fig.width=5,fig.height=3,fig.align = "center"------------
f = function(x, y) x ^ 3 + 3 * y - y ^ 3 - 3 * x
y = x ; z <- expand.grid(x, y)
g = f(z[ , 1], z[ , 2])
NNS.reg(z, g, order = "max", plot = FALSE, ncores = 1)

## ----nonlinear_class,fig.width=5,fig.height=3,fig.align = "center", message = FALSE----
NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", location = "topleft", ncores = 1)$equation

## ----nonlinear_class2,fig.width=5,fig.height=3,fig.align = "center", message = FALSE, echo=FALSE----
a = NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", location = "topleft", ncores = 1, plot = FALSE)$equation

## ----nonlinear class threshold,fig.width=5,fig.height=3,fig.align = "center"----
NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", threshold = .75, location = "topleft", ncores = 1)$equation

## ----nonlinear class threshold 2,fig.width=5,fig.height=3,fig.align = "center", echo=FALSE----
a = NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", threshold = .75, location = "topleft", ncores = 1, plot = FALSE)$equation

## ----final,fig.width=5,fig.height=3,fig.align = "center"----------------------
NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", threshold = .75, point.est = iris[1 : 10, 1 : 4], location = "topleft", ncores = 1)$Point.est

## ----class,fig.width=5,fig.height=3,fig.align = "center", message=FALSE-------
NNS.reg(iris[ , 1 : 4], iris[ , 5], type = "CLASS", point.est = iris[1 : 10, 1 : 4], location = "topleft", ncores = 1)$Point.est

## ----stack,fig.width=5,fig.height=3,fig.align = "center", message=FALSE, eval=FALSE----
# NNS.stack(IVs.train = iris[ , 1 : 4],
#           DV.train = iris[ , 5],
#           IVs.test = iris[1 : 10, 1 : 4],
#           dim.red.method = "cor",
#           obj.fn = expression( mean(round(predicted) == actual) ),
#           objective = "max", type = "CLASS",
#           folds = 1, ncores = 1)

## ----stackevalres, eval = FALSE-----------------------------------------------
# Folds Remaining = 0
# Current NNS.reg(... , threshold = 0.80 ) MAX Iterations Remaining = 1
# Current NNS.reg(... , threshold = 0.40 ) MAX Iterations Remaining = 0
# Current NNS.reg(. , n.best = 1 ) MAX Iterations Remaining = 12
# Current NNS.reg(. , n.best = 2 ) MAX Iterations Remaining = 11
# Current NNS.reg(. , n.best = 3 ) MAX Iterations Remaining = 10
# Current NNS.reg(. , n.best = 4 ) MAX Iterations Remaining = 9
# Current NNS.reg(. , n.best = 5 ) MAX Iterations Remaining = 8
# $OBJfn.reg
# [1] 0.9733333
# 
# $NNS.reg.n.best
# [1] 1
# 
# $probability.threshold
# [1] 0.547
# 
# $OBJfn.dim.red
# [1] 0.9666667
# 
# $NNS.dim.red.threshold
# [1] 0.8
# 
# $reg
#  [1] 1 1 1 1 1 1 1 1 1 1
# 
# $reg.pred.int
# NULL
# 
# $dim.red
#  [1] 1 1 1 1 1 1 1 1 1 1
# 
# $dim.red.pred.int
# NULL
# 
# $stack
#  [1] 1 1 1 1 1 1 1 1 1 1
# 
# $pred.int
# NULL

## ----stack2, message = FALSE,fig.width=5,fig.height=3,fig.align = "center",results='hide', eval = FALSE----
# set.seed(123)
# x = rnorm(100); y = rnorm(100)
# 
# nns.params = NNS.stack(IVs.train = cbind(x, x),
#                         DV.train = y,
#                         method = 1, ncores = 1)

## ----stack2optim, echo = FALSE------------------------------------------------
set.seed(123)
x = rnorm(100); y = rnorm(100)

nns.params = list()
nns.params$NNS.reg.n.best = 100

## ----stack2res, fig.width=5,fig.height=3,fig.align = "center",results='hide'----
NNS.reg(cbind(x, x), y, 
        n.best = nns.params$NNS.reg.n.best,
        point.est = cbind(x, x), 
        residual.plot = TRUE,  
        ncores = 1, confidence.interval = .95)

## ----smooth, fig.width=5,fig.height=3,fig.align = "center",results='hide'-----
NNS.reg(x, y, smooth = TRUE)

## ----uniimpute, eval=FALSE----------------------------------------------------
# set.seed(123)
# 
# # Univariate predictor with nonlinear signal
# n <- 400
# x <- sort(runif(n, -3, 3))
# y <- sin(x) + 0.2 * x^2 + rnorm(n, 0, 0.25)
# 
# # Induce ~25% MCAR missingness in y
# miss <- rbinom(n, 1, 0.25) == 1
# y_mis <- y
# y_mis[miss] <- NA
# 
# # ---- Increasing dimensions trick ----
# # Duplicate x so the distance operates in a 2D space: cbind(x, x).
# # This sharpens nearest-neighbor selection even in a nominally univariate setting.
# x2_train <- cbind(x[!miss], x[!miss])
# x2_miss  <- cbind(x[miss],  x[miss])
# 
# # 1-NN donor imputation with NNS.reg
# y_hat_uni <- NNS::NNS.reg(
#   x         = x2_train,             # predictors (duplicated x)
#   y         = y[!miss],             # observed responses
#   point.est = x2_miss,              # rows to impute
#   order     = "max",                # dependence-maximizing order
#   n.best    = 1,                    # 1-NN donor
#   plot      = FALSE
# )$Point.est
# 
# # Fill back
# y_completed_uni <- y_mis
# y_completed_uni[miss] <- y_hat_uni
# 
# # Plot observed vs imputed (NNS 1-NN)
# plot(x, y, pch = 1, col = "steelblue", cex = 1.5, lwd = 2,
#      xlab = "x", ylab = "y", main = "NNS 1-NN Imputation")
# points(x[miss], y_hat_uni, col = "red", pch = 15, cex = 1.3)
# 
# legend("topleft",
#        legend = c("Observed", "Imputed (NNS 1-NN)"),
#        col    = c("steelblue", "red"),
#        pch    = c(1, 15),
#        pt.lwd = c(2, NA),
#        bty    = "n")

## ----multiimpute, eval=FALSE--------------------------------------------------
# set.seed(123)
# 
# # Multivariate predictors with nonlinear & interaction structure
# n <- 800
# X <- cbind(
#   x1 = rnorm(n),
#   x2 = runif(n, -2, 2),
#   x3 = rnorm(n, 0, 1)
# )
# 
# f <- function(x1, x2, x3) 1.1*x1 - 0.8*x2 + 0.5*x3 + 0.6*x1*x2 - 0.4*x2*x3 + 0.3*sin(1.3*x1)
# y <- f(X[,1], X[,2], X[,3]) + rnorm(n, 0, 0.4)
# 
# # Induce ~30% MCAR missingness in y
# miss <- rbinom(n, 1, 0.30) == 1
# y_mis <- y
# y_mis[miss] <- NA
# 
# # Training (observed) vs rows to impute
# X_obs <- X[!miss, , drop = FALSE]
# y_obs <- y[!miss]
# X_mis <- X[ miss, , drop = FALSE]
# 
# # 1-NN donor imputation with NNS.reg
# y_hat_mv <- NNS::NNS.reg(
#   x         = X_obs,     # all observed predictors
#   y         = y_obs,     # observed responses
#   point.est = X_mis,     # rows to impute
#   order     = "max",     # dependence-maximizing order
#   n.best    = 1,         # 1-NN donor
#   plot      = FALSE
# )$Point.est
# 
# # Completed vector
# y_completed_mv <- y_mis
# y_completed_mv[miss] <- y_hat_mv
# 
# # Plot observed vs imputed (multivariate, NNS 1-NN)
# plot(seq_along(y), y,
#      pch = 1, col = "steelblue", cex = 1.5, lwd = 2,
#      xlab = "Observation index", ylab = "y",
#      main = "NNS 1-NN Multivariate Imputation")
# 
# # Overlay imputed values
# points(which(miss), y_hat_mv, pch = 15, col = "red", cex = 1.2)
# 
# # Legend
# legend("topleft",
#        legend = c("Observed", "Imputed (NNS 1-NN)"),
#        col    = c("steelblue", "red"),
#        pch    = c(1, 15),
#        pt.lwd = c(2, NA),
#        bty    = "n")

