## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2, message=FALSE, warning=FALSE-------------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

## ----linear-------------------------------------------------------------------
x = seq(-5, 5, .05); y = x ^ 3

for(i in 1 : 4){NNS.part(x, y, order = i, Voronoi = TRUE, obs.req = 0)}

## ----x part,results='hide'----------------------------------------------------
for(i in 1 : 4){NNS.part(x, y, order = i, type = "XONLY", Voronoi = TRUE)}

## ----res2, echo=FALSE---------------------------------------------------------
NNS.part(x,y,order = 4, type = "XONLY")

## ----depreg},results='hide'---------------------------------------------------
for(i in 1 : 3){NNS.part(x, y, order = i, obs.req = 0, Voronoi = TRUE) ; NNS.reg(x, y, order = i, ncores = 1)}

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"------------------
NNS.reg(x, y, ncores = 1)

## ----nonlinear multi,fig.width=5,fig.height=3,fig.align = "center"------------
f= function(x, y) x ^ 3 + 3 * y - y ^ 3 - 3 * x
y = x ; z = expand.grid(x, y)
g = f(z[ , 1], z[ , 2])
NNS.reg(z, g, order = "max", ncores = 1)

## ----nonlinear_class,fig.width=5,fig.height=3,fig.align = "center", message = FALSE----
NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", location = "topleft", ncores = 1)$equation

## ----nonlinear_class2,fig.width=5,fig.height=3,fig.align = "center", message = FALSE, echo=FALSE----
a=NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", location = "topleft", ncores = 1, plot = FALSE)$equation

## ----nonlinear class threshold,fig.width=5,fig.height=3,fig.align = "center"----
NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", threshold = .75, location = "topleft", ncores = 1)$equation

## ----nonlinear class threshold 2,fig.width=5,fig.height=3,fig.align = "center", echo=FALSE----
a=NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", threshold = .75, location = "topleft", ncores = 1, plot = FALSE)$equation

## ----final,fig.width=5,fig.height=3,fig.align = "center"----------------------
NNS.reg(iris[ , 1 : 4], iris[ , 5], dim.red.method = "cor", threshold = .75, point.est = iris[1 : 10, 1 : 4], location = "topleft", ncores = 1)$Point.est

## ----class,fig.width=5,fig.height=3,fig.align = "center", message=FALSE-------
NNS.reg(iris[ , 1 : 4], iris[ , 5], type = "CLASS", point.est = iris[1 : 10, 1 : 4], location = "topleft", ncores = 1)$Point.est

## ----stack,fig.width=5,fig.height=3,fig.align = "center", message=FALSE-------
NNS.stack(IVs.train = iris[ , 1 : 4], 
          DV.train = iris[ , 5], 
          IVs.test = iris[1 : 10, 1 : 4],
          dim.red.method = "cor",
          obj.fn = expression( mean(round(predicted) == actual) ),
          objective = "max", type = "CLASS", 
          folds = 1, ncores = 1)

## ----stack2, message = FALSE,fig.width=5,fig.height=3,fig.align = "center",results='hide'----
set.seed(123)
x <- rnorm(100); y <- rnorm(100)

nns.params <- NNS.stack(IVs.train = cbind(x, x),
                        DV.train = y,
                        method = 1, ncores = 1)

NNS.reg(cbind(x, x), y, 
        n.best = nns.params$NNS.reg.n.best,
        point.est = cbind(x, x), ncores = 1)

