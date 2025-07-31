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

## ----linear,fig.width=5,fig.height=3,fig.align = "center"---------------------
x = seq(0, 3, .01) ; y = 2 * x

## ----linear1,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(x, y, Voronoi = TRUE, order = 3)

## ----res1---------------------------------------------------------------------
cor(x, y)
NNS.dep(x, y)

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center", results='hide'----
x = seq(0, 3, .01) ; y = x ^ 10

## ----nonlinear1,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(x, y, Voronoi = TRUE, order = 3)

## ----res2a--------------------------------------------------------------------
cor(x, y)
NNS.dep(x, y)

## ----nonlinear_sin,fig.width=5,fig.height=3,fig.align = "center", results='hide'----
x = seq(0, 12*pi, pi/100) ; y = sin(x)

## ----nonlinear1_sin,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(x, y, Voronoi = TRUE, order = 3, obs.req = 0)

## ----res2_sin-----------------------------------------------------------------
cor(x, y)
NNS.dep(x, y)

## ----asym1--------------------------------------------------------------------
cor(x, y)
NNS.dep(x, y, asym = TRUE)

## ----asym2--------------------------------------------------------------------
cor(y, x)
NNS.dep(y, x, asym = TRUE)

## ----dependence,fig.width=5,fig.height=3,fig.align = "center"-----------------
set.seed(123)
df = data.frame(x = runif(10000, -1, 1), y = runif(10000, -1, 1))
df = subset(df, (x ^ 2 + y ^ 2 <= 1 & x ^ 2 + y ^ 2 >= 0.95))

## ----circle1,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(df$x, df$y, Voronoi = TRUE, order = 3, obs.req = 0)

## ----res3---------------------------------------------------------------------
NNS.dep(df$x, df$y)

## ----permutations-------------------------------------------------------------
## p-values for [NNS.dep]
set.seed(123)
x = seq(-5, 5, .1); y = x^2 + rnorm(length(x))

## ----perm1,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(x, y, Voronoi = TRUE, order = 3)

## ----permutattions_res,fig.width=5,fig.height=3,fig.align = "center"----------
NNS.dep(x, y, p.value = TRUE, print.map = TRUE)

## ----multi, warning=FALSE-----------------------------------------------------
set.seed(123)
x = rnorm(1000); y = rnorm(1000); z = rnorm(1000)
NNS.copula(cbind(x, y, z), plot = TRUE, independence.overlay = TRUE)

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")

