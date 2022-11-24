## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2,message=FALSE,warning = FALSE-------------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

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

## ----dependence,fig.width=5,fig.height=3,fig.align = "center"-----------------
set.seed(123)
df <- data.frame(x = runif(10000, -1, 1), y = runif(10000, -1, 1))
df <- subset(df, (x ^ 2 + y ^ 2 <= 1 & x ^ 2 + y ^ 2 >= 0.95))

## ----circle1,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(df$x, df$y, Voronoi = TRUE, order = 3, obs.req = 0)

## ----res3---------------------------------------------------------------------
NNS.dep(df$x, df$y)

## ----permutations-------------------------------------------------------------
## p-values for [NNS.dep]
x <- seq(-5, 5, .1); y <- x^2 + rnorm(length(x))

## ----perm1,fig.width=5,fig.height=3,fig.align = "center", results='hide', echo=FALSE----
NNS.part(x, y, Voronoi = TRUE, order = 3)

## ----permutattions_res,fig.width=5,fig.height=3,fig.align = "center"----------
NNS.dep(x, y, p.value = TRUE, print.map = TRUE)

## ----multi, warning=FALSE-----------------------------------------------------
set.seed(123)
x <- rnorm(1000); y <- rnorm(1000); z <- rnorm(1000)
NNS.copula(cbind(x, y, z), plot = TRUE, independence.overlay = TRUE)

## ----multisim-----------------------------------------------------------------
# Add variable x to original data to avoid total independence (example only)
original.data <- cbind(x, y, z, x)

# Determine dependence structure
dep.structure <- apply(original.data, 2, function(x) LPM.ratio(1, x, x))
  
# Generate new data with different mean and sd (or distribution type or dimensions)
new.data <- sapply(1:ncol(original.data), function(x) rnorm(dim(original.data)[1], mean = 10, sd = 20))

# Apply dependence structure to new data
new.dep.data <- sapply(1:ncol(original.data), function(x) LPM.VaR(dep.structure[,x], 1, new.data[,x]))

## ----comparison, warning=FALSE------------------------------------------------
NNS.copula(original.data)
NNS.copula(new.dep.data)

