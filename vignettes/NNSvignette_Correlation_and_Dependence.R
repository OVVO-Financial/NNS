## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2,message=FALSE------------------------------------------------
require(NNS)
require(knitr)
require(rgl)
require(data.table)
require(plyr)

## ----linear,fig.width=5,fig.height=3,fig.align = "center"----------------
x = seq(0, 3, .01) ; y = 2 * x

cor(x, y)
NNS.dep(x, y, print.map = TRUE, order = 3)

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"-------------
x=seq(0, 3, .01) ; y = x ^ 10

cor(x, y)
NNS.dep(x, y, print.map = TRUE, order = 3)

## ----dependence,fig.width=5,fig.height=3,fig.align = "center"------------
set.seed(123)
df <- data.frame(x = runif(10000, -1, 1), y = runif(10000, -1, 1))
df <- subset(df, (x ^ 2 + y ^ 2 <= 1 & x ^ 2 + y ^ 2 >= 0.95))
NNS.dep(df$x, df$y, print.map = TRUE)

