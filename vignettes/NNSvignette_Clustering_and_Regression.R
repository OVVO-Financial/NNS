## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(NNS)
require(knitr)
require(rgl)

## ----linear--------------------------------------------------------------
x=seq(-5,5,.05); y=x^3

NNS.part(x,y,Voronoi = T)

## ----x part--------------------------------------------------------------
NNS.part(x,y,Voronoi = T,type="XONLY",order=3)

## ----depreg,results='hide'-----------------------------------------------
for(i in 1:3){NNS.part(x,y,order=i,Voronoi = T);NNS.reg(x,y,order=i)}

## ----nonlinear-----------------------------------------------------------
NNS.reg(x,y,order=4,noise.reduction = 'off')

## ----nonlinear multi,results='hide'--------------------------------------
f= function(x,y) x^3+3*y-y^3-3*x
y=x; z=expand.grid(x,y)
g=f(z[,1],z[,2])
NNS.reg(z,g,order='max')

## ----iris point.est------------------------------------------------------
NNS.reg(iris[,1:4],iris[,5],point.est=iris[1:10,1:4],type="CLASS")$Point.est

## ----nonlinear class-----------------------------------------------------
NNS.reg(iris[,1:4],iris[,5],dim.red=TRUE)$equation

## ----nonlinear class threshold-------------------------------------------
NNS.reg(iris[,1:4],iris[,5],dim.red=TRUE,threshold=.75)$equation

## ------------------------------------------------------------------------
NNS.reg(iris[,1:4],iris[,5],dim.red=TRUE,threshold=.75,point.est=iris[1:10,1:4])$Point.est

