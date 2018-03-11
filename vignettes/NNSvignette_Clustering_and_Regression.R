## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(NNS)
require(knitr)
require(rgl)

## ----linear,results='hide'-----------------------------------------------
x=seq(-5,5,.05); y=x^3

for(i in 1:4){NNS.part(x,y,order=i,noise.reduction = 'off',Voronoi = T)}

## ----x part,results='hide'-----------------------------------------------
for(i in 1:4){NNS.part(x,y,order=i,type="XONLY",Voronoi = T)}

## ----depreg,results='hide'-----------------------------------------------
for(i in 1:3){NNS.part(x,y,order=i,Voronoi = T);NNS.reg(x,y,order=i)}

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"-------------
NNS.reg(x,y,order=4,noise.reduction = 'off')

## ----nonlinear multi,results='hide',fig.width=5,fig.height=3,fig.align = "center"----
f= function(x,y) x^3+3*y-y^3-3*x
y=x; z=expand.grid(x,y)
g=f(z[,1],z[,2])
NNS.reg(z,g,order='max')

## ----iris point.est,fig.width=5,fig.height=3,fig.align = "center"--------
NNS.reg(iris[,1:4],iris[,5],point.est=iris[1:10,1:4],type="CLASS",location = 'topleft')$Point.est

## ----nonlinear class,fig.width=5,fig.height=3,fig.align = "center"-------
NNS.reg(iris[,1:4],iris[,5],dim.red.method="cor",location = 'topleft')$equation

## ----nonlinear class threshold,fig.width=5,fig.height=3,fig.align = "center"----
NNS.reg(iris[,1:4],iris[,5],dim.red.method="cor",threshold=.75,location = 'topleft')$equation

## ----final,fig.width=5,fig.height=3,fig.align = "center"-----------------
NNS.reg(iris[,1:4],iris[,5],dim.red.method="cor",threshold=.75,point.est=iris[1:10,1:4],location = 'topleft')$Point.est

