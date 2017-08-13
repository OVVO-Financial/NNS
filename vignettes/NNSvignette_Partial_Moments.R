## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(NNS)

## ----mean----------------------------------------------------------------
set.seed(123); x=rnorm(100); y=rnorm(100)

mean(x)
UPM(1,0,x)-LPM(1,0,x)

## ----variance------------------------------------------------------------
var(x)

# Sample Variance:
UPM(2,mean(x),x)+LPM(2,mean(x),x)

# Population Variance:
(UPM(2,mean(x),x)+LPM(2,mean(x),x))*(length(x)/(length(x)-1))

# Variance is also the co-variance of itself:
(Co.LPM(1,1,x,x,mean(x),mean(x))+Co.UPM(1,1,x,x,mean(x),mean(x))-D.LPM(1,1,x,x,mean(x),mean(x))-D.UPM(1,1,x,x,mean(x),mean(x)))*(length(x)/(length(x)-1))

## ----stdev---------------------------------------------------------------
sd(x)
((UPM(2,mean(x),x)+LPM(2,mean(x),x))*(length(x)/(length(x)-1)))^.5

## ----covariance----------------------------------------------------------
cov(x,y)
(Co.LPM(1,1,x,y,mean(x),mean(y))+Co.UPM(1,1,x,y,mean(x),mean(y))-D.LPM(1,1,x,y,mean(x),mean(y))-D.UPM(1,1,x,y,mean(x),mean(y)))*(length(x)/(length(x)-1))

## ----cov_dec-------------------------------------------------------------
PM.matrix(LPM.degree = 1,UPM.degree = 1,target = 'mean', variable = cbind(x,y), pop.adj = TRUE)

## ----pearson-------------------------------------------------------------
cor(x,y)
cov.xy=(Co.LPM(1,1,x,y,mean(x),mean(y))+Co.UPM(1,1,x,y,mean(x),mean(y))-D.LPM(1,1,x,y,mean(x),mean(y))-D.UPM(1,1,x,y,mean(x),mean(y)))*(length(x)/(length(x)-1))
sd.x=((UPM(2,mean(x),x)+LPM(2,mean(x),x))*(length(x)/(length(x)-1)))^.5
sd.y=((UPM(2,mean(y),y)+LPM(2,mean(y),y))*(length(y)/(length(y)-1)))^.5
cov.xy/(sd.x*sd.y)

## ----cdfs----------------------------------------------------------------
P=ecdf(x)
P(0);P(1)
LPM(0,0,x);LPM(0,1,x)

# Vectorized targets:
LPM(0,c(0,1),x)

plot(ecdf(x))
points(sort(x),LPM(0,sort(x),x),col='red')
legend('left',legend=c('ecdf','LPM.CDF'),fill=c('black','red'),border=NA,bty='n')

# Joint CDF:
Co.LPM(0,0,x,y,0,0)

# Vectorized targets:
Co.LPM(0,0,x,y,c(0,1),c(0,1))

# Continuous CDF:
plot(sort(x),LPM.ratio(1,sort(x),x),type = 'l',col='blue',lwd=3,xlab="x")

## ----pdfs----------------------------------------------------------------
NNS.PDF(degree=1,x)

## ----numerical integration-----------------------------------------------
x=seq(0,1,.001);y=x^2
UPM(1,0,y)-LPM(1,0,y)

