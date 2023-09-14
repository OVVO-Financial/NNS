## Partial Moments Equivalences
Below are some basic equivalences demonstrating partial moments' role as the elements of variance.

### Why is this relevant?  
The additional information generated from partial moments permits a level of analysis simply not possible with traditional summary statistics.
There is further introductory material on partial moments and their extension into nonlinear analysis & behavioral finance applications available at:

https://www.linkedin.com/pulse/elements-variance-fred-viole


## Installation
```r
require(devtools); install_github('OVVO-Financial/NNS', ref = "NNS-Beta-Version")
```

### Mean
A difference between the upside area and the downside area of f(x).
```r
set.seed(123); x = rnorm(100); y = rnorm(100)

> mean(x)
[1] 0.09040591
> UPM(1,0,x)-LPM(1,0,x)
[1] 0.09040591
```
### Variance
A sum of the squared upside area and the squared downside area.
```r
> var(x)
[1] 0.8332328
# Sample Variance:
> UPM(2,mean(x),x)+LPM(2,mean(x),x)
[1] 0.8249005
# Population Variance:
> (UPM(2,mean(x),x)+LPM(2,mean(x),x))*(length(x)/(length(x)-1))
[1] 0.8332328
# Variance is also the co-variance of itself:
> (Co.LPM(1,x,x,mean(x),mean(x))+Co.UPM(1,x,x,mean(x),mean(x))-D.LPM(1,1,x,x,mean(x),mean(x))-D.UPM(1,1,x,x,mean(x),mean(x)))*(length(x)/(length(x)-1))
[1] 0.8332328
```

### The first 4 moments are returned with the function `NNS.moments`. For sample statistics, set `population = FALSE`.
```r
> NNS.moments(x)
$mean
[1] 0.09040591

$variance
[1] 0.8332328

$skewness
[1] 0.06049948
 
$kurtosis
[1] -0.161053

> NNS.moments(x, population = FALSE)
$mean
[1] 0.09040591

$variance
[1] 0.8249005
 
$skewness
[1] 0.06235774

$kurtosis
[1] -0.1069186
```

### Standard Deviation
```r
> sd(x)
[1] 0.9128159
> ((UPM(2,mean(x),x)+LPM(2,mean(x),x))*(length(x)/(length(x)-1)))^.5
[1] 0.9128159
```
### Covariance
```r
> cov(x,y)
[1] -0.04372107
> (Co.LPM(1,x,y,mean(x),mean(y))+Co.UPM(1,x,y,mean(x),mean(y))-D.LPM(1,1,x,y,mean(x),mean(y))-D.UPM(1,1,x,y,mean(x),mean(y)))*(length(x)/(length(x)-1))
[1] -0.04372107
```
### Covariance Elements and Covariance Matrix
The covariance matrix $(\Sigma)$ is equal to the sum of the co-partial moments matrices less the divergent partial moments matrices.

$$\Sigma = CLPM + CUPM - DLPM - DUPM $$


```r
> cov.mtx = PM.matrix(LPM_degree = 1, UPM_degree = 1, target = 'mean', variable = cbind(x,y), pop_adj = TRUE)
> cov.mtx
$cupm
          x         y
x 0.4299250 0.1033601
y 0.1033601 0.5411626

$dupm
          x         y
x 0.0000000 0.1469182
y 0.1560924 0.0000000

$dlpm
          x         y
x 0.0000000 0.1560924
y 0.1469182 0.0000000

$clpm
          x         y
x 0.4033078 0.1559295
y 0.1559295 0.3939005

$cov.matrix
            x           y
x  0.83323283 -0.04372107
y -0.04372107  0.93506310


# Reassembled Covariance Matrix
> cov.mtx$cupm + cov.mtx$clpm - cov.mtx$dupm - cov.mtx$dlpm
            x           y
x  0.83323283 -0.04372107
y -0.04372107  0.93506310


# Standard Covariance Matrix
> cov(cbind(x,y))
            x           y
x  0.83323283 -0.04372107
y -0.04372107  0.93506310
```


### Pearson Correlation
```r
> cor(x,y)
[1] -0.04953215
> cov.xy = (Co.LPM(1,x,y,mean(x),mean(y))+Co.UPM(1,x,y,mean(x),mean(y))-D.LPM(1,1,x,y,mean(x),mean(y))-D.UPM(1,1,x,y,mean(x),mean(y)))*(length(x)/(length(x)-1))
> sd.x = ((UPM(2,mean(x),x)+LPM(2,mean(x),x))*(length(x)/(length(x)-1)))^.5
> sd.y = ((UPM(2,mean(y),y)+LPM(2,mean(y),y))*(length(y)/(length(y)-1)))^.5
> cov.xy/(sd.x*sd.y)
[1] -0.04953215
```
### Skewness
A normalized difference between upside area and downside area.
```r
> library(PerformanceAnalytics)
> skewness(x)
[1] 0.06049948
> ((UPM(3,mean(x),x)-LPM(3,mean(x),x))/(UPM(2,mean(x),x)+LPM(2,mean(x),x))^(3/2))
[1] 0.06049948
```
### UPM/LPM - a more intuitive measure of skewness.  (Upside area / Downside area)
```r
> UPM(2,mean(x),x)/LPM(2,mean(x),x)
[1] 1.065997
```
### Kurtosis
A normalized sum of upside area and downside area.
```r
> library(PerformanceAnalytics)
> kurtosis(x)
[1] -0.161053
> ((UPM(4,mean(x),x)+LPM(4,mean(x),x))/(UPM(2,mean(x),x)+LPM(2,mean(x),x))^2)-3
[1] -0.161053
```
### CDFs
```r
> P = ecdf(x)
> P(0); P(1)
[1] 0.48
[1] 0.83
> LPM(0,0,x); LPM(0,1,x)
[1] 0.48
[1] 0.83
# Vectorized targets:
> LPM(0,c(0,1),x)
[1] 0.48 0.83
# Joint CDF:
> Co.LPM(0,x,y,0,0)
[1] 0.28
# Vectorized targets:
> Co.LPM(0,x,y,c(0,1),c(0,1))
[1] 0.28 0.73

# Alternatively via NNS.CDF()
> NNS.CDF(x)
```
### Copulas
```r
# Transform x and y so that they are uniform
u_x = LPM.ratio(0, x, x)
u_y = LPM.ratio(0, y, y)

# Value of copula at c(.5, .5)
Co.LPM(0, u_x, u_y, .5, .5)
[1] 0.26
```
### Numerical Integration - [UPM(1,0,f(x))-LPM(1,0,f(x))]=[F(b)-F(a)]/[b-a]
```r
# x is uniform sample over interval [a,b]; y = f(x)
> x = seq(0,1,.001); y = x^2
> UPM(1,0,y)-LPM(1,0,y)
[1] 0.3335
```

### Bayes' Theorem
```r
See the following example explaining Bayes' Theorem and partial moments: 
```
https://github.com/OVVO-Financial/NNS/blob/NNS-Beta-Version/examples/Bayes'%20Theorem%20From%20Partial%20Moments.pdf
