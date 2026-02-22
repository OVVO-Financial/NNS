## Partial Moments Equivalences

Why is it necessary to parse variance with partial moments? The additional information generated from partial moments permits a level of analysis that is not possible with traditional summary statistics alone.

Below are basic equivalences that demonstrate partial moments as the elements of variance.

### Installation
```r
install.packages("NNS")
```

### Setup
```r
library(NNS)
set.seed(123)
x <- rnorm(100)
y <- rnorm(100)
```

### Mean
A difference between the upside area and the downside area of `f(x)`.

```r
mean(x)
UPM(1, 0, x) - LPM(1, 0, x)
```

### Variance
A sum of the squared upside area and squared downside area.

```r
# Sample variance (base R)
var(x)

# Sample variance from partial moments
(UPM(2, mean(x), x) + LPM(2, mean(x), x)) * (length(x) / (length(x) - 1))

# Population adjustment of sample variance (base R)
var(x) * ((length(x) - 1) / length(x))

# Population variance from partial moments
UPM(2, mean(x), x) + LPM(2, mean(x), x)

# Variance as covariance with itself
Co.LPM(1, x, x, mean(x), mean(x)) +
  Co.UPM(1, x, x, mean(x), mean(x)) -
  D.LPM(1, 1, x, x, mean(x), mean(x)) -
  D.UPM(1, 1, x, x, mean(x), mean(x))
```

### First 4 Moments
The first four moments are returned with `NNS.moments()`. For sample statistics, set `population = FALSE`.

```r
NNS.moments(x)
NNS.moments(x, population = FALSE)
```

### Standard Deviation
```r
sd(x)
((UPM(2, mean(x), x) + LPM(2, mean(x), x)) * (length(x) / (length(x) - 1))) ^ 0.5
sqrt(NNS.moments(x, population = FALSE)$variance)
```

### Covariance
```r
cov(x, y)
(Co.LPM(1, x, y, mean(x), mean(y)) +
  Co.UPM(1, x, y, mean(x), mean(y)) -
  D.LPM(1, 1, x, y, mean(x), mean(y)) -
  D.UPM(1, 1, x, y, mean(x), mean(y))) *
  (length(x) / (length(x) - 1))
```

### Covariance Elements and Covariance Matrix
The covariance matrix $(\Sigma)$ is equal to the sum of the co-partial moments matrices less the divergent partial moments matrices.

$$\Sigma = CLPM + CUPM - DLPM - DUPM$$

```r
cov.mtx <- PM.matrix(
  LPM_degree = 1,
  UPM_degree = 1,
  target = "mean",
  variable = cbind(x, y),
  pop_adj = TRUE
)

cov.mtx

# Reassembled covariance matrix
cov.mtx$clpm + cov.mtx$cupm - cov.mtx$dlpm - cov.mtx$dupm

# Standard covariance matrix
cov(cbind(x, y))
```

### Pearson Correlation
```r
cor(x, y)

cov.xy <- (Co.LPM(1, x, y, mean(x), mean(y)) +
  Co.UPM(1, x, y, mean(x), mean(y)) -
  D.LPM(1, 1, x, y, mean(x), mean(y)) -
  D.UPM(1, 1, x, y, mean(x), mean(y))) *
  (length(x) / (length(x) - 1))

sd.x <- ((UPM(2, mean(x), x) + LPM(2, mean(x), x)) *
  (length(x) / (length(x) - 1))) ^ 0.5

sd.y <- ((UPM(2, mean(y), y) + LPM(2, mean(y), y)) *
  (length(y) / (length(y) - 1))) ^ 0.5

cov.xy / (sd.x * sd.y)
```

### Skewness
A normalized difference between upside area and downside area.

```r
PerformanceAnalytics::skewness(x)

(UPM(3, mean(x), x) - LPM(3, mean(x), x)) /
  (UPM(2, mean(x), x) + LPM(2, mean(x), x)) ^ (3 / 2)

NNS.moments(x)$skewness
```

### UPM / LPM Ratio
A more intuitive skewness measure: upside area divided by downside area.

```r
UPM(2, mean(x), x) / LPM(2, mean(x), x)
```

### Kurtosis
A normalized sum of upside area and downside area.

```r
PerformanceAnalytics::kurtosis(x)

((UPM(4, mean(x), x) + LPM(4, mean(x), x)) /
  (UPM(2, mean(x), x) + LPM(2, mean(x), x)) ^ 2) - 3

NNS.moments(x)$kurtosis
```

### CDFs
```r
P <- ecdf(x)
P(0)
P(1)

LPM(0, 0, x)
LPM(0, 1, x)

# Vectorized targets
LPM(0, c(0, 1), x)

# Joint CDF
Co.LPM(0, x, y, 0, 0)

# Vectorized targets
Co.LPM(0, x, y, c(0, 1), c(0, 1))

# Continuous CDF and variants
NNS.CDF(x, 1)
NNS.CDF(x, 1, target = mean(x))
NNS.CDF(x, 1, type = "survival")
```

### Copulas
```r
# Transform x and y to uniform marginals
u_x <- LPM.ratio(0, x, x)
u_y <- LPM.ratio(0, y, y)

# Value of copula at c(0.5, 0.5)
Co.LPM(0, u_x, u_y, 0.5, 0.5)
```

### Numerical Integration
Partial moments are asymptotic area approximations of `f(x)` akin to trapezoidal and Simpson's rules.

$$[UPM(1,0,f(x)) - LPM(1,0,f(x))] \asymp \frac{F(b)-F(a)}{b-a}$$
$$[UPM(1,0,f(x)) - LPM(1,0,f(x))] \cdot (b-a) \asymp F(b)-F(a)$$

```r
x <- seq(0, 1, 0.001)
y <- x ^ 2

(UPM(1, 0, y) - LPM(1, 0, y)) * (1 - 0)
```

For total area, not just the definite integral, sum the partial moments and multiply by `(b - a)`:

$$[UPM(1,0,f(x)) + LPM(1,0,f(x))] \cdot (b-a) \asymp \left|\int_a^b f(x)dx\right|$$

### Bayes' Theorem
For example, when ascertaining the probability of an increase in $A$ given an increase in $B$, the `Co.UPM(degree, x, y, target_x, target_y)` target parameters are set to `target_x = 0` and `target_y = 0` and the `UPM(degree, target, variable)` target parameter is also set to `target = 0`.

$$P(A|B)=\frac{Co.UPM(0,A,B,0,0)}{UPM(0,0,B)}$$

### References
- [Partial Moments as a Unifying Primitive: From Distribution and Moments to Expected Utility and Nonparametric Regression](https://doi.org/10.2139/ssrn.6249658)
- [Nonlinear Nonparametric Statistics: Using Partial Moments](https://github.com/OVVO-Financial/NNS/blob/NNS-Beta-Version/examples/index.md)
- [Cumulative Distribution Functions and UPM/LPM Analysis](https://doi.org/10.2139/ssrn.2148482)
- [Continuous CDFs and ANOVA with NNS](https://doi.org/10.2139/ssrn.3007373)
- [f(Newton)](https://doi.org/10.2139/ssrn.2186471)
- [Bayes' Theorem From Partial Moments](https://doi.org/10.2139/ssrn.3457377)
