# Getting Started with NNS: Overview

``` r

# Prereqs (uncomment if needed):
# install.packages("NNS")

library(NNS)
```

## Orientation

**Goal.** A complete, hands‑on curriculum for Nonlinear Nonparametric
Statistics (NNS) using **partial moments**. Each section blends
narrative intuition, precise math, and executable code.

**Structure.** 1. Foundations — partial moments & variance decomposition
2. Descriptive & distributional tools 3. Dependence & nonlinear
association 4. Normalization & Rescaling 5. Hypothesis testing, ANOVA &
Stochastic Superiority 6. Regression, boosting, stacking & causality 7.
Time series & forecasting 8. Simulation (max‑entropy) & Monte Carlo 9.
Portfolio & stochastic dominance

**Notation.** For a random variable $`X`$ and threshold/target $`t`$,
the population $`n`$‑th **partial moments** are defined as:

``` math
\operatorname{LPM}(n,t,X) 
= \int_{-\infty}^{t} (t-x)^{n} \, dF_X(x),
\qquad
\operatorname{UPM}(n,t,X) 
= \int_{t}^{\infty} (x-t)^{n} \, dF_X(x).
```

The **empirical** estimators replace $`F_X`$ with the empirical CDF
$`\hat F_n`$ (or, equivalently, use indicator functions):

``` math
\widehat{\operatorname{LPM}}_n(t;X) = \frac{1}{n} \sum_{i=1}^n (t-x_i)^n \, \mathbf{1}_{\{x_i \le t\}},
\qquad
\widehat{\operatorname{UPM}}_n(t;X) = \frac{1}{n} \sum_{i=1}^n (x_i-t)^n \, \mathbf{1}_{\{x_i > t\}}.
```

These correspond to integrals over the measurable subsets
$`\{X \le t\}`$ and $`\{X > t\}`$ in a $`\sigma`$‑algebra; the empirical
sums are discrete analogues of Lebesgue integrals.

------------------------------------------------------------------------

## 1. Foundations — Partial Moments & Variance Decomposition

### 1.1 Why partial moments

- Classical variance treats upside and downside symmetrically. Partial
  moments separate them, allowing **asymmetric risk/reward** analysis
  around a chosen target $`t`$ (often the mean or a benchmark).
- At $`t=\mu_X`$:
  ``` math
  \operatorname{Var}(X) = \operatorname{UPM}(2,\mu_X,X) + \operatorname{LPM}(2,\mu_X,X)\quad\text{(exact empirical identity)}.
  ```
  This **is not** the same as splitting conditional variances around a
  threshold; partial moments use a *global* reference, preserving the
  between‑group contribution.

### 1.2 Core functions and headers

- `LPM(degree, target, variable)`
- `UPM(degree, target, variable)`

### 1.3 Code: variance decomposition & CDF

``` r

set.seed(42)

# Normal sample
y <- rnorm(3000)
mu <- mean(y)
L2 <- LPM(2, mu, y); U2 <- UPM(2, mu, y)
cat(sprintf("LPM2 + UPM2 = %.6f vs var(y)=%.6f\n", (L2+U2)*(length(y) / (length(y) - 1)), var(y)))
```

    ## LPM2 + UPM2 = 1.011889 vs var(y)=1.011889

``` r

# Empirical CDF via LPM.ratio(0, t, x)
for (t in c(-1,0,1)) {
  cdf_lpm <- LPM.ratio(0, t, y)
  cat(sprintf("CDF at t=%+.1f : LPM.ratio=%.4f | empirical=%.4f\n", t, cdf_lpm, mean(y<=t)))
}
```

    ## CDF at t=-1.0 : LPM.ratio=0.1633 | empirical=0.1633
    ## CDF at t=+0.0 : LPM.ratio=0.5043 | empirical=0.5043
    ## CDF at t=+1.0 : LPM.ratio=0.8480 | empirical=0.8480

``` r

# Asymmetry on a skewed distribution
z <- rexp(3000)-1; mu_z <- mean(z)
cat(sprintf("Skewed z: LPM2=%.4f, UPM2=%.4f (expect imbalance)\n", LPM(2,mu_z,z), UPM(2,mu_z,z)))
```

    ## Skewed z: LPM2=0.2780, UPM2=0.7682 (expect imbalance)

**Interpretation.** The equality `LPM2 + UPM2 == var(x)` (Bessel
adjustment used) holds because deviations are measured against the
*global* mean. `LPM.ratio(0, t, x)` constructs an empirical CDF directly
from partial‑moment counts.

------------------------------------------------------------------------

## 2. Descriptive & Distributional Tools

### 2.1 Higher moments from partial moments

Define asymmetric analogues of skewness/kurtosis using
$`\operatorname{UPM}_3`$, $`\operatorname{LPM}_3`$ (and degree 4),
yielding robust tail diagnostics without parametric assumptions.

**Header.**

- `NNS.moments(x)`

``` r

M <- NNS.moments(y)
M
```

    ## $mean
    ## [1] -0.0114498
    ## 
    ## $variance
    ## [1] 1.011552
    ## 
    ## $skewness
    ## [1] -0.007412142
    ## 
    ## $kurtosis
    ## [1] 0.06723772

### 2.2 Mode estimation (no bin‑or‑bandwidth angst)

**Header.**

- `NNS.mode(x)`

``` r

set.seed(23)
multimodal <- c(rnorm(1500,-2,.5), rnorm(1500,2,.5))
NNS.mode(multimodal,multi = TRUE)
```

    ## [1] -2.012070  2.024336

### 2.3 CDF tables via LPM ratios

**Headers.**

- `LPM.ratio(degree = 0, target, variable)` (empirical CDF when
  `degree=0`)
- `UPM.ratio(degree = 0, target, variable)`
- `LPM.VaR(p, degree, variable)` (quantiles via partial‑moment CDFs)
- `UPM.VaR(p, degree, variable)`

``` r

qgrid <- LPM.VaR(seq(0.05,0.95,.1),0,z) # equivalent to quantile(z,probs = seq(0.05,0.95,by=0.1))
CDF_tbl <- data.frame(threshold = as.numeric(qgrid), CDF = LPM.ratio(0,qgrid,z))
CDF_tbl
```

    ##      threshold  CDF
    ## 1  -0.94052127 0.05
    ## 2  -0.83748109 0.15
    ## 3  -0.71317882 0.25
    ## 4  -0.57443327 0.35
    ## 5  -0.41017671 0.45
    ## 6  -0.20424962 0.55
    ## 7   0.06850182 0.65
    ## 8   0.41462712 0.75
    ## 9   0.94307172 0.85
    ## 10  2.09633977 0.95

------------------------------------------------------------------------

## 3. Dependence & Nonlinear Association

### 3.1 Why move beyond Pearson $`r`$

Pearson captures linear monotone relationships. Many structures
(U‑shapes, saturation, asymmetric tails) produce near‑zero $`r`$ despite
strong dependence. Partial‑moment dependence metrics respond to such
structure.

**Headers.**

- `Co.LPM(degree_lpm, x, y, target_x, target_y, degree_y)` /
  `Co.UPM(...)` (co‑partial moments)
- `PM.matrix(LPM_degree, UPM_degree, target=NULL, variable, pop_adj=TRUE)`
- `NNS.dep(x, y)` (scalar dependence coefficient)
- `NNS.copula(X, target=NULL, continuous=TRUE, plot=FALSE, independence.overlay=FALSE)`

### 3.2 Code: nonlinear dependence

``` r

set.seed(1)
x <- runif(2000,-1,1)
y <- x^2 + rnorm(2000, sd=.05)
cat(sprintf("Pearson r = %.4f\n", cor(x,y)))
```

    ## Pearson r = 0.0006

``` r

cat(sprintf("NNS.dep  = %.4f\n", NNS.dep(x,y)$Dependence))
```

    ## NNS.dep  = 0.9073

``` r

X <- data.frame(a=x, b=y, c=x*y + rnorm(2000, sd=.05))
pm <- PM.matrix(1, 1, target = "means", variable=X, pop_adj=TRUE)
pm
```

    ## $cupm
    ##            a          b          c
    ## a 0.17384174 0.05668152 0.10450858
    ## b 0.05668152 0.05566363 0.04414923
    ## c 0.10450858 0.04414923 0.07529373
    ## 
    ## $dupm
    ##              a          b            c
    ## a 0.0000000000 0.05675501 0.0005598221
    ## b 0.0143108307 0.00000000 0.0036839026
    ## c 0.0004239566 0.04430691 0.0000000000
    ## 
    ## $dlpm
    ##              a           b            c
    ## a 0.0000000000 0.014310831 0.0004239566
    ## b 0.0567550147 0.000000000 0.0443069142
    ## c 0.0005598221 0.003683903 0.0000000000
    ## 
    ## $clpm
    ##            a           b           c
    ## a 0.16803827 0.014485430 0.102709867
    ## b 0.01448543 0.037120650 0.003051617
    ## c 0.10270987 0.003051617 0.074865823
    ## 
    ## $cov.matrix
    ##              a             b            c
    ## a 0.3418800141  0.0001011068  0.206234664
    ## b 0.0001011068  0.0927842833 -0.000789973
    ## c 0.2062346637 -0.0007899730  0.150159552

``` r

cop <- NNS.copula(X, continuous=TRUE, plot=FALSE)
cop
```

    ## [1] 0.5217301

### 3.3 Code: copula

``` r

# Data
set.seed(123); x = rnorm(100); y = rnorm(100); z = expand.grid(x, y)

# Plot
rgl::plot3d(z[,1], z[,2], Co.LPM(0, z[,1], z[,2], z[,1], z[,2]), col = "red")

# Uniform values
u_x = LPM.ratio(0, x, x); u_y = LPM.ratio(0, y, y); z = expand.grid(u_x, u_y)

# Plot
rgl::plot3d(z[,1], z[,2], Co.LPM(0, z[,1], z[,2], z[,1], z[,2]), col = "blue")
```

**Interpretation.** `NNS.dep` remains high for curved relationships;
`PM.matrix` collects co‑partial moments across variables; `NNS.copula`
summarizes higher‑dimensional dependence using partial‑moment ratios.
Copulas are returned and evaluated via `Co.LPM` functions.

------------------------------------------------------------------------

## 4. Normalization and Rescaling

NNS provides two main tools for scaling data while preserving rank
structure and distributional shape. Both operate via deterministic
affine transformations.

### 4.1 Normalization

[`NNS.norm()`](https://OVVO-Financial.github.io/NNS/reference/NNS.norm.md)
rescales variables to a common magnitude while preserving distributional
structure. The method can be **linear** (all variables forced to have
the same mean) or **nonlinear** (using dependence weights to produce a
more nuanced scaling). In the nonlinear case, the degree of association
between variables influences the final normalized values.

**Header.**

- `NNS.norm(x, linear=TRUE, chart.type = NULL)`

``` r

# NNS.norm scales by ratios of column means, so the columns need positive means.
A <- rnorm(100, mean = 5, sd = 1)
B <- rnorm(100, mean = 20, sd = 5)
C <- rnorm(100, mean = 10, sd = 1)
D <- rnorm(100, mean = 50, sd = 10)

X <- data.frame(A, B, C, D)

# Linear scaling
lin_norm <- NNS.norm(X, linear = TRUE, chart.type=NULL, location=NULL)
```

**Interpretation.**
[`NNS.norm()`](https://OVVO-Financial.github.io/NNS/reference/NNS.norm.md)
brings variables to a common scale without distorting their
distributional shape. Linear mode equalizes means; nonlinear mode
additionally weights each variable by its dependence with others, so
more correlated variables exert greater influence on the final scaling.

### 4.2 Risk‑neutral rescale (pricing context)

[`NNS.rescale()`](https://OVVO-Financial.github.io/NNS/reference/NNS.rescale.md)
performs one‑dimensional affine transformations.

**Header.**

- `NNS.rescale(x, a, b, method=c("minmax","riskneutral"), T=NULL, type=c("Terminal","Discounted"))`

``` r

px <- 100 + cumsum(rnorm(260, sd = 1))
rn <- NNS.rescale(px, a=100, b=0.03, method="riskneutral", T=1, type="Terminal")
c( target = 100*exp(0.03*1), mean_rn = mean(rn) )
```

    ##   target  mean_rn 
    ## 103.0455 103.0455

**Interpretation.** `riskneutral` shifts the mean to match
$`S_0 e^{rT}`$ (Terminal) or $`S_0`$ (Discounted), preserving
distributional shape.

------------------------------------------------------------------------

## 5. Hypothesis Testing, ANOVA & Stochastic Superiority

### 5.1 Concept

Instead of distributional assumptions, compare groups via **LPM‑based
CDFs**. Output is a *degree of certainty* (not a p‑value) for equality
of populations or means.

**Header.**

- `NNS.ANOVA(control, treatment, means.only=FALSE, medians=FALSE, confidence.interval=.95, tails=c("Both","left","right"), pairwise=FALSE, plot=TRUE, robust=FALSE)`
- `NNS.SS(x, y, ...)`

### 5.2 Code: two‑sample & multi‑group

``` r

ctrl <- rnorm(200, 0, 1)
trt  <- rnorm(180, 0.35, 1.2)
NNS.ANOVA(control=ctrl, treatment=trt, means.only=FALSE, plot=FALSE)
```

    ## $Control
    ## [1] 0.05568255
    ## 
    ## $Treatment
    ## [1] 0.2771257
    ## 
    ## $Grand_Statistic
    ## [1] 0.1605767
    ## 
    ## $Control_CDF
    ## [1] 0.5670595
    ## 
    ## $Treatment_CDF
    ## [1] 0.4385169
    ## 
    ## $Certainty
    ## [1] 0.6905098
    ## 
    ## $Effect_Size_LB
    ##        2.5% 
    ## -0.07055716 
    ## 
    ## $Effect_Size_UB
    ##     97.5% 
    ## 0.5317766 
    ## 
    ## $Confidence_Level
    ## [1] 0.95

``` r

A <- list(g1=rnorm(150,0.0,1.1), g2=rnorm(150,0.2,1.0), g3=rnorm(150,-0.1,0.9))
NNS.ANOVA(control=A, means.only=TRUE, plot=FALSE)
```

    ## Certainty 
    ## 0.6876008

**Math sketch.** For each quantile/threshold $`t`$, compare CDFs built
from `LPM.ratio(0, t, •)` (possibly with one‑sided tails). Aggregate
across $`t`$ to a certainty score.

### 5.3 Stochastic Superiority

Stochastic superiority asks a different question than equality of means
or equality of distributions. Rather than testing whether two samples
came from the same population, or whether they share the same mean or
median, stochastic superiority measures the probability that a random
draw from one distribution exceeds a random draw from another.

For two random variables $`X`$ and $`Y`$, the stochastic superiority
probability is:

``` math
P(X > Y)
```

and with ties accounted for, the tie-adjusted stochastic superiority
measure is:

``` math
P^* = P(X > Y) + \frac{1}{2} P(X = Y)
```

A value of $`P^* = 0.5`$ indicates no directional advantage, values
above $`0.5`$ favor $`X`$, and values below $`0.5`$ favor $`Y`$.

This differs from stochastic dominance. Stochastic superiority is a
pairwise exceedance probability, while stochastic dominance requires one
distribution to be preferred to another over the entire shared support.

Below is an example comparing two distributions with unequal means.

``` r

set.seed(123)
x = rnorm(1000, mean = 0, sd = 1)
y = rnorm(1000, mean = 1, sd = 1)

NNS.SS(x, y)
```

    ## $p_gt
    ## [1] 0.233915
    ## 
    ## $p_tie
    ## [1] 0
    ## 
    ## $p_star
    ## [1] 0.233915

Since $`y`$ was generated with a higher mean, the stochastic superiority
probability for $`x`$ relative to $`y`$ should be less than $`0.5`$,
indicating that a draw from $`x`$ is less likely to exceed a draw from
$`y`$.

We can also obtain confidence intervals for the tie-adjusted superiority
probability using maximum entropy bootstrap replicates.

``` r
NNS.SS(x, y, confidence.interval = TRUE, reps = 999, ci = 0.95)[1:5]

$p_gt
[1] 0.233915

$p_tie
[1] 0

$p_star
[1] 0.233915

$lower
[1] 0.2105631

$upper
[1] 0.2537789
```

This provides an interpretable effect size for directional comparison
between two distributions without requiring identical distributions or
equal variances.

For discrete variables, ties may occur with positive probability, and
the reported `p_tie` and `p_star` values reflect that adjustment
explicitly.

``` r

set.seed(123)
x = sample(1:5, 100, replace = TRUE)
y = sample(1:5, 100, replace = TRUE)

NNS.SS(x, y)
```

    ## $p_gt
    ## [1] 0.3982
    ## 
    ## $p_tie
    ## [1] 0.1992
    ## 
    ## $p_star
    ## [1] 0.4978

------------------------------------------------------------------------

## 6. Regression, Boosting, Stacking & Causality

### 6.1 Philosophy

`NNS.reg` learns **partitioned** relationships using partial‑moment
weights — linear where appropriate, nonlinear where needed — avoiding
fragile global parametric forms.

**Headers.**

- `NNS.reg(x, y, order=NULL, smooth=TRUE, ncores=1, ...)` →
  `$Fitted.xy`, `$Point.est`, …
- `NNS.boost(IVs.train, DV.train, IVs.test, epochs, learner.trials, status, balance, type, folds)`
- `NNS.stack(IVs.train, DV.train, IVs.test, type, balance, ncores, folds)`
- `NNS.caus(x, y)` (directional causation index for each direction, in
  ``` math
  0, 1
  ```
  )

### 6.2 Code: classification via regression + ensembles

``` r

# Example 1: Nonlinear regression
set.seed(123)
x_train <- runif(1000, -2, 2)
y_train <- sin(pi * x_train) + rnorm(1000, sd = 0.2)

x_test <- seq(-2, 2, length.out = 100)

NNS.reg(x = x_train, y = y_train, order = NULL, point.est = x_test)
```

![](NNSvignette_01_Overview_files/figure-html/unnamed-chunk-11-1.png)

    ## $R2
    ## [1] 0.9303023
    ## 
    ## $SE
    ## [1] 0.1971037
    ## 
    ## $Prediction.Accuracy
    ## NULL
    ## 
    ## $equation
    ## NULL
    ## 
    ## $x.star
    ## NULL
    ## 
    ## $derivative
    ##   Coefficient X.Lower.Range X.Upper.Range
    ## 1   2.5807500     -1.998139     -1.971809
    ## 2   3.9622780     -1.971809     -1.905718
    ## 3   3.8135391     -1.905718     -1.824605
    ## 4   0.6167479     -1.824605     -1.766196
    ## 5   2.1353953     -1.766196     -1.699572
    ## 6   0.8895008     -1.699572     -1.619703
    ## ---
    ## 60    2.001269      1.616540      1.669831
    ## 61    1.941175      1.669831      1.785845
    ## 62    2.868900      1.785845      1.846424
    ## 63    3.281862      1.846424      1.899405
    ## 64    2.491275      1.899405      1.963857
    ## 65    1.143921      1.963857      1.997618
    ## --- [ 65 rows x 3 cols ]; showing first and last 6. Use as.data.frame(x) for all. ---
    ## 
    ## $Point.est
    ##   [1]  0.017018681  0.138163442  0.298255482  0.454341905  0.608424293
    ##   [6]  0.677391810  0.715401475  0.801680074  0.859546521  0.895485946
    ##  [11]  0.941696907  0.995115180  0.988182153  0.969416723  0.977297922
    ##  [16]  0.954922821  0.913217684  0.892635178  0.770713048  0.646015031
    ##  [21]  0.541633189  0.457878622  0.323449082  0.150791607  0.077047336
    ##  [26]  0.003303065 -0.119433168 -0.268456262 -0.531786681 -0.518318419
    ##  [31] -0.671710663 -0.687771741 -0.768570199 -0.875828428 -0.988267823
    ##  [36] -0.975685627 -1.019584380 -1.028329777 -0.930342309 -0.915902057
    ##  [41] -0.873524550 -0.868204406 -0.836485641 -0.738944553 -0.641403466
    ##  [46] -0.525528161 -0.407420178 -0.264103654 -0.107021974  0.068133988
    ##  [51]  0.187065526  0.252315480  0.317565434  0.436013361  0.605309144
    ##  [56]  0.719889572  0.748101263  0.793617096  0.906453801  0.957472631
    ##  [61]  0.977868082  0.998263532  1.010496658  1.005245980  0.999995302
    ##  [66]  0.926157445  0.861276999  0.856777312  0.678223624  0.572049622
    ##  [71]  0.517111026  0.394714412  0.231889082  0.124984707  0.054956469
    ##  [76] -0.018741586 -0.172500238 -0.287018344 -0.393118160 -0.488942880
    ##  [81] -0.639500452 -0.802935287 -0.808798298 -0.856488778 -0.951329071
    ##  [86] -1.068485641 -0.997254747 -0.946701782 -0.951745343 -1.052665001
    ##  [91] -0.920193073 -0.839750588 -0.761319271 -0.682887954 -0.593199102
    ##  [96] -0.477283958 -0.348003866 -0.231046717 -0.130389154 -0.078428675
    ## 
    ## $pred.int
    ## NULL
    ## 
    ## $regression.points
    ##           x          y
    ## 1 -1.998139 0.02182248
    ## 2 -1.971809 0.08977377
    ## 3 -1.905718 0.35164124
    ## 4 -1.824605 0.66097084
    ## 5 -1.766196 0.69699463
    ## 6 -1.699572 0.83926179
    ## ---
    ##           x           y
    ## 61 1.669831 -0.85321621
    ## 62 1.785845 -0.62801195
    ## 63 1.846424 -0.45421675
    ## 64 1.899405 -0.28034153
    ## 65 1.963857 -0.11977316
    ## 66 1.997618 -0.08115337
    ## --- [ 66 rows x 2 cols ]; showing first and last 6. Use as.data.frame(x) for all. ---
    ## 
    ## $Fitted.xy
    ##            x          y      y.hat NNS.ID   gradient   residuals
    ## 1 -0.8496899 -0.5752368 -0.5475488     18 -0.7278556  0.02768805
    ## 2  1.1532205 -0.6617217 -0.3971628     50 -2.3716618  0.26455893
    ## 3 -0.3640923 -0.7048691 -0.8562335     30 -0.5794835 -0.15136436
    ## 4  1.5320696 -0.8447168 -0.9488136     57 -0.1248281 -0.10409685
    ## 5  1.7618691 -0.9820881 -0.6745538     61  1.9411751  0.30753433
    ## 6 -1.8177740  0.5226887  0.6651838      4  0.6167479  0.14249509
    ## ---
    ##               x          y      y.hat NNS.ID    gradient    residuals
    ## 995   0.2618615  0.5142993  0.7475673     37  0.69823937  0.233267974
    ## 996   1.3184955 -0.7988901 -0.8079711     53  0.15421526 -0.009081002
    ## 997   0.5684553  1.1554781  1.0022569     40 -0.12995428 -0.153221166
    ## 998  -0.4340050 -0.7748325 -0.9168735     28  0.09950728 -0.142040965
    ## 999   0.8383194  0.7041960  0.5034640     44 -1.35973026 -0.200732008
    ## 1000 -1.5647037  0.9467853  0.9830204      7  1.32210225  0.036235038
    ## --- [ 1000 rows x 6 cols ]; showing first and last 6. Use as.data.frame(x) for all. ---
    ## 
    ## $class.levels
    ## NULL

``` r

# Simple train/test for boosting & stacking
test.set = 141:150
 
boost <- NNS.boost(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              epochs = 10, learner.trials = 10, 
              status = FALSE, balance = TRUE,
              type = "CLASS")


mean(boost$results == as.numeric(iris[test.set,5]))
# [1] 1


boost$feature.weights; boost$feature.frequency

stacked <- NNS.stack(IVs.train = iris[-test.set, 1:4], 
                     DV.train = iris[-test.set, 5],
                     IVs.test = iris[test.set, 1:4],
                     type = "CLASS", balance = TRUE,
                     ncores = 1, folds = 2)
mean(stacked$stack == as.numeric(iris[test.set,5]))
# [1] 1
```

### 6.3 Code: directional causality

``` r

NNS.caus(mtcars$hp,  mtcars$mpg)  # hp -> mpg
```

    ## C(x--->y) C(y--->x) 
    ##   0.48695   0.51305

``` r

NNS.caus(mtcars$mpg, mtcars$hp)   # hp -> mpg
```

    ## C(x--->y) C(y--->x) 
    ##   0.51305   0.48695

**Interpretation.** Each call returns `C(x--->y)` and `C(y--->x)`, which
sum to 1; a value above 0.5 favours that direction, and swapping the
arguments swaps the two values. The method compares lagged pairs on the
copula scale (set `tau` for time series), so it is unaffected by the
units or scale of either variable.

------------------------------------------------------------------------

## 7. Time Series & Forecasting

**Headers.**

- `NNS.ARMA`
- `NNS.ARMA.optim`
- `NNS.seas`
- `NNS.VAR`

``` r

# Univariate nonlinear ARMA
set.seed(42)
z <- as.numeric(scale(sin(1:480/8) + rnorm(480, sd=.35)))

# Seasonality detection (prints a summary)
seasonal_period <- NNS.seas(z, plot = FALSE)
head(seasonal_period$all.periods)
```

    ##   Period Coefficient.of.Variation Variable.Coefficient.of.Variation
    ## 1     50                0.5090898                          1.150504
    ## 2     51                0.5532376                          1.150504
    ## 3     49                0.6885932                          1.150504
    ## 4     52                0.7495139                          1.150504
    ## 5     48                0.9370540                          1.150504
    ## 6     53                0.9374740                          1.150504

``` r

# Validate seasonal periods
NNS.ARMA.optim(z, h = 48, seasonal.factor = seasonal_period$periods, plot = TRUE, ncores = 1)
```

    ## [1] "CURRNET METHOD: lin"
    ## [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.425337970282028"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49, 51 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.314272381262221"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49, 51, 48 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.288043942135049"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49, 51, 48, 50 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.272024984385753"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49, 51, 48, 50, 47 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.270145105918314"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49, 51, 48, 50, 47, 53 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.267559815481421"
    ## [1] "NNS.ARMA(... method =  'lin' , seasonal.factor =  c( 49, 51, 48, 50, 47, 53, 25 ) ...)"
    ## [1] "CURRENT lin OBJECTIVE FUNCTION = 0.266123583980802"
    ## [1] "BEST method = 'lin', seasonal.factor = c( 49, 51, 48, 50, 47, 53, 25 )"
    ## [1] "BEST lin OBJECTIVE FUNCTION = 0.266123583980802"
    ## [1] "CURRNET METHOD: nonlin"
    ## [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
    ## [1] "NNS.ARMA(... method =  'nonlin' , seasonal.factor =  c( 49, 51, 48, 50, 47, 53, 25 ) ...)"
    ## [1] "CURRENT nonlin OBJECTIVE FUNCTION = 0.326722962745325"
    ## [1] "BEST method = 'nonlin' PATH MEMBER = c( 49, 51, 48, 50, 47, 53, 25 )"
    ## [1] "BEST nonlin OBJECTIVE FUNCTION = 0.326722962745325"
    ## [1] "CURRNET METHOD: both"
    ## [1] "COPY LATEST PARAMETERS DIRECTLY FOR NNS.ARMA() IF ERROR:"
    ## [1] "NNS.ARMA(... method =  'both' , seasonal.factor =  c( 49, 51, 48, 50, 47, 53, 25 ) ...)"
    ## [1] "CURRENT both OBJECTIVE FUNCTION = 0.275058604442101"
    ## [1] "BEST method = 'both' PATH MEMBER = c( 49, 51, 48, 50, 47, 53, 25 )"
    ## [1] "BEST both OBJECTIVE FUNCTION = 0.275058604442101"

![](NNSvignette_01_Overview_files/figure-html/unnamed-chunk-14-1.png)

    ## $periods
    ## [1] 49 51 48 50 47 53 25
    ## 
    ## $weights
    ## NULL
    ## 
    ## $obj.fn
    ## [1] 0.2661236
    ## 
    ## $method
    ## [1] "lin"
    ## 
    ## $shrink
    ## [1] FALSE
    ## 
    ## $nns.regress
    ## [1] FALSE
    ## 
    ## $bias.shift
    ## [1] 0
    ## 
    ## $bias.validation
    ## $bias.validation$policy
    ## [1] "none"
    ## 
    ## $bias.validation$estimated.shift
    ## [1] 0
    ## 
    ## $bias.validation$applied.shift
    ## [1] 0
    ## 
    ## $bias.validation$estimate.rows
    ## integer(0)
    ## 
    ## $bias.validation$acceptance.rows
    ## integer(0)
    ## 
    ## $bias.validation$raw.score
    ## [1] NA
    ## 
    ## $bias.validation$adjusted.score
    ## [1] NA
    ## 
    ## $bias.validation$accepted
    ## [1] FALSE
    ## 
    ## $bias.validation$status
    ## [1] "disabled"
    ## 
    ## 
    ## $raw.results
    ##  [1] -0.85581903 -0.96308748 -1.05600863 -1.24786301 -1.48681514 -1.30478168
    ##  [7] -1.35340984 -1.50006345 -1.41520981 -1.28993182 -1.22970244 -1.23421646
    ## [13] -1.17041069 -0.97209723 -0.88818032 -0.80061095 -0.65535580 -0.41468950
    ## [19] -0.40327882 -0.21161959 -0.31816376  0.01664552  0.43839280  0.30972663
    ## [25]  0.54213425  1.02400803  0.89891408  1.20583966  1.19311818  1.22864610
    ## [31]  1.27941996  1.35512214  1.08791716  1.34974950  1.33208959  1.30592055
    ## [37]  1.14524440  0.99653070  0.87842882  0.89786929  0.70893198  0.50263175
    ## [43]  0.26363731  0.21469814  0.06487646 -0.01096599 -0.16398874 -0.38813956
    ## 
    ## $validation.forecast
    ##  [1] -1.22924456 -1.52633220 -1.58439314 -1.64281360 -1.67245706 -1.40210102
    ##  [7] -1.59901938 -1.18889701 -1.39930879 -1.04702588 -1.00833507 -0.70279084
    ## [13] -0.82632734 -0.75501899 -0.32326480 -0.26948061 -0.26263844  0.29209784
    ## [19]  0.08072486  0.39790523  0.58602896  0.79118310  0.74897065  1.03701832
    ## [25]  0.80100746  1.10567655  1.00379912  0.99212370  1.21675681  1.37274012
    ## [31]  1.37655032  1.30884422  1.08905609  1.39257997  0.95762346  1.21645460
    ## [37]  1.16912616  0.89707781  0.81863549  0.35806587  0.18424184  0.13116453
    ## [43] -0.10026677 -0.20158338 -0.24994234 -0.26800697 -0.46057159 -1.01119665
    ## [49] -1.15545281 -0.99970402 -1.52865708 -1.56002748 -1.54290235 -1.61980183
    ## [55] -1.62866972 -1.68989695 -1.56650570 -1.38016720 -1.26930860 -1.22063576
    ## [61] -0.93891683 -0.99452769 -0.57757598 -0.71606699 -0.32339222 -0.42968508
    ## [67]  0.02819257  0.09906550  0.20705355  0.56192701  0.69847624  0.59829469
    ## [73]  0.78863121  1.07734609  1.07256100  1.03670644  1.00739471  1.19759711
    ## [79]  1.46737309  1.25270063  1.38305218  1.18766938  1.55729234  1.24891923
    ## [85]  1.41340818  0.95535328  1.26485338  0.93030190  0.67978857  0.16507589
    ## [91]  0.36466568 -0.04346508  0.17445007 -0.22449025 -0.23301587 -0.41706525
    ## 
    ## $errors
    ##  [1] -0.01998592 -0.22845283 -0.18710865 -0.61282561 -0.98294926 -0.39029217
    ##  [7] -0.43995297  0.21278349 -0.41773642  0.33224733  0.15511793 -0.27099786
    ## [13] -0.22161157 -0.19461416  0.31663663 -0.31932342 -0.65337189  0.73762111
    ## [19] -0.05344162  0.10754886  0.42513154  0.67317003  0.04228204 -0.09371272
    ## [25] -0.81416736  0.92486742  0.19627017 -0.05501108  0.73752763  0.45667485
    ## [31]  0.46877913  0.39260494  0.87242097  0.14464405  0.16617241 -0.45747776
    ## [37] -0.06977411  0.79736201  0.15054503  0.14713314 -0.54523292 -0.06283345
    ## [43]  0.14251743 -0.24775643 -0.00728412  0.27334440 -0.75917428 -0.01713549
    ## [49] -0.05988593  0.69704361 -0.41331683 -0.78277562  0.63949861 -0.21519330
    ## [55]  0.22400013 -1.01770498 -0.86993989 -0.25418868 -0.69170088  0.32441878
    ## [61] -0.30548791  0.22847339  0.11587307  0.12580313  0.36627610 -0.26417854
    ## [67]  0.87367425  0.30769749  0.06333232  0.51819261  0.26895872  0.21097657
    ## [73] -0.61577514 -0.04762451 -1.18775805  0.37617612 -0.46665845 -0.92483070
    ## [79]  0.10919559  0.49149282  0.57792556 -0.04574932 -0.25346112 -0.25855530
    ## [85]  1.28212825  0.26034097  0.36024003  0.79914461 -0.29461140  0.75998628
    ## [91] -0.06874751 -0.07870108  0.34459264  0.01324222 -0.10576509 -0.17013463
    ## 
    ## $results
    ##  [1] -0.85581903 -0.96308748 -1.05600863 -1.24786301 -1.48681514 -1.30478168
    ##  [7] -1.35340984 -1.50006345 -1.41520981 -1.28993182 -1.22970244 -1.23421646
    ## [13] -1.17041069 -0.97209723 -0.88818032 -0.80061095 -0.65535580 -0.41468950
    ## [19] -0.40327882 -0.21161959 -0.31816376  0.01664552  0.43839280  0.30972663
    ## [25]  0.54213425  1.02400803  0.89891408  1.20583966  1.19311818  1.22864610
    ## [31]  1.27941996  1.35512214  1.08791716  1.34974950  1.33208959  1.30592055
    ## [37]  1.14524440  0.99653070  0.87842882  0.89786929  0.70893198  0.50263175
    ## [43]  0.26363731  0.21469814  0.06487646 -0.01096599 -0.16398874 -0.38813956
    ## 
    ## $lower.pred.int
    ##  [1] -1.729023295 -1.836291749 -1.929212899 -2.121067283 -2.360019409
    ##  [6] -2.177985951 -2.226614114 -2.373267719 -2.288414075 -2.163136087
    ## [11] -2.102906706 -2.107420731 -2.043614959 -1.845301503 -1.761384585
    ## [16] -1.673815224 -1.528560074 -1.287893765 -1.276483091 -1.084823861
    ## [21] -1.191368026 -0.856558753 -0.434811472 -0.563477644 -0.331070018
    ## [26]  0.150803764  0.025709810  0.332635393  0.319913909  0.355441830
    ## [31]  0.406215689  0.481917874  0.214712887  0.476545234  0.458885322
    ## [36]  0.432716279  0.272040131  0.123326430  0.005224547  0.024665018
    ## [41] -0.164272292 -0.370572522 -0.609566963 -0.658506128 -0.808327813
    ## [46] -0.884170256 -1.037193011 -1.261343833
    ## 
    ## $upper.pred.int
    ##  [1]  0.105335772 -0.001932682 -0.094853832 -0.286708216 -0.525660342
    ##  [6] -0.343626884 -0.392255047 -0.538908652 -0.454055008 -0.328777020
    ## [11] -0.268547639 -0.273061664 -0.209255892 -0.010942436  0.072974482
    ## [16]  0.160543843  0.305798993  0.546465302  0.557875976  0.749535206
    ## [21]  0.642991040  0.977800314  1.399547595  1.270881423  1.503289049
    ## [26]  1.985162831  1.860068877  2.166994460  2.154272976  2.189800897
    ## [31]  2.240574756  2.316276941  2.049071954  2.310904301  2.293244389
    ## [36]  2.267075346  2.106399198  1.957685497  1.839583614  1.859024085
    ## [41]  1.670086775  1.463786545  1.224792104  1.175852939  1.026031254
    ## [46]  0.950188811  0.797166056  0.573015234

**Notes.** NNS seasonality uses coefficient of variation instead of
ACF/PACFs, and NNS ARMA blends multiple seasonal periods into the linear
or nonlinear regression forecasts.

------------------------------------------------------------------------

## 8. Simulation & Bootstrap & Risk‑Neutral Rescaling

### 8.1 Maximum entropy bootstrap (shape‑preserving)

**Header.**

- `NNS.meboot(x, reps=999, rho=NULL, type="spearman", drift=TRUE, ...)`

``` r

x_ts <- cumsum(rnorm(350, sd=.7))
mb <- NNS.meboot(x_ts, reps=5, rho = 1)
dim(mb["replicates", ]$replicates)
```

    ## [1] 350   5

### 8.2 Monte Carlo over the full correlation space

**Header.**

- `NNS.MC(x, reps=30, lower_rho=-1, upper_rho=1, by=.01, exp=1, type="spearman", ...)`

``` r

mc <- NNS.MC(x_ts, reps=5, lower_rho=-1, upper_rho=1, by=.5, exp=1)
length(mc$ensemble); names(mc$replicates)
```

    ## [1] 350

    ## [1] "rho = 1"    "rho = 0.5"  "rho = 0"    "rho = -0.5" "rho = -1"

``` r

head(mc$replicates$`rho = 0`)
```

    ##      Replicate 1 Replicate 2 Replicate 3 Replicate 4 Replicate 5
    ## [1,]    8.294844    6.022924    8.272604    5.594026    6.874941
    ## [2,]    4.725953   12.091955    9.629211    4.656695    4.514312
    ## [3,]   -3.662776    6.613393   11.215198    6.987755    3.488336
    ## [4,]   -2.017371    4.635251    5.903151    9.134349   12.628982
    ## [5,]   -7.888070    7.055852    9.099748    9.483581   13.727086
    ## [6,]   -3.961955    5.813247   10.437696    9.425552   13.587576

------------------------------------------------------------------------

## 9. Portfolio & Stochastic Dominance

Stochastic dominance orders uncertain prospects for broad classes of
risk‑averse utilities; partial moments supply practical, nonparametric
estimators.

**Headers.**

- `NNS.FSD.uni(x, y)`
- `NNS.SSD.uni(x, y)`
- `NNS.TSD.uni(x, y)`
- `NNS.SD.cluster(R)`
- `NNS.SD.efficient.set(R)`

``` r

set.seed(42)
RA <- rnorm(240, 0.005, 0.03)
RB <- rnorm(240, 0.003, 0.02)
RC <- rnorm(240, 0.006, 0.04)

NNS.FSD.uni(RA, RB)
```

    ## [1] 0

``` r

NNS.SSD.uni(RA, RB)
```

    ## [1] 0

``` r

NNS.TSD.uni(RA, RB)
```

    ## [1] 0

``` r

Rmat <- cbind(A=RA, B=RB, C=RC)
try(NNS.SD.cluster(Rmat, degree = 1))
```

    ## $Clusters
    ## $Clusters$Cluster_1
    ## [1] "A" "C" "B"

``` r

try(NNS.SD.efficient.set(Rmat, degree = 1))
```

    ## Checking 1 of 2Checking 2 of 2

    ## [1] "A" "C" "B"

------------------------------------------------------------------------

## Appendix A — Measure‑theoretic sketch (why partial moments are rigorous)

Let $`(\Omega, \mathcal{F}, \mathbb{P})`$ be a probability space,
$`X: \Omega\to\mathbb{R}`$ measurable. For any fixed $`t\in\mathbb{R}`$,
the sets $`\{X\le t\}`$ and $`\{X>t\}`$ are in $`\mathcal{F}`$ because
they are preimages of Borel sets. The **population** partial moments are

``` math
\operatorname{LPM}(k,t,X) = \int_{-\infty}^{t} (t-x)^k\, dF_X(x),
\qquad
\operatorname{UPM}(k,t,X) = \int_{t}^{\infty} (x-t)^k\, dF_X(x).
```

The **empirical** versions correspond to replacing $`F_X`$ with the
empirical measure $`\mathbb{P}_n`$ (or CDF $`\hat F_n`$):

``` math
\widehat{\operatorname{LPM}}_k(t;X) = \int_{(-\infty,t]} (t-x)^k\, d\mathbb{P}_n(x),
\qquad
\widehat{\operatorname{UPM}}_k(t;X) = \int_{(t,\infty)} (x-t)^k\, d\mathbb{P}_n(x).
```

Centering at $`t=\mu_X`$ yields the variance decomposition identity in
Section 1.

------------------------------------------------------------------------

## Appendix B — Quick Reference (Grouped by Topic)

### Overall Theory

- [Nonlinear Nonparametric Statistics: Using Partial
  Moments](https://ovvo-financial.github.io/NNS/book/)

### 1. Partial Moments & Ratios

- `LPM(degree, target, variable)` — lower partial moment of order
  `degree` at `target`.
- `UPM(degree, target, variable)` — upper partial moment of order
  `degree` at `target`.
- `LPM.ratio(degree, target, variable)`; `UPM.ratio(...)` — normalized
  shares; `degree=0` gives CDF.
- `LPM.VaR(p, degree, variable)` — partial-moment quantile at
  probability `p`.
- `Co.LPM(degree_lpm, x, y, target_x, target_y, degree_y)` — co-lower
  partial moment between two variables.
- `Co.UPM(degree_upm, x, y, target_x, target_y, degree_y)` — co-upper
  partial moment between two variables.
- `D.LPM(degree, target, variable)` — divergent lower partial moment
  (away from `target`).
- `D.UPM(degree, target, variable)` — divergent upper partial moment
  (away from `target`).
- `NNS.CDF(x, target = NULL, points = NULL, plot = TRUE/FALSE)` — CDF
  from partial moments.
- `NNS.moments(x)` — mean/var/skew/kurtosis via partial moments.

### 2. Descriptive Statistics & Distributions

- `NNS.mode(x, multi = FALSE)` — nonparametric mode(s).
- `PM.matrix(l_degree, u_degree, target, variable, pop_adj)` —
  co-/divergent partial-moment matrices.
- `NNS.gravity(x, w = NULL)` — partial-moment weighted location (gravity
  center).

See NNS Vignette: [Getting Started with NNS: Partial
Moments](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_02_Partial_Moments.md)

### 3. Dependence & Association

- `NNS.dep(x, y)` — nonlinear dependence coefficient.
- `NNS.copula(X, target, continuous, plot, independence.overlay)` —
  dependence from co-partial moments.

See NNS Vignette: [Getting Started with NNS: Correlation and
Dependence](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_03_Correlation_and_Dependence.md)

### 4. Normalization & Rescaling

- `NNS.norm(x, linear=FALSE)` — normalization retaining target moments.
- `NNS.rescale(x, a, b, method=c("minmax","riskneutral"), T=NULL, type=c("Terminal","Discounted"))`
  — risk-neutral or min–max rescaling.

See NNS Vignette: [Getting Started with NNS: Normalization and
Rescaling](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_04_Normalization_and_Rescaling.md)

### 5. Hypothesis Testing

- `NNS.ANOVA(control, treatment, ...)` — certainty of equality
  (distributions or means).
- `NNS.SS(x, y, ...)` — stochastic superiority between two variables.

See NNS Vignette: [Getting Started with NNS: Comparing
Distributions](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_06_Comparing_Distributions.md)

### 6. Regression, Classification & Causality

- `NNS.part(x, y, ...)` — partition analysis for variable segmentation.
- `NNS.reg(x, y, ...)` — partition-based regression/classification
  (`$Fitted.xy`, `$Point.est`).
- `NNS.boost(IVs, DV, ...)`, `NNS.stack(IVs, DV, ...)` — ensembles using
  `NNS.reg` base learners.
- `NNS.caus(x, y)` — directional causation index.

See NNS Vignette: [Getting Started with NNS: Clustering and
Regression](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_07_Clustering_and_Regression.md)

See NNS Vignette: [Getting Started with NNS:
Classification](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_08_Classification.md)

### 7. Differentiation & Slope Measures

- `dy.dx(x, y)` — numerical derivative of `y` with respect to `x` via
  `NNS.reg`.
- `dy.d_(x, Y, var)` — partial derivative of multivariate `Y` w.r.t.
  `var`.
- `NNS.diff(x, y)` — derivative via secant projections.

### 8. Time Series & Forecasting

- `NNS.ARMA(...)`, `NNS.ARMA.optim(...)` — nonlinear ARMA modeling.
- `NNS.seas(...)` — detect seasonality.
- `NNS.VAR(...)` — nonlinear VAR modeling.
- `NNS.nowcast(x, h, ...)` — near-term nonlinear forecast.

See NNS Vignette: [Getting Started with NNS:
Forecasting](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_09_Forecasting.md)

### 9. Simulation & Bootstrap

- `NNS.meboot(...)` — maximum entropy bootstrap.
- `NNS.MC(...)` — Monte Carlo over correlation space.

See NNS Vignette: [Getting Started with NNS: Sampling and
Simulation](https://OVVO-Financial.github.io/NNS/articles/NNSvignette_05_Sampling.md)

### 10. Portfolio Analysis & Stochastic Dominance

- `NNS.FSD.uni(x, y)`, `NNS.SSD.uni(x, y)`, `NNS.TSD.uni(x, y)` —
  univariate stochastic dominance tests.
- `NNS.SD.cluster(R)`, `NNS.SD.efficient.set(R)` — dominance-based
  portfolio sets.

For complete references, please see the Vignettes linked above and their
specific referenced materials.
