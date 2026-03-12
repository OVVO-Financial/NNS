# NNS Co-Partial Moments Higher Dimension Dependence

Determines higher dimension dependence coefficients based on co-partial
moment matrices ratios.

## Usage

``` r
NNS.copula(
  X,
  target = NULL,
  continuous = TRUE,
  plot = FALSE,
  independence.overlay = FALSE
)
```

## Arguments

- X:

  a numeric matrix or data frame.

- target:

  numeric; Typically the mean of Variable X for classical statistics
  equivalences, but does not have to be. (Vectorized) `(target = NULL)`
  (default) will set the target as the mean of every variable.

- continuous:

  logical; `TRUE` (default) Generates a continuous measure using degree
  1
  [PM.matrix](https://OVVO-Financial.github.io/NNS/reference/PM.matrix.md),
  while discrete `FALSE` uses degree 0
  [PM.matrix](https://OVVO-Financial.github.io/NNS/reference/PM.matrix.md).

- plot:

  logical; `FALSE` (default) Generates a 3d scatter plot with regression
  points.

- independence.overlay:

  logical; `FALSE` (default) Creates and overlays independent
  [Co.LPM](https://OVVO-Financial.github.io/NNS/reference/Co.LPM.md) and
  [Co.UPM](https://OVVO-Financial.github.io/NNS/reference/Co.UPM.md)
  regions to visually reference the difference in dependence from the
  data.frame of variables being analyzed. Under independence, the light
  green and red shaded areas would be occupied by green and red data
  points respectively.

## Value

Returns a multivariate dependence value \[0,1\].

## References

Viole, F. (2016) "Beyond Correlation: Using the Elements of Variance for
Conditional Means and Probabilities"
[doi:10.2139/ssrn.2745308](https://doi.org/10.2139/ssrn.2745308) .

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(1000) ; y <- rnorm(1000) ; z <- rnorm(1000)
A <- data.frame(x, y, z)
NNS.copula(A, target = colMeans(A), plot = TRUE, independence.overlay = TRUE)

### Target 0
NNS.copula(A, target = rep(0, ncol(A)), plot = TRUE, independence.overlay = TRUE)
} # }
```
