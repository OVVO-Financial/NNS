# NNS CDF

This function generates an empirical CDF using partial moment ratios
[LPM.ratio](https://OVVO-Financial.github.io/NNS/reference/LPM.ratio.md),
and resulting survival, hazard and cumulative hazard functions.

## Usage

``` r
NNS.CDF(variable, degree = 0, target = NULL, type = "CDF", plot = TRUE)
```

## Arguments

- variable:

  a numeric vector or data.frame of \>= 2 variables for joint CDF.

- degree:

  numeric; `(degree = 0)` (default) is frequency, the empirical
  distribution of `variable`; `(degree = 1)` is area, the share of
  deviation area below each point, which is not \\P(X \le t)\\.

- target:

  numeric; `NULL` (default) Must lie within support of each variable
  (one value per column for multivariate `variable`).

- type:

  options("CDF", "survival", "hazard", "cumulative hazard"); `"CDF"`
  (default) Selects type of function to return. For multivariate
  `variable` the joint survival is \\P(X_1 \> t_1, \ldots, X_d \> t_d)\\
  at `degree = 0`, by inclusion-exclusion over the column subsets for
  \\d \le 3\\ and as the share of observations above \\t\\ in every
  column for larger \\d\\, and the joint upper partial moment
  [Co.UPM_nD](https://OVVO-Financial.github.io/NNS/reference/Co.UPM_nD.md)
  at `degree > 0`, normalised like the CDF. The multivariate cumulative
  hazard is the transform \\-\log S(t)\\ and requires `degree = 0`; it
  is not automatically an integral of a scalar multivariate hazard.
  Multivariate `type = "hazard"` is unsupported and raises an error: a
  joint density-over-survival or a coordinate-hazard definition requires
  a separate estimator, and dividing the joint CDF by survival does not
  estimate either. At `degree > 0` with two or more columns, each
  target's co-partial moment is normalised by its own total, so the
  "CDF" and "survival" values are generalised co-partial-moment ratios,
  not a joint distribution function: they can assign a negative
  probability to a rectangle. At `degree = 0`, and in one dimension at
  any degree, they are proper distribution functions; at positive degree
  the univariate induced law describes deviation-area balance rather
  than the frequency law of the original variable.

- plot:

  logical; plots CDF.

## Value

Returns:

- `"Function"` a data.frame (class `NNS.data.frame`) containing the
  observations and the resulting function of the variable.

- `"target.value"` value from the `target` argument.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"
[doi:10.2139/ssrn.3007373](https://doi.org/10.2139/ssrn.3007373)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)
NNS.CDF(x)

## Empirical CDF (degree = 0)
NNS.CDF(x)

## Area-weighted CDF (degree = 1): the share of deviation area below each point
NNS.CDF(x, 1)

## Joint CDF
x <- rnorm(5000) ; y <- rnorm(5000)
A <- cbind(x,y)

NNS.CDF(A, 0)

## Joint CDF with target
NNS.CDF(A, 0, target = rep(0, ncol(A)))
} # }
```
