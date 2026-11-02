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

  numeric; `(degree = 0)` (default) is frequency, `(degree = 1)` is
  area.

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
  at `degree > 0`, normalised like the CDF. The hazard and cumulative
  hazard are built from it.

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

## Continuous CDF (degree = 1)
NNS.CDF(x, 1)

## Joint CDF
x <- rnorm(5000) ; y <- rnorm(5000)
A <- cbind(x,y)

NNS.CDF(A, 0)

## Joint CDF with target
NNS.CDF(A, 0, target = rep(0, ncol(A)))
} # }
```
