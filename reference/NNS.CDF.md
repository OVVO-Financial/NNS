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

  numeric; `NULL` (default) Must lie within support of each variable.

- type:

  options("CDF", "survival", "hazard", "cumulative hazard"); `"CDF"`
  (default) Selects type of function to return for bi-variate analysis.
  Multivariate analysis is restricted to `"CDF"`.

- plot:

  logical; plots CDF.

## Value

Returns:

- `"Function"` a data.table containing the observations and resulting
  CDF of the variable.

- `"target.value"` value from the `target` argument.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

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
