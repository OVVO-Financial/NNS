# UPM VaR

Generates an upside value at risk (VaR) quantile based on the Upper
Partial Moment ratio.

## Usage

``` r
UPM.VaR(percentile, degree, x)
```

## Arguments

- percentile:

  numeric \[0, 1\]; The percentile for right-tail VaR.

- degree:

  integer; `(degree = 0)` gives the interpolated empirical quantile
  (`quantile(x, 1 - percentile, type = 7)`), which on tied data need not
  be an exact inverse of `UPM.ratio(0, ., x)`. Use it for probability
  statements. `(degree = 1)` gives an area-weighted (severity-weighted)
  threshold, the (1 - `percentile`)-expectile of `x` (Newey and Powell
  1987): `percentile` of the area of deviations, not of the
  observations, lies above it, so it is not a quantile of `x`. Higher
  degrees weight large deviations more.

- x:

  a numeric vector.

## Value

Returns the type-7 empirical quantile (`degree = 0`), whose strict or
inclusive sample tail fraction may differ from `percentile` because of
interpolation or ties, or the threshold with `percentile` of the
directional deviation area above it (`degree > 0`).

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
set.seed(123)
x <- rnorm(100)

## For 5th percentile, right-tail
UPM.VaR(0.05, 0, x)
#>      95% 
#> 1.566526 
```
