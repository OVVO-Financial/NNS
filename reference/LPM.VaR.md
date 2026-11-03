# LPM VaR

Generates a value at risk (VaR) quantile based on the Lower Partial
Moment ratio.

## Usage

``` r
LPM.VaR(percentile, degree, x)
```

## Arguments

- percentile:

  numeric \[0, 1\]; The percentile for left-tail VaR.

- degree:

  integer; `(degree = 0)` gives the interpolated empirical quantile
  (`quantile(x, percentile, type = 7)`), which on tied data need not be
  an exact inverse of `LPM.ratio(0, ., x)`. Use it for probability
  statements, such as a VaR exceeded with a given frequency.
  `(degree = 1)` gives an area-weighted (severity-weighted) threshold,
  the `percentile`-expectile of `x` (Newey and Powell 1987):
  `percentile` of the area of deviations, not of the observations, lies
  below it, so it is not a quantile of `x` (for exponential data the
  degree-1 5% point has about a quarter of the observations below it).
  Higher degrees weight large deviations more.

- x:

  a numeric vector.

## Value

Returns the type-7 empirical quantile (`degree = 0`), whose strict or
inclusive sample tail fraction may differ from `percentile` because of
interpolation or ties, or the threshold with `percentile` of the
directional deviation area below it (`degree > 0`).

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)

## For 5th percentile, left-tail
LPM.VaR(0.05, 0, x)
} # }
```
