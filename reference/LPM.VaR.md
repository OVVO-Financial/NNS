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

  integer; `(degree = 0)` for discrete distributions, `(degree = 1)` for
  continuous distributions. With `degree = 0` the result is the
  interpolated empirical quantile (`quantile(x, percentile, type = 7)`),
  which on tied data need not be an exact inverse of
  `LPM.ratio(0, ., x)`.

- x:

  a numeric vector.

## Value

Returns a numeric value representing the point at which `"percentile"`
of the area of `x` is below.

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
