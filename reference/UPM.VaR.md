# UPM VaR

Generates an upside value at risk (VaR) quantile based on the Upper
Partial Moment ratio

## Usage

``` r
UPM.VaR(percentile, degree, x)
```

## Arguments

- percentile:

  numeric \[0, 1\]; The percentile for right-tail VaR (vectorized).

- degree:

  integer; `(degree = 0)` for discrete distributions, `(degree = 1)` for
  continuous distributions.

- x:

  a numeric vector.

## Value

Returns a numeric value representing the point at which `"percentile"`
of the area of `x` is above.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

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
