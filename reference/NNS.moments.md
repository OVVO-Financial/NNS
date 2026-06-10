# NNS moments

This function returns the first 4 moments of the distribution.

## Usage

``` r
NNS.moments(x, population = TRUE)
```

## Arguments

- x:

  a numeric vector.

- population:

  logical; `TRUE` (default) Performs the population adjustment.
  Otherwise returns the sample statistic.

## Value

Returns:

- `"$mean"` mean of the distribution.

- `"$variance"` variance of the distribution.

- `"$skewness"` skewness of the distribution.

- `"$kurtosis"` excess kurtosis.

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
NNS.moments(x)
} # }
```
