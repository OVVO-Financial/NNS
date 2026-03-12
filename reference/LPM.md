# Lower Partial Moment

This function generates a univariate lower partial moment for any degree
or target.

## Usage

``` r
LPM(degree, target, variable, excess_ret = FALSE)
```

## Arguments

- degree:

  numeric; `(degree = 0)` is frequency, `(degree = 1)` is area.

- target:

  numeric; Set to `target = mean(variable)` for classical equivalences,
  but does not have to be. (Vectorized)

- variable:

  a numeric vector. [data.frame](https://rdrr.io/r/base/data.frame.html)
  or [list](https://rdrr.io/r/base/list.html) type objects are not
  permissible.

- excess_ret:

  logical; `FALSE` (default)

## Value

LPM of variable

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
set.seed(123)
x <- rnorm(100)
LPM(0, mean(x), x)
#> [1] 0.51
```
