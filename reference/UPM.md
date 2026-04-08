# Upper Partial Moment

This function generates a univariate upper partial moment for any degree
or target.

## Usage

``` r
UPM(degree, target, variable, excess_ret = FALSE)
```

## Arguments

- degree:

  numeric; `(degree = 0)` is frequency, `(degree = 1)` is area.

- target:

  numeric; Set to `target = mean(variable)` for classical equivalences,
  but does not have to be. When `excess_ret = FALSE`, this can be a
  scalar or a vectorized target for the standard partial moment
  calculation. When `excess_ret = TRUE`, it is interpreted element-wise
  as the benchmark/threshold relative to `variable`.

- variable:

  a numeric vector. [data.frame](https://rdrr.io/r/base/data.frame.html)
  or [list](https://rdrr.io/r/base/list.html) type objects are not
  permissible.

- excess_ret:

  logical; `FALSE` (default). If `TRUE`, switches from the standard
  vectorized-target partial moment to an element-wise excess-deviation
  calculation. For `UPM`, this computes `pmax(variable - target, 0)`
  raised to `degree` and averaged. In this mode, `target` must have
  length 1 or the same length as `variable`.

## Value

UPM of variable

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
UPM(0, mean(x), x)
#> [1] 0.49
```
