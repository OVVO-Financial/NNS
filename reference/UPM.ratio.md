# Upper Partial Moment Ratio

This function generates a standardized univariate upper partial moment
of any non‑negative degree for a given target.

## Usage

``` r
UPM.ratio(degree, target, variable)
```

## Arguments

- degree:

  numeric; degree = 0 gives frequency, degree = 1 gives area.

- target:

  numeric vector; threshold(s). There is no default: pass `NULL` to use
  `mean(variable)`.

- variable:

  numeric vector or data‑frame column to evaluate.

## Value

Numeric vector of standardized upper partial moments. For `degree > 0`,
a target at which every observation equals the target (0/0) returns the
`degree = 0` ratio.

## References

Viole, F. & Nawrocki, D. (2013) \*Nonlinear Nonparametric Statistics:
Using Partial Moments\* (ISBN:1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
  set.seed(123)
  x <- rnorm(100)
  UPM.ratio(0, mean(x), x)
#> [1] 0.49
if (FALSE) { # \dontrun{
  plot(sort(x), UPM.ratio(0, sort(x), x))
} # }
```
