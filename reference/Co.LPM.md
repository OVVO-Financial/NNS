# Co‑Lower Partial Moment

Computes the co‑lower partial moment (lower‑left quadrant 4) between two
equal‑length numeric vectors at any degree and target.

## Usage

``` r
Co.LPM(degree_lpm, x, y, target_x, target_y, degree_y = NULL)
```

## Arguments

- degree_lpm:

  numeric; degree for x ("degree_x"). degree = 0 gives frequency, degree
  = 1 gives area.

- x:

  numeric vector of observations.

- y:

  numeric vector of the same length as x.

- target_x:

  numeric vector; thresholds for x. There is no default: pass `NULL` to
  use `mean(x)`.

- target_y:

  numeric vector; thresholds for y. There is no default: pass `NULL` to
  use `mean(y)`. `target_x` and `target_y` must be non-empty and the
  same length, or one of them length 1.

- degree_y:

  numeric; optional degree for y. If omitted, \`degree_lpm\` is used for
  both x and y.

## Value

Numeric vector of co‑LPM values.

## References

Viole, F. & Nawrocki, D. (2013) \*Nonlinear Nonparametric Statistics:
Using Partial Moments\* (ISBN:1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
  set.seed(123)
  x <- rnorm(100); y <- rnorm(100)
  Co.LPM(0, x, y, mean(x), mean(y))
#> [1] 0.31
```
