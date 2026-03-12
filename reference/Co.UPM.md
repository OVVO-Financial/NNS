# Co‑Upper Partial Moment

Computes the co‑upper partial moment (upper‑right quadrant 1) between
two equal‑length numeric vectors at any degree and target.

## Usage

``` r
Co.UPM(degree_upm, x, y, target_x, target_y, degree_y = NULL)
```

## Arguments

- degree_upm:

  numeric; degree for x ("degree_x"). degree = 0 gives frequency, degree
  = 1 gives area.

- x:

  numeric vector of observations.

- y:

  numeric vector of the same length as x.

- target_x:

  numeric vector; thresholds for x (defaults to mean(x)).

- target_y:

  numeric vector; thresholds for y (defaults to mean(y)).

- degree_y:

  numeric; optional degree for y. If omitted, \`degree_upm\` is used for
  both x and y.

## Value

Numeric vector of co‑UPM values.

## References

Viole, F. & Nawrocki, D. (2013) \*Nonlinear Nonparametric Statistics:
Using Partial Moments\* (ISBN:1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
  set.seed(123)
  x <- rnorm(100); y <- rnorm(100)
  Co.UPM(0, x, y, mean(x), mean(y))
#> [1] 0.27
```
