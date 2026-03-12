# Divergent‑Lower Partial Moment

Computes the divergent lower partial moment (lower‑right quadrant 3)
between two equal‑length numeric vectors.

## Usage

``` r
D.LPM(degree_lpm, degree_upm, x, y, target_x, target_y)
```

## Arguments

- degree_lpm:

  numeric; LPM degree = 0 gives frequency, = 1 gives area.

- degree_upm:

  numeric; UPM degree = 0 gives frequency, = 1 gives area.

- x:

  numeric vector of observations.

- y:

  numeric vector of the same length as x.

- target_x:

  numeric vector; thresholds for x (defaults to mean(x)).

- target_y:

  numeric vector; thresholds for y (defaults to mean(y)).

## Value

Numeric vector of divergent LPM values.

## References

Viole, F. & Nawrocki, D. (2013) \*Nonlinear Nonparametric Statistics:
Using Partial Moments\* (ISBN:1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
  set.seed(123)
  x <- rnorm(100); y <- rnorm(100)
  D.LPM(0, 0, x, y, mean(x), mean(y))
#> [1] 0.22
```
