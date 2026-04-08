# Partial Derivative dy/dx

Returns the numerical partial derivative of `y` wrt `x` for a point of
interest.

## Usage

``` r
dy.dx(x, y, eval.point = NULL)
```

## Arguments

- x:

  a numeric vector.

- y:

  a numeric vector.

- eval.point:

  numeric or ("overall"); `x` point to be evaluated, must be provided.
  Defaults to `(eval.point = NULL)`. Set to `(eval.point = "overall")`
  to find an overall partial derivative estimate (1st derivative only).

## Value

Returns a `data.table` of eval.point along with both 1st and 2nd
derivative.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
[doi:10.1007/s10614-017-9713-5](https://doi.org/10.1007/s10614-017-9713-5)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
x <- seq(0, 2 * pi, pi / 100) ; y <- sin(x)
dy.dx(x, y, eval.point = 1.75)

# First derivative
dy.dx(x, y, eval.point = 1.75)[ , first.derivative]

# Second derivative
dy.dx(x, y, eval.point = 1.75)[ , second.derivative]

# Vector of derivatives
dy.dx(x, y, eval.point = c(1.75, 2.5))
} # }
```
