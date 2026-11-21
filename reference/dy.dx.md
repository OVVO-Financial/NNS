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

  Finite numeric points, or `"overall"` for the average first
  derivative.

## Value

For numeric points, a `data.frame` with the first and second derivatives
in the supplied order, including repeated points. For `"overall"`, a
single first derivative.

## Details

`eval.point = "overall"` returns the mean of the `gradient` column of
`NNS.reg(x, y)$Fitted.xy`. Each observation receives the slope of the
unsmoothed regression segment that contains it, so the result is an
observation-weighted average of segment slopes, not an average of
smoothed derivatives. Very short segments between nearly coincident
regression points can have extreme slopes;
`mean(dy.dx(x, y, eval.point = x)$first.derivative)` averages the
smoothed pointwise estimates instead.

Numeric evaluation points use central differences of smoothed
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
estimates over five steps from 2 to `floor(sqrt(n))` typical spacings of
`x`. The left point is bounded below by `min(x)` and the right point
above by `max(x)`. At a support endpoint with at least three distinct
observations of `x`, the second difference uses the midpoint of the
one-sided interval. At coordinates outside that range, the differences
use the extrapolated fit. `second.derivative` is a second difference of
those estimates, and second differences amplify noise: for `y = x^2`
plus noise with an R-squared near 0.99 and n = 400, estimates of the
true value 2 at three interior points had a root mean squared error of
about 5, and about one in five had the wrong sign. Treat it as a
descriptive diagnostic on noisy data.

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
dy.dx(x, y, eval.point = 1.75)$first.derivative

# Second derivative
dy.dx(x, y, eval.point = 1.75)$second.derivative

# Vector of derivatives
dy.dx(x, y, eval.point = c(1.75, 2.5))
} # }
```
