# NNS Numerical Differentiation

Determines numerical derivative of a given univariate function using
projected secant lines on the y-axis. These projected points infer
finite steps `h`, in the finite step method.

## Usage

``` r
NNS.diff(
  f,
  point,
  h = abs(point) * 0.1 + 0.01,
  tol = 1e-10,
  max.iter = NULL,
  digits = 12,
  print.trace = FALSE,
  plot = FALSE
)
```

## Arguments

- f:

  an expression or call or a formula with no lhs.

- point:

  numeric; Point to be evaluated for derivative of a given function `f`.

- h:

  numeric \[0, ...\]; Initial step for secant projection. Defaults to
  `(h = abs(point) * 0.1 + 0.01)`.

- tol:

  numeric; Sets the tolerance for the stopping condition of the inferred
  `h`. Defaults to `(tol = 1e-10)`.

- max.iter:

  integer; `NULL` (default) Maximum number of bisection iterations.
  `NULL` sets the limit to `100L`. For noisy functions the bisection may
  stall before `tol` is reached; `max.iter` provides a hard upper bound.

- digits:

  numeric; Sets the number of digits specification of the output.
  Defaults to `(digits = 12)`.

- print.trace:

  logical; `FALSE` (default) Displays each iteration, lower y-intercept,
  upper y-intercept and inferred `h`.

- plot:

  logical; plots range, secant lines and y-intercept convergence.

## Value

Returns a matrix of values, intercepts, derivatives, inferred step sizes
for multiple methods of estimation.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
f <- function(x) sin(x) / x
NNS.diff(f, 4.1)

## Noisy function with explicit iteration cap
f_noisy <- function(x) sin(x) + rnorm(1, 0, 0.001)
NNS.diff(f_noisy, 1.0, max.iter = 100)
} # }
```
