# Partial Derivative dy/d\_\[wrt\]

Returns the numerical partial derivative of `y` with respect to \[wrt\]
any regressor for a point of interest. Finite difference method is used
with
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
estimates as `f(x + h)` and `f(x - h)` values.

## Usage

``` r
dy.d_(x, y, wrt, eval.points = "obs", mixed = FALSE, messages = TRUE)
```

## Arguments

- x:

  a numeric matrix or data frame.

- y:

  a numeric vector with compatible dimensions to `x`.

- wrt:

  integer; Selects the regressor to differentiate with respect to
  (vectorized).

- eval.points:

  numeric or options: ("obs", "apd", "mean", "median", "last");
  Regressor points to be evaluated.

  - Numeric values must be in matrix or data.frame form to be evaluated
    for each regressor, otherwise, a vector of points will evaluate only
    at the `wrt` regressor. See examples for use cases.

  - Set to `(eval.points = "obs")` (default) to find the average partial
    derivative at every observation of the variable with respect to *for
    specific tuples of given observations.*

  - Set to `(eval.points = "apd")` to find the average partial
    derivative at every observation of the variable with respect to
    *over the entire distribution of other regressors.*

  - Set to `(eval.points = "mean")` to find the partial derivative at
    the mean of value of every variable.

  - Set to `(eval.points = "median")` to find the partial derivative at
    the median value of every variable.

  - Set to `(eval.points = "last")` to find the partial derivative at
    the last observation of every value (relevant for time-series data).

- mixed:

  logical; `FALSE` (default) If mixed derivative is to be evaluated, set
  `(mixed = TRUE)`.

- messages:

  logical; `TRUE` (default) Prints status messages.

## Value

Returns column-wise matrix of wrt regressors:

- `dy.d_(...)[, wrt]$First` the 1st derivative

- `dy.d_(...)[, wrt]$Second` the 2nd derivative

- `dy.d_(...)[, wrt]$Mixed` the mixed derivative (for two independent
  variables only).

## Note

For binary regressors, it is suggested to use
`eval.points = seq(0, 1, .05)` for a better resolution around the
midpoint.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

Vinod, H. and Viole, F. (2020) "Comparing Old and New Partial Derivative
Estimates from Nonlinear Nonparametric Regressions"
[doi:10.2139/ssrn.3681104](https://doi.org/10.2139/ssrn.3681104)

## Author

Fred Viole, OVVO Financial Systems

## Examples
