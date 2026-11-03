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
    at the `wrt` regressor, averaging over the observed rows of the
    other regressors (at most `max(21, 10000 / k)` of them for `k`
    points, taken evenly through the data). See examples for use cases.

  - Set to `(eval.points = "obs")` (default) to find the average partial
    derivative at every observation of the variable with respect to *for
    specific tuples of given observations.*

  - Set to `(eval.points = "apd")` to find the average partial
    derivative at every observation of the variable with respect to
    *over the entire distribution of other regressors*, represented by
    their observed rows.

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

## Details

Derivatives are finite-difference summaries of the fitted
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
surface, averaged over four local probability windows derived from
degree-1 partial-moment quantiles. Window sizes depend on the number of
observations; they are capped before constructing a local coordinate
half-width. All evaluation points and perturbed coordinates must stay
inside the observed coordinate ranges. Interior unequal spacings use
polynomial finite-difference weights; support endpoints use one-sided
three-point formulas. This avoids averaging artificial full-support
extrapolation bands into a local derivative.

The coordinate-range check does not establish that every perturbed tuple
lies in a densely observed part of the joint support. Sparse cells,
neighbour changes, fitted-surface bias, and correlated predictors can
affect the result, including at the boundary. These are numerical
derivatives of an estimated surface, with no universal recovery or
coverage guarantee. A two-valued coordinate returns its finite contrast
divided by the coordinate difference and an undefined (NA) second
derivative. Constant coordinates return NA. Mixed derivatives require
two predictors and explicit matrix evaluation points.

## Note

Binary regressors provide a finite contrast, not an identified
continuous derivative between the two levels.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Vinod, H. and Viole, F. (2020) "Comparing Old and New Partial Derivative
Estimates from Nonlinear Nonparametric Regressions"
[doi:10.2139/ssrn.3681104](https://doi.org/10.2139/ssrn.3681104)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123) ; x_1 <- runif(1000) ; x_2 <- runif(1000) ; y <- x_1 ^ 2 * x_2 ^ 2
B <- cbind(x_1, x_2)

## To find derivatives of y wrt 1st regressor for specific points of both regressors
dy.d_(B, y, wrt = 1, eval.points = t(c(.5, median(x_2))))

## To find average partial derivative of y wrt 1st regressor,
## only supply 1 value in [eval.points], or a vector of [eval.points]:
dy.d_(B, y, wrt = 1, eval.points = .5)

dy.d_(B, y, wrt = 1, eval.points = fivenum(B[,1]))


## To find average partial derivative of y wrt 1st regressor,
## for every observation of 1st regressor:
apd <- dy.d_(B, y, wrt = 1, eval.points = "apd")
plot(B[,1], apd[,1]$First)

## Central 95% range across estimated partial derivatives (descriptive)
### Lower CI
LPM.VaR(.025, 0, apd[,1]$First)

### Upper CI
UPM.VaR(.025, 0, apd[,1]$First)
} # }
```
