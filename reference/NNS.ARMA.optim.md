# NNS ARMA Optimizer

Wrapper function for optimizing any combination of a given
`seasonal.factor` vector in
[NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md).
Minimum sum of squared errors (forecast-actual) is used to determine
optimum across all
[NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
methods.

## Usage

``` r
NNS.ARMA.optim(
  variable,
  h = NULL,
  training.set = NULL,
  seasonal.factor,
  lin.only = FALSE,
  negative.values = FALSE,
  obj.fn = expression(mean((predicted - actual)^2)/(NNS::Co.LPM(1, predicted, actual,
    target_x = mean(predicted), target_y = mean(actual)) + NNS::Co.UPM(1, predicted,
    actual, target_x = mean(predicted), target_y = mean(actual)))),
  objective = "min",
  linear.approximation = TRUE,
  ncores = NULL,
  pred.int = 0.95,
  print.trace = TRUE,
  plot = FALSE
)
```

## Arguments

- variable:

  a numeric vector.

- h:

  integer; `NULL` (default) Number of periods to forecast out of sample.
  If `NULL`, `h = length(variable) - training.set`.

- training.set:

  integer; `NULL` (default) Sets the number of variable observations as
  the training set. See `Note` below for recommended uses.

- seasonal.factor:

  integers; Multiple frequency integers considered for
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  model, i.e. `(seasonal.factor = c(12, 24, 36))`.

- lin.only:

  logical; `FALSE` (default) For fast optimization of the linear
  regression method. More robust than `lin.only = TRUE`.

- negative.values:

  logical; `FALSE` (default) If the variable can be negative, set to
  `(negative.values = TRUE)`. It will automatically select
  `(negative.values = TRUE)` if the minimum value of the `variable` is
  negative.

- obj.fn:

  expression;
  `expression(cor(predicted, actual, method = "spearman") / sum((predicted - actual)^2))`
  (default) Rank correlation / sum of squared errors is the default
  objective function. Any `expression(...)` using the specific terms
  `predicted` and `actual` can be used.

- objective:

  options: ("min", "max") `"max"` (default) Select whether to minimize
  or maximize the objective function `obj.fn`.

- linear.approximation:

  logical; `TRUE` (default) Uses the best linear output from `NNS.reg`
  to generate a nonlinear and mixture regression for comparison. `FALSE`
  is a more exhaustive search over the objective space.

- ncores:

  integer; value specifying the number of cores to be used in the
  parallelized procedure. If NULL (default), the number of cores to be
  used is equal to the number of cores of the machine - 1.

- pred.int:

  numeric \[0, 1\]; 0.95 (default) Returns the associated prediction
  intervals for the final estimate. Constructed using the maximum
  entropy bootstrap
  [NNS.meboot](https://OVVO-Financial.github.io/NNS/reference/NNS.meboot.md)
  on the final estimates.

- print.trace:

  logical; `TRUE` (default) Prints current iteration information.
  Suggested as backup in case of error, best parameters to that point
  still known and copyable!

- plot:

  logical; `FALSE` (default)

## Value

Returns a list containing:

- `$period` a vector of optimal seasonal periods

- `$weights` the optimal weights of each seasonal period between an
  equal weight or NULL weighting

- `$obj.fn` the objective function value

- `$method` the method identifying which
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  method was used.

- `$shrink` whether to use the `shrink` parameter in
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md).

- `$nns.regress` whether to smooth the variable via
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  before forecasting.

- `$bias.shift` a numerical result of the overall bias of the optimum
  objective function result. To be added to the final result when using
  the
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  with the derived parameters.

- `$errors` a vector of model errors from internal calibration.

- `$results` a vector of length `h`.

- `$lower.pred.int` a vector of lower prediction intervals per forecast
  point.

- `$upper.pred.int` a vector of upper prediction intervals per forecast
  point.

## Note

- Typically, `(training.set = 0.8 * length(variable))` is used for
  optimization. Smaller samples could use
  `(training.set = 0.9 * length(variable))` (or larger) in order to
  preserve information.

- The number of combinations will grow prohibitively large, they should
  be kept as small as possible. `seasonal.factor` containing an element
  too large will result in an error. Please reduce the maximum
  `seasonal.factor`.

- Set `(ncores = 1)` if routine is used within a parallel architecture.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
## Nonlinear NNS.ARMA period optimization using 2 yearly lags on AirPassengers monthly data
if (FALSE) { # \dontrun{
nns.optims <- NNS.ARMA.optim(AirPassengers[1:132], training.set = 120,
seasonal.factor = seq(12, 24, 6))

## To predict out of sample using best parameters:
NNS.ARMA.optim(AirPassengers[1:132], h = 12, seasonal.factor = seq(12, 24, 6))

## Incorporate any objective function from external packages (such as \code{Metrics::mape})
NNS.ARMA.optim(AirPassengers[1:132], h = 12, seasonal.factor = seq(12, 24, 6),
obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")
} # }
```
