# NNS ARMA

Autoregressive model incorporating nonlinear regressions of component
series.

## Usage

``` r
NNS.ARMA(
  variable,
  h = 1,
  training.set = NULL,
  seasonal.factor = TRUE,
  weights = NULL,
  best.periods = 1,
  modulo = NULL,
  mod.only = TRUE,
  negative.values = FALSE,
  method = "nonlin",
  dynamic = FALSE,
  shrink = FALSE,
  plot = TRUE,
  seasonal.plot = TRUE,
  pred.int = NULL
)
```

## Arguments

- variable:

  a numeric vector.

- h:

  integer; 1 (default) Number of periods to forecast.

- training.set:

  numeric; `NULL` (default) Sets the number of variable observations

  `(variable[1 : training.set])` to monitor performance of forecast over
  in-sample range.

- seasonal.factor:

  logical or integer(s); `TRUE` (default) Automatically selects the best
  seasonal lag from the seasonality test,
  [NNS.seas](https://OVVO-Financial.github.io/NNS/reference/NNS.seas.md)'s
  `best.period` without its permutation test. For a
  significance-filtered lag, run `NNS.seas(variable, p.value = TRUE)`
  first and pass its `best.period` (or `periods`) here; `best.period` is
  1 when no period is significant. To use weighted average of all
  seasonal lags set to `(seasonal.factor = FALSE)`. Otherwise, directly
  input known frequency integer lag to use, i.e.
  `(seasonal.factor = 12)` for monthly data. Multiple frequency integers
  can also be used, i.e. `(seasonal.factor = c(12, 24, 36))`

- weights:

  numeric or `"equal"`; `NULL` (default) sets the weights of the
  `seasonal.factor` vector when specified as integers. If
  `(weights = NULL)` each `seasonal.factor` is weighted on its
  [NNS.seas](https://OVVO-Financial.github.io/NNS/reference/NNS.seas.md)
  result and number of observations it contains, else an `"equal"`
  weight is used. Numeric weights must have one value per seasonal lag
  and are treated as relative weights (rescaled to sum to 1). With
  `(dynamic = TRUE)` the lags are re-selected at every step, so only
  `NULL` or `"equal"` weights are accepted.

- best.periods:

  integer; 1 (default) used in conjunction with
  `(seasonal.factor = FALSE)`, uses the `best.periods` number of
  detected seasonal lags instead of `ALL` lags when
  `(seasonal.factor = FALSE, best.periods = NULL)`.

- modulo:

  integer(s); NULL (default) Adds the multiples of `modulo` on either
  side of each detected period (period 7 with `modulo = 5` gives 5 and
  10).

- mod.only:

  logical; `TRUE` (default) Limits the number of seasonal periods
  returned to the specified `modulo`.

- negative.values:

  logical; `FALSE` (default) If the variable can be negative, set to
  `(negative.values = TRUE)`. If there are negative values within the
  training data `variable[1 : training.set]`, `negative.values` will
  automatically be detected.

- method:

  options: ("lin", "nonlin", "both", "means"); `"nonlin"` (default) To
  select the regression type of the component series, select
  `(method = "both")` where both linear and nonlinear estimates are
  generated. To use a nonlinear regression, set to
  `(method = "nonlin")`; to use a linear regression set to
  `(method = "lin")`. Means for each subset are returned with
  `(method = "means")`.

- dynamic:

  logical; `FALSE` (default) To update the seasonal factor with each
  forecast point, set to `(dynamic = TRUE)`. The default is
  `(dynamic = FALSE)` to retain the original seasonal factor from the
  inputted variable for all ensuing `h`.

- shrink:

  logical; `FALSE` (default) Ensembles forecasts with
  `method = "means"`.

- plot:

  logical; `TRUE` (default) Returns the plot of all periods exhibiting
  seasonality and the `variable` level reference in upper panel. Lower
  panel returns original data and forecast.

- seasonal.plot:

  logical; `TRUE` (default) Adds the seasonality plot above the
  forecast. Will be set to `FALSE` if no seasonality is detected or
  `seasonal.factor` is set to an integer value.

- pred.int:

  numeric \[0, 1\]; `NULL` (default) Plots and returns lower and upper
  bands around the forecasts. The forecast path itself is resampled with
  the maximum entropy bootstrap
  ([NNS.MC](https://OVVO-Financial.github.io/NNS/reference/NNS.MC.md)
  over `rho` from -1 to 1), the bounds at each step are
  [LPM.VaR](https://OVVO-Financial.github.io/NNS/reference/LPM.VaR.md)
  and
  [UPM.VaR](https://OVVO-Financial.github.io/NNS/reference/UPM.VaR.md)
  of the replicates at probability `(1 - pred.int) / 2`, widened by the
  spread of the per-lag linear estimates; a one-point (`h = 1`) or
  constant path has no dispersion to resample, so its band is the path
  widened by that spread alone (zero width with a single seasonal lag);
  [NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md)
  gives intervals from validation errors. The band therefore describes
  the dispersion of the forecast path, not the uncertainty of the
  outcome at each horizon: it does not widen with the horizon, and it is
  not calibrated to contain future values with probability `pred.int`.
  [NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md)
  derives its interval from validation errors instead.

## Value

Returns a vector of forecasts of length `(h)` if no `pred.int`
specified. Else, returns a `data.table` with the forecasts as well as
the lower and upper forecast-path bands per forecast point (see
`pred.int`).

## Note

For monthly data series, increased accuracy may be realized from forcing
seasonal factors to multiples of 12. For example, if the best periods
reported are: {37, 47, 71, 73} use `(seasonal.factor = c(36, 48, 72))`.

`(seasonal.factor = FALSE)` can be a very computationally expensive
exercise due to the number of seasonal periods detected.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Viole, F. (2019) "Forecasting Using NNS"
[doi:10.2139/ssrn.3382300](https://doi.org/10.2139/ssrn.3382300)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r

## Nonlinear NNS.ARMA using AirPassengers monthly data and 12 period lag
if (FALSE) { # \dontrun{
NNS.ARMA(AirPassengers, h = 45, training.set = 100, seasonal.factor = 12, method = "nonlin")

## Linear NNS.ARMA using AirPassengers monthly data and 12, 24, and 36 period lags
NNS.ARMA(AirPassengers, h = 45, training.set = 120, seasonal.factor = c(12, 24, 36), method = "lin")

## Nonlinear NNS.ARMA using AirPassengers monthly data and 2 best periods lag
NNS.ARMA(AirPassengers, h = 45, training.set = 120, seasonal.factor = FALSE, best.periods = 2)
} # }
```
