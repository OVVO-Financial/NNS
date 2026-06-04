# NNS VAR

Nonparametric vector autoregressive model incorporating
[NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
estimates of variables into
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md) for
a multi-variate time-series forecast.

## Usage

``` r
NNS.VAR(
  variables,
  h,
  tau = 1,
  dim.red.method = "cor",
  naive.weights = TRUE,
  obj.fn = expression(mean((predicted - actual)^2)/(NNS::Co.LPM(1, predicted, actual,
    target_x = mean(predicted), target_y = mean(actual)) + NNS::Co.UPM(1, predicted,
    actual, target_x = mean(predicted), target_y = mean(actual)))),
  objective = "min",
  status = TRUE,
  ncores = NULL,
  nowcast = FALSE
)
```

## Arguments

- variables:

  a numeric matrix or data.frame of contemporaneous time-series to
  forecast.

- h:

  integer; 1 (default) Number of periods to forecast. `(h = 0)` will
  return just the interpolated and extrapolated values.

- tau:

  positive integer \[ \> 0\]; 1 (default) Number of lagged observations
  to consider for the time-series data. Vector for single lag for each
  respective variable or list for multiple lags per each variable.

- dim.red.method:

  options: ("cor", "NNS.dep", "NNS.caus", "all") method for reducing
  regressors via
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md).
  `(dim.red.method = "cor")` (default) uses standard linear correlation
  for dimension reduction in the lagged variable matrix.
  `(dim.red.method = "NNS.dep")` uses
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  for nonlinear dependence weights, while
  `(dim.red.method = "NNS.caus")` uses
  [NNS.caus](https://OVVO-Financial.github.io/NNS/reference/NNS.caus.md)
  for causal weights. `(dim.red.method = "all")` averages all methods
  for further feature engineering.

- naive.weights:

  logical; `TRUE` (default) Equal weights applied to univariate and
  multivariate outputs in ensemble. `FALSE` will apply weights based on
  the number of relevant variables detected.

- obj.fn:

  expression;
  `expression(mean((predicted - actual)^2)) / (Sum of NNS Co-partial moments)`
  (default) MSE / co-movements is the default objective function. Any
  `expression(...)` using the specific terms `predicted` and `actual`
  can be used.

- objective:

  options: ("min", "max") `"min"` (default) Select whether to minimize
  or maximize the objective function `obj.fn`.

- status:

  logical; `TRUE` (default) Prints status update message in console.

- ncores:

  integer; value specifying the number of cores to be used in the
  parallelized subroutine
  [NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md).
  If NULL (default), the number of cores to be used is equal to the
  number of cores of the machine - 1.

- nowcast:

  logical; `FALSE` (default) internal call for frequency alignment in
  downstream nowcasting applications.

## Value

Returns the following matrices of forecasted variables:

- `"interpolated_and_extrapolated"` Returns a `data.frame` of the linear
  interpolated and
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  extrapolated values to replace `NA` values in the original `variables`
  argument. This is required for working with variables containing
  different frequencies, e.g. where `NA` would be reported for
  intra-quarterly data when indexed with monthly periods.

- `"relevant_variables"` Returns the relevant variables from the
  dimension reduction step.

- `"univariate"` Returns the univariate
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  forecasts.

- `"multivariate"` Returns the multi-variate
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  forecasts.

- `"ensemble"` Returns the ensemble of both `"univariate"` and
  `"multivariate"` forecasts.

## Note

- `"Error in { : task xx failed -}"` should be re-run with
  `NNS.VAR(..., ncores = 1)`.

- Not recommended for factor variables, even after transformed to
  numeric.
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  is better suited for factor or binary regressor extrapolation.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric
Vector Autoregression Using NNS"
[doi:10.2139/ssrn.3489550](https://doi.org/10.2139/ssrn.3489550)

Viole, F. (2020) "NOWCASTING with NNS"
[doi:10.2139/ssrn.3589816](https://doi.org/10.2139/ssrn.3589816)

Viole, F. (2019) "Forecasting Using NNS"
[doi:10.2139/ssrn.3382300](https://doi.org/10.2139/ssrn.3382300)

Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
[doi:10.1007/s10614-017-9713-5](https://doi.org/10.1007/s10614-017-9713-5)

Vinod, H. and Viole, F. (2018) "Clustering and Curve Fitting by Line
Segments"
[doi:10.20944/preprints201801.0090.v1](https://doi.org/10.20944/preprints201801.0090.v1)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r

 if (FALSE) { # \dontrun{
 ####################################################
 ### Standard Nonparametric Vector Autoregression ###
 ####################################################

 set.seed(123)
 x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
 A <- cbind(x = x, y = y, z = z)

 ### Using lags 1:4 for each variable
 NNS.VAR(A, h = 12, tau = 4, status = TRUE)

 ### Using lag 1 for variable 1, lag 3 for variable 2 and lag 3 for variable 3
 NNS.VAR(A, h = 12, tau = c(1,3,3), status = TRUE)

 ### Using lags c(1,2,3) for variables 1 and 3, while using lags c(4,5,6) for variable 2
 NNS.VAR(A, h = 12, tau = list(c(1,2,3), c(4,5,6), c(1,2,3)), status = TRUE)

 ### PREDICTION INTERVALS
 # Store NNS.VAR output
 nns_estimate <- NNS.VAR(A, h = 12, tau = 4, status = TRUE)

 # Create bootstrap replicates using NNS.meboot
 replicates <- NNS.meboot(nns_estimate$ensemble[,1], rho = seq(-1,1,.25))["replicates",]
 replicates <- do.call(cbind, replicates)

 # Apply UPM.VaR and LPM.VaR for desired prediction interval...95 percent illustrated
 # Tail percentage used in first argument per {LPM.VaR} and {UPM.VaR} functions
 lower_CIs <- apply(replicates, 1, function(z) LPM.VaR(0.025, 0, z))
 upper_CIs <- apply(replicates, 1, function(z) UPM.VaR(0.025, 0, z))

 # View results
 cbind(nns_estimate$ensemble[,1], lower_CIs, upper_CIs)


 #########################################
 ### NOWCASTING with Mixed Frequencies ###
 #########################################

 library(Quandl)
 econ_variables <- Quandl(c("FRED/GDPC1", "FRED/UNRATE", "FRED/CPIAUCSL"),type = 'ts',
                          order = "asc", collapse = "monthly", start_date = "2000-01-01")

 ### Note the missing values that need to be imputed
 head(econ_variables)
 tail(econ_variables)


 NNS.VAR(econ_variables, h = 12, tau = 12, status = TRUE)
 } # }
```
