# NNS meboot

Adapted maximum entropy bootstrap routine from `meboot`
<https://cran.r-project.org/package=meboot>.

## Usage

``` r
NNS.meboot(
  x,
  reps = 999,
  rho = NULL,
  type = "spearman",
  drift = TRUE,
  target_drift = NULL,
  target_drift_scale = NULL,
  trim = 0.1,
  xmin = NULL,
  xmax = NULL,
  reachbnd = TRUE,
  expand.sd = TRUE,
  force.clt = TRUE,
  scl.adjustment = FALSE,
  sym = FALSE,
  elaps = FALSE,
  digits = 6,
  colsubj,
  coldata,
  coltimes,
  ...
)
```

## Arguments

- x:

  vector of data.

- reps:

  numeric; number of replicates to generate.

- rho:

  numeric \[-1,1\] (vectorized); A `rho` must be provided, otherwise a
  blank list will be returned.

- type:

  options("spearman", "pearson", "NNScor", "NNSdep");
  `type = "spearman"`(default) dependence metric desired.

- drift:

  logical; `drift = TRUE` (default) preserves the drift of the original
  series.

- target_drift:

  numerical; `target_drift = NULL` (default) Specifies the desired drift
  when `drift = TRUE`, i.e. a risk-free rate of return.

- target_drift_scale:

  numerical; instead of calculating a `target_drift`, provide a scalar
  to the existing drift when `drift = TRUE`.

- trim:

  numeric \[0,1\]; The mean trimming proportion, defaults to
  `trim = 0.1`.

- xmin:

  numeric; the lower limit for the left tail.

- xmax:

  numeric; the upper limit for the right tail.

- reachbnd:

  logical; If `TRUE` potentially reached bounds (xmin = smallest value -
  trimmed mean and xmax = largest value + trimmed mean) are given when
  the random draw happens to be equal to 0 and 1, respectively.

- expand.sd:

  logical; If `TRUE` the standard deviation in the ensemble is expanded.
  See `expand.sd` in `meboot::meboot`.

- force.clt:

  logical; If `TRUE` the ensemble is forced to satisfy the central limit
  theorem. See `force.clt` in `meboot::meboot`.

- scl.adjustment:

  logical; If `TRUE` scale adjustment is performed to ensure that the
  population variance of the transformed series equals the variance of
  the data.

- sym:

  logical; If `TRUE` an adjustment is performed to ensure that the ME
  density is symmetric.

- elaps:

  logical; If `TRUE` elapsed time during computations is displayed.

- digits:

  integer; 6 (default) number of digits to round output to.

- colsubj:

  numeric; the column in `x` that contains the individual index. It is
  ignored if the input data `x` is not a `pdata.frame` object.

- coldata:

  numeric; the column in `x` that contains the data of the variable to
  create the ensemble. It is ignored if the input data `x` is not a
  `pdata.frame` object.

- coltimes:

  numeric; an optional argument indicating the column that contains the
  times at which the observations for each individual are observed. It
  is ignored if the input data `x` is not a `pdata.frame` object.

- ...:

  possible argument `fiv` to be passed to `expand.sd`.

## Value

Returns the following row names in a matrix:

- x original data provided as input.

- replicates maximum entropy bootstrap replicates.

- ensemble average observation over all replicates.

- xx sorted order stats (xx\[1\] is minimum value).

- z class intervals limits.

- dv deviations of consecutive data values.

- dvtrim trimmed mean of dv.

- xmin data minimum for ensemble=xx\[1\]-dvtrim.

- xmax data x maximum for ensemble=xx\[n\]+dvtrim.

- desintxb desired interval means.

- ordxx ordered x values.

- kappa scale adjustment to the variance of ME density.

- elaps elapsed time.

## Note

Vectorized `rho` and `drift` parameters will not vectorize both
simultaneously. Also, do not specify `target_drift = NULL`.

## References

- Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank
  Correlations in Maximum Entropy Bootstrap and Improved Monte Carlo
  Simulations.
  [doi:10.2139/ssrn.3621614](https://doi.org/10.2139/ssrn.3621614)

- Vinod, H.D. (2013), Maximum Entropy Bootstrap Algorithm Enhancements.
  [doi:10.2139/ssrn.2285041](https://doi.org/10.2139/ssrn.2285041)

- Vinod, H.D. (2006), Maximum Entropy Ensembles for Time Series
  Inference in Economics, *Journal of Asian Economics*, **17**(6), pp.
  955-978.

- Vinod, H.D. (2004), Ranking mutual funds using unconventional utility
  theory and stochastic dominance, *Journal of Empirical Finance*,
  **11**(3), pp. 353-377.

## Examples

``` r
if (FALSE) { # \dontrun{
# To generate an orthogonal rank correlated time-series to AirPassengers
boots <- NNS.meboot(AirPassengers, reps = 100, rho = 0, xmin = 0)

# Verify correlation of replicates ensemble to original
cor(boots["ensemble",]$ensemble, AirPassengers, method = "spearman")

# Plot all replicates
matplot(boots["replicates",]$replicates , type = 'l')

# Plot ensemble
lines(boots["ensemble",]$ensemble, lwd = 3)

# Plot original
lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red")

### Vectorized drift with a single rho
boots <- NNS.meboot(AirPassengers, reps = 10, rho = 0, xmin = 0, target_drift = c(1,7))
matplot(do.call(cbind, boots["replicates", ]), type = "l")
lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red")

### Vectorized rho with a single target drift
boots <- NNS.meboot(AirPassengers, reps = 10, rho = c(0, .5, 1), xmin = 0, target_drift = 3)
matplot(do.call(cbind, boots["replicates", ]), type = "l")
lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red")

### Vectorized rho with a single target drift scale
boots <- NNS.meboot(AirPassengers, reps = 10, rho = c(0, .5, 1), xmin = 0, target_drift_scale = 0.5)
matplot(do.call(cbind, boots["replicates", ]), type = "l")
lines(1:length(AirPassengers), AirPassengers, lwd = 3, col = "red") 
} # }
```
