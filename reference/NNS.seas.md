# NNS Seasonality Test

Seasonality test based on the coefficient of variation for the variable
and lagged component series. A linear trend is removed first
(`detrend = TRUE`), keeping the mean of the variable, so that a trend
does not mask the period. For each lag `m`, the variable is split into
its `m` component series `variable[seq(j, n, m)]`, `j = 1, ..., m`, and
lag `m` is reported when the length-weighted mean of the component
coefficients of variation is below the coefficient of variation of the
variable. Dispersion is the Gini mean difference, the degree-1 lower
plus upper partial moment about each observation in turn, which is
unbiased for any distribution at any component length. The coefficient
of variation is used when the (detrended) variable stays on one side of
zero; otherwise dispersions are compared directly. Lags are scanned up
to `max.period`, by default `n / 5`, so every component series has at
least 5 observations. Periods added through `modulo` are kept up to
`n / 2`. On a series without seasonality about half of the lags pass by
chance, so treat the result as candidate periods and validate them, for
example with
[NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md),
or set `p.value = TRUE`. A result of 1 signifies no seasonality present.
A multiple of the true period scores as well in expectation, and its
shorter components make its score noisier, so the lowest score often
falls on a multiple. `best.period` is therefore the fundamental period:
the smallest period that divides the lowest-scoring lag, or nearly
divides it for a non-integer period, and whose score exceeds the lowest
by at most two jackknife standard errors of their difference, from
paired delete-one scores. It is listed first in `all.periods` and
`periods`, followed by the other periods in order of score.

## Usage

``` r
NNS.seas(
  variable,
  modulo = NULL,
  mod.only = TRUE,
  plot = TRUE,
  detrend = TRUE,
  max.period = NULL,
  p.value = FALSE,
  nperm = 100L,
  dispersion = c("auto", "absolute")
)
```

## Arguments

- variable:

  a numeric vector.

- modulo:

  integer(s); NULL (default) Adds the multiples of `modulo` on either
  side of each detected period (period 7 with `modulo = 5` gives 5 and
  10).

- mod.only:

  logical; `TRUE` (default) Limits the number of seasonal periods
  returned to the specified `modulo`.

- plot:

  logical; `TRUE` (default) Returns the plot of all periods exhibiting
  seasonality and the variable level reference.

- detrend:

  logical; `TRUE` (default) removes the variable's linear trend, keeping
  its mean, before the lags are scored; an exactly linear variable is
  then constant and has no seasonal period.

- max.period:

  integer; `NULL` (default) scans lags up to `n / 5`. A smaller value
  caps the scan, whose cost grows with `n` times the number of lags.

- p.value:

  logical; `FALSE` (default). If `TRUE`, each period's dispersion ratio
  is compared with the best ratio over all lags in `nperm` shuffled
  copies of the (detrended) variable, which remove the original temporal
  ordering.
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)'s
  automatic lag selection does not run this test; for a
  significance-filtered lag, pass
  `NNS.seas(variable, p.value = TRUE)$best.period` (or `$periods`) to
  its `seasonal.factor`. `all.periods` gains a `p.value` column, and
  `periods` and `best.period` keep only the periods with
  `p.value <= 0.05` (`1` if none); each copy contributes its best lag,
  accounting for scanning multiple lags under an exchangeable null. A
  serially dependent nonseasonal series need not satisfy that null. Uses
  the random number generator.

- nperm:

  integer; 100 (default) number of shuffled copies when
  `p.value = TRUE`.

- dispersion:

  `"auto"` (default) retains relative CV scoring on one-sided data and
  absolute Gini dispersion otherwise. `"absolute"` always compares Gini
  dispersions, by centering the series first; use this
  translation-invariant option for interval-scale data such as
  temperature.

## Value

Returns a matrix of all periods exhibiting less coefficient of variation
than the variable with `"all.periods"`; and the fundamental period with
`"best.period"` (see Description); as well as a vector of `"periods"`
for easy call into
[NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md).
If no seasonality is detected, `all.periods` holds a single row with
`Period = 1` and a coefficient of variation equal to the variable's
(`Period = 0` for fewer than 5 observations). The coefficients of
variation are those of the detrended variable when `detrend = TRUE`.

## Note

The historical column names `Coefficient.of.Variation` and
`Variable.Coefficient.of.Variation` are retained in absolute mode; their
values then contain absolute Gini dispersions.

`best.period` need not be the lowest-scoring lag:
`with(all.periods, Period[which.min(Coefficient.of.Variation)])` returns
that lag. The jackknife behind the two-standard-error rule treats
observations as exchangeable and does not adjust for serial dependence.
The minimum five observations per component is a reliability default;
`modulo` can supply longer domain-informed candidates. Shuffling removes
serial dependence as well as seasonality, so the optional permutation
comparison is not calibrated against every dependent nonseasonal
process.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)

## To call strongest period based on coefficient of variation:
NNS.seas(x, plot = FALSE)$best.period

## Using modulos for logical seasonal inference:
NNS.seas(x, modulo = c(2,3,5,7), plot = FALSE)

## Keep only periods that beat shuffled copies of the series:
NNS.seas(x, plot = FALSE, p.value = TRUE)$best.period
} # }
```
