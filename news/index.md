# Changelog

## NNS 13.3 (development changes)

CRAN release: 2026-10-30

### Statistical targets and interpretation

- Preserve the distinction between frequency probabilities (degree
  zero), directional deviation-area shares (positive degrees), and their
  induced laws. Correct one-sided probability bounds, divergent-moment
  labels, finite-moment conditions, prediction-interval claims, and the
  scope of dependence/causation.
- ANOVA retains its separate mean, median, and distribution-profile
  targets. Median comparisons now handle CDF jumps, and degenerate
  area-balance cases use explicit agreement semantics. Certainty remains
  a descriptive score; inferential calibration must match the target and
  sampling design.
- `NNS.FSD(type = "continuous")` keeps the same numerical degree-one
  criterion but labels it Omega/expectile order. Its result strings now
  use `OMEGA`, not `FSD`; callers comparing strings should update.
  Discrete FSD strings and numeric dominance indicators/set/cluster
  structures are unchanged.

### Regression, classification, and derivatives

- `NNS.reg` and `NNS.stack` accept `type = "REG"` to select numeric
  regression explicitly. Forecasting, derivative, and hazard-smoothing
  helpers use this contract even when a response has only two values;
  they retain fractional predictions on the original scale.
- Nominal classification retains class probabilities in predictor cells
  and combines those probabilities before choosing labels. Numeric class
  labels are not interpolation coordinates. Probability outputs and
  descriptive prediction sets expose the target; factor/character output
  codes use `class.levels`.
- Binary numeric coordinates retain distinct terminal cells in both R
  and native regression paths, including dummy columns in mixed
  predictor designs.
- Stack’s multivariate component searches selected predictor
  coordinates; its dimension-reduction component separately uses the
  synthetic index. Complexity selection uses both response orientations
  without changing the underlying directional statistics. These
  corrections can change fitted predictions.
- Distance kernels use consistent exponential rates and sample rank
  scales across scalar, bulk, and parallel paths. Partial derivatives
  use admissible local support and explicit boundary differences.
  Unsupported multivariate `NNS.CDF(type = "hazard")` requests now fail
  instead of returning CDF/survival as if it were density/survival.
- Stack’s multivariate component keeps the same coordinates in every
  fold, so a predictor that is constant in one fold no longer stops the
  search. A declared factor level with no training rows is a valid test
  value, a single test row with a factor predictor is accepted, and
  categorical `NNS.reg` honours an explicit `n.best`.
- Nominal classification reuses each fold’s cells and sorted distances
  for every neighbour count and skips fitted classes in point-only
  calls. It is still slower than the former class-code regression: for
  three classes, two predictors and `folds = 2`, `NNS.stack` takes 1.3 s
  and `NNS.boost` 1.3 s at n = 1000 (0.9 s and 0.12 s before), and 5.8 s
  and 8.0 s at n = 5000 (2.6 s and 1.2 s before).
- Gravity retains its location blend with corrected quartiles and
  observed bin centers. Mode retains distinct continuous and class
  behavior. Changes to these defaults can affect regression,
  optimization, and forecasts.

### Forecasting and scenarios

- ARMA validation support is selected from training observations only.
  Nonlinear shrinkage follows the same intended mean blend as linear
  shrinkage, and a single shrunk lag refits recursively. Long-horizon
  runtime and path bands can therefore differ from earlier releases.
- The optimizer’s objective describes its returned internal validation
  forecast. `NNS.ARMA.optim` no longer adds the gravity bias shift by
  default (`bias.correction = "none"`). The former rule estimated the
  shift on the validation window and accepted it whenever it improved
  the objective on that same window; it is removed. On rolling-origin
  forecasts of nine series the shift made held-out MASE worse on eight,
  under the former rule and under `bias.correction = "split"`, which
  estimates it on the first half of the window and accepts it on the
  second and remains an option. Base-model tuning still requires
  independent outer evaluation. Raw/corrected forecasts and bias
  metadata are available for assessment.
- `NNS.ARMA.optim` and the nonlinear `NNS.ARMA` methods run faster. The
  nonlinear step’s weighted rise and run were means of a vector about
  k^3/3 long for k regression points; they are now computed in closed
  form, which agrees with the former value to rounding. The optimizer
  computes each validation forecast once, partitions no longer build
  their data frames through R, and small partial-moment kernels run on
  the calling thread; these leave results bit-identical. On a synthetic
  monthly series of 8000 observations the optimizer took 13 s, against
  126 s in 13.3 and 210 s before this change; the AirPassengers example
  takes 0.10 s, against 0.30 s and 0.39 s.
- Long series forecast faster with the same results. When every lag is
  at least `h`, no forecast reaches a lag’s component series, so
  `NNS.ARMA`’s linear method forecasts several lags directly, as it did
  one lag, and combines them step by step as before; `NNS.ARMA.optim`’s
  linear search makes these calls. `NNS.seas` scores its lags on several
  threads, and its fundamental-period jackknife sorts all components at
  once. The tides example (222,000 six-minute observations, 200
  candidate periods) took 11 min on one core with 13.0 and 13 min before
  this change; it now takes 1.0 min, or 0.7 min with the default cores.
- `NNS.VAR` fills the missing values of mixed-frequency panels as
  documented: values before a variable’s last observation are
  interpolated linearly between its neighbouring observations, and
  trailing values are `NNS.ARMA` forecasts from the completed series.
  When a variable had trailing `NA`s, every missing value, interior ones
  included, came from an `NNS.stack` fit on the time index, which lags a
  trending series: quarterly real GDP in a monthly panel was filled up
  to 3% below its observed quarters in 2021-2022. On one- and two-month
  trailing gaps in 12 FRED series (origins 2016-2022), `NNS.ARMA` had
  the lowest median error of five fills and the index fit the highest,
  six to ten times the others for quarterly series.
- VAR keeps its two-stage forecast design and uses prefix-generated
  predictors for chronological validation by default.
  `validation = "observed"` is the compatibility option. Known
  contemporaneous outcomes and imputed future targets are not used as
  deployment-validation inputs or ground truth.
- `NNS.seas` reports the fundamental period as `best.period`: the
  smallest period that divides the lowest-scoring lag (or nearly divides
  it, for a non-integer period) and scores within two jackknife standard
  errors of it, from paired delete-one differences. The lowest score
  often fell on a multiple of the true period; a monthly sine with noise
  now returns 12 instead of 36 to 108, and `NNS.ARMA`’s automatic lag
  follows. The period is listed first in `all.periods` and `periods`.
  `NNS.ARMA` does not run the permutation test; pass
  `NNS.seas(x, p.value = TRUE)$best.period` as `seasonal.factor` for a
  significance-filtered lag.
- Seasonality adds absolute dispersion for interval-scale data while
  retaining the historical automatic CV/Gini default. Selected
  predictive periods need not be fundamental physical periods;
  permutation interpretations require an appropriate null.
- Monte Carlo exponent grids stay within requested bounds and reject
  invalid exponents. Constant/singleton superiority groups and exactly
  zero bootstrap residuals are handled explicitly. Bootstrap residual
  expansion is separately configurable. The default forced-mean scale is
  an iid reference, not a universal serial-dependence uncertainty
  calibration.

### Numerical, documentation, and reproducibility repairs

- Fix mixed-degree co-partial moments, validate paired lengths before
  native parallel work, and make integer-degree VaR root stopping scale
  invariant.
- Fix target hazard evaluation, stochastic-dominance cluster ordering,
  residual clusters, efficient-set order dependence, and caller
  column-name mutation.
- Replace unsupported example rankings and inference claims with
  explicit targets and corrected executable workflows. Forecast examples
  reserve test outcomes; fitted-contrast uncertainty includes model
  refitting. Add revised editable companions for the historical PDFs,
  preserve originals, and record that the original editable sources were
  not recovered.
