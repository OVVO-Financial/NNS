# NNS 13.3 (development changes)

## Statistical targets and interpretation

* Preserve the distinction between frequency probabilities (degree zero),
  directional deviation-area shares (positive degrees), and their induced laws.
  Correct one-sided probability bounds, divergent-moment labels, finite-moment
  conditions, prediction-interval claims, and the scope of dependence/causation.
* ANOVA retains its separate mean, median, and distribution-profile targets.
  Median comparisons now handle CDF jumps, and degenerate area-balance cases
  use explicit agreement semantics. Certainty remains a descriptive score;
  inferential calibration must match the target and sampling design.
* `NNS.FSD(type = "continuous")` keeps the same numerical degree-one criterion
  but labels it Omega/expectile order. Its result strings now use `OMEGA`, not
  `FSD`; callers comparing strings should update. Discrete FSD strings and
  numeric dominance indicators/set/cluster structures are unchanged.

## Regression, classification, and derivatives

* `NNS.reg` and `NNS.stack` accept `type = "REG"` to select numeric
  regression explicitly. Forecasting, derivative, and hazard-smoothing helpers
  use this contract even when a response has only two values; they retain
  fractional predictions on the original scale.
* Nominal classification retains class probabilities in predictor cells and
  combines those probabilities before choosing labels. Numeric class labels are
  not interpolation coordinates. Probability outputs and descriptive prediction
  sets expose the target; factor/character output codes use `class.levels`.
* Binary numeric coordinates retain distinct terminal cells in both R and
  native regression paths, including dummy columns in mixed predictor designs.
* Stack's multivariate component searches selected predictor coordinates; its
  dimension-reduction component separately uses the synthetic index. Complexity
  selection uses both response orientations without changing the underlying
  directional statistics. These corrections can change fitted predictions.
* Distance kernels use consistent exponential rates and sample rank scales
  across scalar, bulk, and parallel paths. Partial derivatives use admissible
  local support and explicit boundary differences. Unsupported multivariate
  `NNS.CDF(type = "hazard")` requests now fail instead of returning CDF/survival
  as if it were density/survival.
* Gravity retains its location blend with corrected quartiles and observed bin
  centers. Mode retains distinct continuous and class behavior. Changes to
  these defaults can affect regression, optimization, and forecasts.

## Forecasting and scenarios

* ARMA validation support is selected from training observations only. Nonlinear
  shrinkage follows the same intended mean blend as linear shrinkage, and a
  single shrunk lag refits recursively. Long-horizon runtime and path bands can
  therefore differ from earlier releases.
* The optimizer's objective describes its returned internal validation forecast.
  The default bias policy estimates a shift on an earlier validation portion and
  accepts it on a later portion; `bias.correction = "legacy"` reproduces the old
  same-window policy and `"none"` disables it. Base-model tuning still requires
  independent outer evaluation. Raw/corrected forecasts and bias metadata are
  available for assessment.
* VAR keeps its two-stage forecast design and uses prefix-generated predictors
  for chronological validation by default. `validation = "observed"` is the
  compatibility option. Known contemporaneous outcomes and imputed future
  targets are not used as deployment-validation inputs or ground truth.
* Seasonality adds absolute dispersion for interval-scale data while retaining
  the historical automatic CV/Gini default. Selected predictive periods need not be fundamental
  physical periods; permutation interpretations require an appropriate null.
* Monte Carlo exponent grids stay within requested bounds and reject invalid
  exponents. Constant/singleton superiority groups and exactly zero bootstrap
  residuals are handled explicitly. Bootstrap residual expansion is separately
  configurable. The default forced-mean scale is an iid reference, not a
  universal serial-dependence uncertainty calibration.

## Numerical, documentation, and reproducibility repairs

* Fix mixed-degree co-partial moments, validate paired lengths before native
  parallel work, and make integer-degree VaR root stopping scale invariant.
* Fix target hazard evaluation, stochastic-dominance cluster ordering, residual
  clusters, efficient-set order dependence, and caller column-name mutation.
* Replace unsupported example rankings and inference claims with explicit
  targets and corrected executable workflows. Forecast examples reserve test
  outcomes; fitted-contrast uncertainty includes model refitting. Add revised
  editable companions for the historical PDFs, preserve originals, and record
  that the original editable sources were not recovered.
