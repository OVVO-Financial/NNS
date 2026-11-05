# NNS Stack

Prediction model using the predictions of the NNS base models
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md) as
features (i.e. meta-features) for the stacked model.

## Usage

``` r
NNS.stack(
  IVs.train,
  DV.train,
  IVs.test = NULL,
  type = NULL,
  obj.fn = expression(sum((predicted - actual)^2)),
  objective = "min",
  optimize.threshold = TRUE,
  dist = NULL,
  CV.size = NULL,
  balance = FALSE,
  ts.test = NULL,
  folds = 5,
  order = NULL,
  method = c(1, 2),
  stack = TRUE,
  dim.red.method = "cor",
  pred.int = NULL,
  status = TRUE,
  ncores = NULL,
  seed = 123L,
  .validation = NULL
)
```

## Arguments

- IVs.train:

  a vector, matrix or data frame of variables of numeric or factor data
  types.

- DV.train:

  a numeric or factor vector with compatible dimensions to
  `(IVs.train)`.

- IVs.test:

  a vector, matrix or data frame of variables of numeric or factor data
  types with compatible dimensions to `(IVs.train)`. If NULL, will use
  `(IVs.train)` as default. Columns are matched to `(IVs.train)` by name
  when the two share the same predictor names, and positionally when
  they share no names at all (as in `cbind(test.x_1, test.x_2)` against
  `cbind(x_1, x_2)`). Names that only partly overlap the training
  predictors are ambiguous and return an error.

- type:

  `NULL` (default) automatically classifies factor, character, logical
  and two-valued numeric responses. Set `type = "CLASS"` for nominal
  classification of any numeric labels, or `type = "REG"` to fit a
  numeric response as regression regardless of its number of distinct
  values. Forecasting calls use `type = "REG"`.

- obj.fn:

  expression; `expression(sum((predicted - actual)^2))` (default) Sum of
  squared errors is the default objective function. Any
  [`expression()`](https://rdrr.io/r/base/expression.html) using the
  specific terms `predicted` and `actual` can be used.

- objective:

  options: ("min", "max") `"min"` (default) Select whether to minimize
  or maximize the objective function `obj.fn`.

- optimize.threshold:

  logical; `TRUE` (default) Optimizes the threshold for the second
  first-observed class in binary classification; multiclass decisions
  use the largest class probability. If `FALSE`, returns 0.5.

- dist:

  options:(NULL, "NNS", "L1", "L2", "FACTOR") the method of distance
  calculation. `dist = NULL` is the default and selects the native
  blended NNS distance; `dist = "NNS"` is an explicit alias for the
  default. `dist = "L2"` selects Euclidean distance; `dist = "L1"`
  selects Manhattan distance; `dist = "FACTOR"` uses a frequency.

- CV.size:

  numeric \[0, 1\]; `NULL` (default) Sets the holdout size for repeated
  holdout cross-validation (one random holdout per fold). With
  `CV.size = NULL` and `folds > 1`, (stratified, for classification)
  k-fold cross-validation is used; with `CV.size = NULL` and
  `folds = 1`, a single holdout of random size between 0.2 and 0.33 is
  used.

- balance:

  logical; `FALSE` (default) Uses both up and down sampling to balance
  the classes. `type="CLASS"` required.

- ts.test:

  integer; NULL (default) Sets the length of the test set for
  time-series data; typically `2*h` parameter value from
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  or double known periods to forecast.

- folds:

  integer; `folds = 5` (default) Select the number of cross-validation
  folds.

- order:

  options: (integer, "max", NULL); `NULL` (default) Sets the order for
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md),
  where `(order = "max")` is the k-nearest neighbors equivalent, which
  is suggested for mixed continuous and discrete (unordered, ordered)
  data.

- method:

  numeric options: (1, 2); Select the NNS method to include in stack.
  `(method = 1)` selects
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md);
  `(method = 2)` selects
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  dimension reduction regression. Defaults to `method = c(1, 2)`, which
  will reduce the dimension first, then find the optimal `n.best`.

- stack:

  logical; `TRUE` (default) Uses dimension reduction to select actual
  predictor coordinates for the multivariate neighbour search; the
  dimension-reduction component separately uses its synthetic index. If
  FALSE, both components use their own full predictor designs.

- dim.red.method:

  options: ("cor", "NNS.dep", "NNS.caus", "equal", "all") method for
  determining synthetic X\* coefficients. `(dim.red.method = "cor")`
  (default) uses Spearman rank correlation for weights.
  `(dim.red.method = "NNS.dep")` uses
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  for nonlinear dependence weights, while
  `(dim.red.method = "NNS.caus")` uses
  [NNS.caus](https://OVVO-Financial.github.io/NNS/reference/NNS.caus.md)
  for directional association weights. `(dim.red.method = "all")`
  averages all methods for further feature engineering.

- pred.int:

  numeric \[0,1\]; `NULL` (default) Returns the associated prediction
  intervals with each `method`.

- status:

  logical; `TRUE` (default) Prints status update message in console.

- ncores:

  integer; value specifying the number of cores to be used in the
  parallelized subroutine
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md).
  If NULL (default), the number of cores to be used is equal to the
  number of cores of the machine - 1.

- .validation:

  Internal forecast-validation contract, normally NULL. A list with
  `IVs.train`, `DV.train`, `IVs.test`, and `DV.test` supplies the single
  chronological validation fold when `ts.test` is set, `folds = 1`, and
  `CV.size = NULL`. It requires exactly `N-ts.test` training rows and
  one through `ts.test` held-out rows; this allows scoring only observed
  future targets. Training prefixes and generated validation predictors
  must be constructed without later data. Used only for regression
  tuning; numeric two-level responses remain regression when this
  contract is supplied with type = NULL. The final model still fits all
  original training observations.

- seed:

  Optional integer random seed used for reproducible resampling, fold
  construction, and stochastic fitting steps. If \`NULL\`, the current
  random-number-generator state is used.

## Value

Returns a vector of fitted values for the dependent variable test set
for all models.

- `"NNS.reg.n.best"` returns the optimum `"n.best"` parameter for the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  multivariate regression; `NA` for a single predictor, whose univariate
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  fit has no `n.best` to tune (use `cbind(x, x)` to search one).
  `"SSE.reg"` returns the SSE for the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  multivariate regression.

- `"OBJfn.reg"` returns the best `obj.fn` value of the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  regression's tuning on the out-of-fold predictions; it is
  selection-optimized, not an estimate of out-of-sample performance.

- `"NNS.dim.red.threshold"` returns the optimum `"threshold"` from the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  dimension reduction regression.

- `"NNS.dim.red.variables"` returns the names of the predictors carried
  by the final dimension reduction regression, ranked by absolute
  coefficient. With `(stack = TRUE)` these are also the only predictors
  used by the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  multivariate regression, so they are exactly the variables behind
  `"stack"`. `NULL` when `method` does not include 2.

- `"NNS.dim.red.coefficients"` returns the named final dimension
  reduction coefficients for every encoded predictor column (zero for
  predictors not selected). `NULL` when `method` does not include 2.

- `"OBJfn.dim.red"` returns the best `obj.fn` value of the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  dimension reduction regression's tuning on the out-of-fold
  predictions; it is selection-optimized, not an estimate of
  out-of-sample performance.

- `"probability.threshold"` returns the optimum probability threshold
  for classification, else 0.5 when set to `FALSE`.

- `"reg"` returns
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  output.

- `"reg.pred.int"` returns the prediction intervals for the regression
  output.

- `"dim.red"` returns
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  dimension reduction regression output.

- `"dim.red.pred.int"` returns the prediction intervals for the
  dimension reduction regression output.

- `"stack"` returns the output of the stacked model.

- `"pred.int"` returns the prediction intervals for the stacked model.

- `"weights"` returns the weights of `"reg"` and `"dim.red"` in
  `"stack"`. For classification they are proportional to the squared
  `"OBJfn.reg"` and `"OBJfn.dim.red"` (inverse squares when
  `objective = "min"`); for regression the weight is selected by
  `obj.fn` on the out-of-fold predictions of the blend.

## Details

Classification tunes and combines class-probability vectors from the NNS
predictor cells, preserving nominal label identity. `probabilities`,
`reg.probabilities`, and `dim.red.probabilities` return the final
probability matrices. Class balancing resamples the training rows of
each fold. The neighbour count and the dimension count are chosen by
optimizing `obj.fn` on the pooled out-of-fold predictions; the default
classification objective is accuracy. As in NNS 13.0, the two models'
probability vectors are then blended with weights proportional to their
squared objective scores (inverse squares when `objective = "min"`), so
models with nearly equal scores are blended nearly equally, and the
binary probability threshold is tuned on the out-of-fold predictions of
that blend. `OBJfn.reg` and `OBJfn.dim.red` are the best scores of that
search on the same predictions, so they are selection-optimized and
overstate held-out performance; assess the final model on an independent
test set. Custom objectives receive original numeric or logical response
values, and factor/character codes as identified by `class.levels`.
Relabelling is invariant for label-invariant objectives such as
accuracy; a custom loss assigning different costs to classes can
deliberately change the result.

Named numeric dimension-reduction coefficients are matched to encoded
predictor names in both regression and classification. Factor columns
use `predictor__level` names; unnamed vectors follow encoded column
order. Classification dimension reduction uses one signed index per
class indicator. Its reported coefficients are a list of class-specific
equations. The tuned dimension count is returned as `NNS.dim.red.count`;
`NNS.dim.red.threshold` is NA because there is no single coefficient
cutoff shared by all class indices. The field `NNS.dim.red.variables` is
the union of their active encoded columns. These equations describe the
synthetic indices, not a complete saved prediction model: reconstruction
also requires training-fitted scaling, the predictor cells, their class
frequencies, and the selected smoothing controls. With `stack = TRUE`,
the multivariate classification component retains the selected number of
actual encoded coordinates ranked by maximum class-indicator dependence.
`Prediction.set` contains classes attaining the requested estimated
probability mass. Binary interval fields contain the range of this set
for compatibility; multiclass interval fields are NULL. The sets do not
assert calibrated coverage.

For regression, all dimension-reduction coefficients, normalizations,
and candidate predictions used in validation are fitted to the
corresponding training fold. The final coefficients and regressions are
refitted to the full supplied training data. Method 1 retains the
selected predictor coordinates, while method 2 evaluates their signed
synthetic index; the components can therefore represent different
regression geometry.

## Note

- Incorporate any objective function from external packages (such as
  `Metrics::mape`) via
  `NNS.stack(..., obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")`

- Like a logistic regression, the `(type = "CLASS")` setting is not
  necessary for target variable of two classes e.g. \[0, 1\]. Numeric
  class values are returned on their original scale; factor and
  character responses use codes with the returned `class.levels`
  mapping.

- Missing data should be handled prior as well using
  [na.omit](https://rdrr.io/r/stats/na.fail.html) or
  [complete.cases](https://rdrr.io/r/stats/complete.cases.html) on the
  full dataset.

## References

Viole, F. (2016) "Classification Using NNS Clustering Analysis"
[doi:10.2139/ssrn.2864711](https://doi.org/10.2139/ssrn.2864711)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
 ## Using 'iris' dataset with a stratified test set of ten rows per species
 ## (rows 141:150 alone are all virginica).
 if (FALSE) { # \dontrun{
 test.set <- c(41:50, 91:100, 141:150)
 NNS.stack(iris[-test.set, 1:4], iris[-test.set, 5], IVs.test = iris[test.set, 1:4],
 type = "CLASS", balance = TRUE)

 ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
 NNS.stack(iris[ , 1:4], iris[ , 5], type = "CLASS")
 } # }
```
