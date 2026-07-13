# NNS Stack

Cross-validated ensemble of the full multivariate and
synthetic-dimension
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
models.

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
  dist = "L2",
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
  seed = 123L
)
```

## Arguments

- IVs.train:

  a vector, matrix, or data frame of numeric, logical, character,
  factor, Date, or date-time predictors.

- DV.train:

  a numeric, logical, character, or factor response with one value per
  row of `IVs.train`.

- IVs.test:

  a vector, matrix, or data frame with the same predictors as
  `IVs.train`. If `NULL`, `IVs.train` is used.

- type:

  `NULL` (default) for regression or `"CLASS"` for classification.
  Factor, character, logical, and two-level numeric responses
  automatically select classification.

- obj.fn:

  an expression using `predicted` and `actual`. Sum of squared errors is
  the regression default. For classification, the untouched default is
  replaced by mean classification accuracy.

- objective:

  one of `"min"` or `"max"`.

- optimize.threshold:

  logical; optimize the class-rounding threshold from out-of-fold
  predictions. If `FALSE`, use 0.5.

- dist:

  distance option. The corrected implementation currently accepts only
  `"L2"`, because the production multivariate `NNS.reg` path does not
  presently implement distinct L1, DTW, or FACTOR estimators.

- CV.size:

  optional validation fraction in `(0, 1)`. If supplied, `folds`
  repeated stratified/random holdouts are used. If `NULL`, disjoint
  k-fold cross-validation is used.

- balance:

  logical; balance only each fitting partition and the final fitting
  data. Validation observations are never resampled.

- ts.test:

  positive integer; validation-block length for chronological
  rolling-origin cross-validation. The final block always contains the
  most recent observations.

- folds:

  positive integer; number of ordinary or rolling-origin folds.

- order:

  integer, `"max"`, or `NULL`; passed unchanged to
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md).

- method:

  any unique combination of `1` and `2`. Method 1 is the full
  multivariate `NNS.reg` model with cross-validated `n.best`; Method 2
  is a synthetic X\* dimension-reduction model.

- stack:

  logical; when both methods are requested, use fold-local and
  training-only X\* as Method 1's input. If `FALSE`, Method 1 uses the
  full independently encoded predictor matrix.

- dim.red.method:

  one of `"cor"`, `"NNS.dep"`, `"NNS.caus"`, `"equal"`, `"all"`, or a
  numeric coefficient vector aligned to the encoded design columns.

- pred.int:

  numeric in `(0, 1)` or `NULL`; prediction interval level for the final
  component fits.

- status:

  logical; print progress messages.

- ncores:

  positive integer or `NULL`; native thread count.

- seed:

  non-negative integer or `NULL`; local random seed. The caller's
  random-number state is restored on exit.

## Value

A list retaining the historical fields:

- `OBJfn.reg`: selected Method 1 out-of-fold objective.

- `NNS.reg.n.best`: selected Method 1 `n.best`.

- `probability.threshold`: threshold optimized directly on the
  out-of-fold weighted ensemble, or 0.5 for regression.

- `OBJfn.dim.red`: selected Method 2 out-of-fold objective.

- `NNS.dim.red.threshold`: full-data coefficient-magnitude cutoff
  corresponding to the selected active-dimension count.

- `reg`, `dim.red`, and `stack`: final predictions.

- component and stacked prediction intervals.

Classification predictions are returned as numeric class codes, matching
the historical NNS.stack interface. Additional fields `weights` and
`class.levels` report the out-of-fold blend and the code-to-label map.

## Note

Categorical encoding and min-max normalization are fitted on each
training partition only and then applied unchanged to its validation
partition. The final transformations are fitted on the complete training
data only; external test observations never affect training scales.

## References

Viole, F. (2016) "Classification Using NNS Clustering Analysis"
[doi:10.2139/ssrn.2864711](https://doi.org/10.2139/ssrn.2864711)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- NNS.stack(
  iris[1:140, 1:4], iris[1:140, 5],
  IVs.test = iris[141:150, 1:4],
  type = "CLASS", balance = TRUE
)
fit$stack
} # }
```
