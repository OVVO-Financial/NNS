# NNS Boost

Ensemble feature-selection method using
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md) as
the base learner.

## Usage

``` r
NNS.boost(
  IVs.train,
  DV.train,
  IVs.test = NULL,
  type = NULL,
  depth = NULL,
  learner.trials = 100,
  epochs = NULL,
  CV.size = NULL,
  balance = FALSE,
  ts.test = NULL,
  threshold = NULL,
  obj.fn = expression(sum((predicted - actual)^2)),
  objective = "min",
  extreme = FALSE,
  features.only = FALSE,
  feature.importance = TRUE,
  pred.int = NULL,
  status = TRUE,
  seed = 123L
)
```

## Arguments

- IVs.train:

  a vector, matrix, or data frame of numeric, logical, character, or
  factor predictors.

- DV.train:

  a numeric, logical, character, or factor response with one value per
  row of `IVs.train`.

- IVs.test:

  a vector, matrix, or data frame with the same predictor columns as
  `IVs.train`. If `NULL`, `IVs.train` is used.

- type:

  `NULL` (default) for regression, or `"CLASS"` for classification.
  Factor, character, and logical responses automatically select
  classification.

- depth:

  integer, `NULL`, or `"max"`; passed to the `order` argument of
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md).

- learner.trials:

  positive integer; maximum number of feature subsets used to estimate
  the learner threshold. If every possible subset can be evaluated
  within this limit, all subsets are evaluated.

- epochs:

  non-negative integer; number of weighted feature subsets evaluated
  after the learner stage. Defaults to `2 * length(DV.train)`. Set to
  zero to use the surviving learner subsets directly.

- CV.size:

  numeric in `(0, 1)`; validation fraction for non-time-series data. If
  `NULL`, one value between 0.2 and 1/3 is drawn under the local seed.

- balance:

  logical; if `TRUE`, down- and up-sampling are applied only to the
  fitting portion of each split. Validation observations are never
  resampled.

- ts.test:

  positive integer smaller than the training sample size; the final
  `ts.test` observations are used as the chronological validation block.

- threshold:

  finite numeric scalar or `NULL`; objective cutoff used to retain
  feature subsets. If `NULL`, the lower quartile is used for a
  minimization objective and the upper quartile for a maximization
  objective.

- obj.fn:

  expression using the names `predicted` and `actual`. Defaults to sum
  of squared errors. For explicit classification, the untouched default
  is replaced by mean classification accuracy.

- objective:

  one of `"min"` or `"max"`; defaults to `"min"`.

- extreme:

  logical; if `TRUE`, use the best learner score rather than a quartile
  cutoff.

- features.only:

  logical; return only feature weights and frequencies.

- feature.importance:

  logical; plot up to the ten most frequently retained features.

- pred.int:

  numeric in `(0, 1)` or `NULL`; prediction interval level passed to the
  final
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  fit.

- status:

  logical; print progress messages.

- seed:

  integer or `NULL`; local random seed. The caller's RNG state is
  restored when the function exits.

## Value

A list containing `results`, `pred.int`, `feature.weights`, and
`feature.frequency`. With `features.only = TRUE`, only the last two
elements are returned.

## Note

- Numeric class labels are returned on their original scale. Factor,
  character, and logical responses retain the historical integer-code
  output.

- Categorical predictors are aligned to training levels. Unseen test
  levels cause an explicit error rather than silent recoding.

- Incorporate an objective from another package with, for example,
  `obj.fn = expression(Metrics::mape(actual, predicted))` and
  `objective = "min"`.

## References

Viole, F. (2016) "Classification Using NNS Clustering Analysis"
[doi:10.2139/ssrn.2864711](https://doi.org/10.2139/ssrn.2864711)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
a <- NNS.boost(
  iris[1:140, 1:4], iris[1:140, 5],
  IVs.test = iris[141:150, 1:4],
  epochs = 100, learner.trials = 100,
  type = "CLASS", balance = TRUE
)

mean(a$results == as.numeric(iris[141:150, 5]))
} # }
```
