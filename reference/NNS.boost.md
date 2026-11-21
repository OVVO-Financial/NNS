# NNS Boost

Ensemble method for classification and regression. Each base learner
fits
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md) on
a sampled subset of predictors.

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
  seed = 123L,
  dist = NULL,
  folds = 5
)
```

## Arguments

- IVs.train:

  a matrix or data frame of variables of numeric or factor data types.

- DV.train:

  a numeric or factor vector with compatible dimensions to
  `(IVs.train)`.

- IVs.test:

  a matrix or data frame of variables of numeric or factor data types
  with compatible dimensions to `(IVs.train)`. If NULL, will use
  `(IVs.train)` as default. Columns are matched to `(IVs.train)` by name
  when the two share the same predictor names, and positionally when
  they share no names at all (as in `cbind(test.x_1, test.x_2)` against
  `cbind(x_1, x_2)`). A CSV UTF-8 byte-order mark in the first header
  (for example, `Male` becoming `X...Male`) is repaired before matching
  names. Other partial overlaps, including misspelled predictors, are
  ambiguous and return an error.

- type:

  With `NULL` (default), Boost selects classification for factor,
  character or logical responses. Use `"CLASS"` for other numeric
  labels.

- depth:

  options: (integer, NULL, "max"); `(depth = NULL)`(default) Specifies
  the `order` parameter in the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  routine, assigning a number of splits in the regressors, analogous to
  tree depth.

- learner.trials:

  integer; 100 (default) Sets the number of trials to obtain an accuracy
  `threshold` level. If the number of all possible feature combinations
  is less than selected value, the minimum of the two values will be
  used.

- epochs:

  integer; `2*length(DV.train)` (default) Number of repeated holdout
  re-evaluations of the learner-trial feature subsets that pass the
  accuracy threshold.

- CV.size:

  numeric \[0, 1\]; `NULL` (default) Sets the cross-validation size.
  Defaults to a random value between 0.2 and 0.33 for a random sampling
  of the training set.

- balance:

  logical; `FALSE` (default) Uses both up and down sampling to balance
  the classes. With `type = NULL`, balancing selects classification.

- ts.test:

  integer; NULL (default) Sets the length of the test set for
  time-series data; typically `2*h` parameter value from
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  or double known periods to forecast.

- threshold:

  numeric \[0, 1\]; `NULL` (default) Probability supplied to
  [LPM.VaR](https://OVVO-Financial.github.io/NNS/reference/LPM.VaR.md)
  over the learner-trial objective distribution to determine the
  objective cutoff for keeping feature combinations. Defaults to 0.80
  when `objective = "max"` and 0.20 when `objective = "min"`. It is not
  a literal objective-score cutoff.

- obj.fn:

  expression; `expression( sum((predicted - actual)^2) )` (default) Sum
  of squared errors is the default objective function. Any
  `expression(...)` using the specific terms `predicted` and `actual`
  can be used. Automatically selects an accuracy measure when
  `(type = "CLASS")`.

- objective:

  options: ("min", "max") `"min"` (default) Select whether to minimize
  or maximize the objective function `obj.fn`.

- extreme:

  logical; `FALSE` (default) Sets the
  [LPM.VaR](https://OVVO-Financial.github.io/NNS/reference/LPM.VaR.md)
  probability to 1 (0) for maximization (minimization) `objective`, i.e.
  the most extreme learner-trial objective value becomes the cutoff.
  Overrides `threshold`.

- features.only:

  logical; `FALSE` (default) Returns only the final feature loadings
  along with the final feature frequencies.

- feature.importance:

  logical; `TRUE` (default) Draws a two-panel diagnostic: the
  learner-trial objective distribution with its
  [LPM.VaR](https://OVVO-Financial.github.io/NNS/reference/LPM.VaR.md)
  cutoff, and the frequency of features used in the final estimate.

- pred.int:

  Prediction interval level strictly between 0 and 1 for regression or
  binary classification. Defaults to `NULL`.

- status:

  logical; `TRUE` (default) Prints status update message in console.

- seed:

  Integer seed for resampling and fitting. The default is `123L`. With
  `seed = NULL`, the fit uses the current RNG state and restores it on
  return.

- dist:

  options:(NULL, "NNS", "L1", "L2", "FACTOR") the method of distance
  calculation passed to delegated
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  and
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)
  calls. `dist = NULL` is the default and selects the native blended NNS
  distance; `dist = "NNS"` is an explicit alias for the default. `"L1"`
  and `"L2"` use Manhattan and Euclidean distance, respectively, on
  range-normalized coordinates. `"FACTOR"` uses the proportion of
  encoded predictor columns that differ.

- folds:

  integer; 5 (default) Number of cross-validation `folds` passed to the
  final
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)
  call.

## Value

With `features.only = FALSE`, returns a list containing the test
predictions `$results`, `$pred.int` (prediction intervals when requested
for regression or binary classification, otherwise `NULL`), the final
feature loadings `$feature.weights`, final feature frequencies
`$feature.frequency`, and (for classification) the class labels
`$class.levels`. Factor and character responses return one-based numeric
codes. Use `class.levels[results]` to recover the labels. Observed
factor labels follow their declared order. Character labels are sorted.
Numeric and logical classes retain their values. The objective function
receives the same values as the predictions.

## Note

- Like a logistic regression, the `(type = "CLASS")` setting is not
  necessary for target variable of two classes e.g. \[0, 1\]. The
  response variable base category should be 1 for classification
  problems.

- Incorporate any objective function from external packages (such as
  `Metrics::mape`) via
  `NNS.boost(..., obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")`

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
 a <- NNS.boost(iris[-test.set, 1:4], iris[-test.set, 5],
 IVs.test = iris[test.set, 1:4],
 epochs = 100, learner.trials = 100,
 type = "CLASS", depth = NULL, balance = TRUE)

 ## Test accuracy
 mean(a$results == as.numeric(iris[test.set, 5]))

 ## Recover the labels
 a$class.levels[a$results]
 } # }
```
