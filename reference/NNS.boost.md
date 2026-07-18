# NNS Boost

Ensemble method for classification using the NNS multivariate regression
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md) as
the base learner instead of trees.

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
  `(IVs.train)` as default.

- type:

  `NULL` (default). To perform a classification of discrete integer
  classes from factor target variable `(DV.train)` with a base category
  of 1, set to `(type = "CLASS")`, else for continuous `(DV.train)` set
  to `(type = NULL)`.

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
  the classes. `type="CLASS"` required.

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

  options: ("min", "max") `"max"` (default) Select whether to minimize
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

  numeric \[0,1\]; `NULL` (default) Returns the associated prediction
  intervals for the final estimate.

- status:

  logical; `TRUE` (default) Prints status update message in console.

- seed:

  Optional integer random seed used for reproducible resampling, fold
  construction, and stochastic fitting steps. If \`NULL\`, the current
  random-number-generator state is used.

- dist:

  options:(NULL, "NNS", "L1", "L2", "FACTOR") the method of distance
  calculation passed to delegated
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  and
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)
  calls. `dist = NULL` is the default and selects the native blended NNS
  distance; `dist = "NNS"` is an explicit alias for the default.

- folds:

  integer; 5 (default) Number of cross-validation `folds` passed to the
  final
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)
  call.

## Value

Returns a vector of fitted values for the dependent variable test set
`$results`, prediction intervals `$pred.int`, the final feature loadings
`$feature.weights`, final feature frequencies `$feature.frequency`, and
(for classification) the class labels `$class.levels`. Classification
results are numeric: a factor or character `DV.train` yields integer
class codes with a base category of 1 (label recoverable as
`class.levels[results]`), and a numeric `DV.train` yields its original
numeric class values.

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
 ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
 if (FALSE) { # \dontrun{
 a <- NNS.boost(iris[1:140, 1:4], iris[1:140, 5],
 IVs.test = iris[141:150, 1:4],
 epochs = 100, learner.trials = 100,
 type = "CLASS", depth = NULL, balance = TRUE)

 ## Test accuracy
 mean(a$results == as.numeric(iris[141:150, 5]))

 ## Recover the labels
 a$class.levels[a$results]
 } # }
```
