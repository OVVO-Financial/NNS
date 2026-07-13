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
  seed = 123L
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

  integer; `2*length(DV.train)` (default) Total number of feature
  combinations to run.

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

  numeric; `NULL` (default) Sets the `obj.fn` threshold to keep feature
  combinations.

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

  logical; `FALSE` (default) Uses the maximum (minimum) `threshold`
  obtained from the `learner.trials`, rather than the upper (lower)
  quintile level for maximization (minimization) `objective`.

- features.only:

  logical; `FALSE` (default) Returns only the final feature loadings
  along with the final feature frequencies.

- feature.importance:

  logical; `TRUE` (default) Plots the frequency of features used in the
  final estimate.

- pred.int:

  numeric \[0,1\]; `NULL` (default) Returns the associated prediction
  intervals for the final estimate.

- status:

  logical; `TRUE` (default) Prints status update message in console.

- seed:

  Optional integer random seed used for reproducible resampling, fold
  construction, and stochastic fitting steps. If \`NULL\`, the current
  random-number-generator state is used.

## Value

Returns a vector of fitted values for the dependent variable test set
`$results`, prediction intervals `$pred.int`, and the final feature
loadings `$feature.weights`, along with final feature frequencies
`$feature.frequency`.

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
 } # }
```
