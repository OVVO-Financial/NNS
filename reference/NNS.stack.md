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

  a vector, matrix or data frame of variables of numeric or factor data
  types.

- DV.train:

  a numeric or factor vector with compatible dimensions to
  `(IVs.train)`.

- IVs.test:

  a vector, matrix or data frame of variables of numeric or factor data
  types with compatible dimensions to `(IVs.train)`. If NULL, will use
  `(IVs.train)` as default.

- type:

  `NULL` (default). To perform a classification of discrete integer
  classes from factor target variable `(DV.train)` with a base category
  of 1, set to `(type = "CLASS")`, else for continuous `(DV.train)` set
  to `(type = NULL)`. Like a logistic regression, this setting is not
  necessary for target variable of two classes e.g. \[0, 1\].

- obj.fn:

  expression; `expression(sum((predicted - actual)^2))` (default) Sum of
  squared errors is the default objective function. Any
  [`expression()`](https://rdrr.io/r/base/expression.html) using the
  specific terms `predicted` and `actual` can be used.

- objective:

  options: ("min", "max") `"min"` (default) Select whether to minimize
  or maximize the objective function `obj.fn`.

- optimize.threshold:

  logical; `TRUE` (default) Will optimize the probability threshold
  value for rounding in classification problems. If `FALSE`, returns
  0.5.

- dist:

  options:("L1", "L2", "DTW", "FACTOR") the method of distance
  calculation; Selects the distance calculation used. `dist = "L2"`
  (default) selects the Euclidean distance and `(dist = "L1")` selects
  the Manhattan distance; `(dist = "DTW")` selects the dynamic time
  warping distance; `(dist = "FACTOR")` uses a frequency.

- CV.size:

  numeric \[0, 1\]; `NULL` (default) Sets the cross-validation size if
  `(IVs.test = NULL)`. Defaults to a random value between 0.2 and 0.33
  for a random sampling of the training set.

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

  logical; `TRUE` (default) Uses dimension reduction output in `n.best`
  optimization, otherwise performs both analyses independently.

- dim.red.method:

  options: ("cor", "NNS.dep", "NNS.caus", "equal", "all") method for
  determining synthetic X\* coefficients. `(dim.red.method = "cor")`
  uses standard linear correlation for weights.
  `(dim.red.method = "NNS.dep")` (default) uses
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  for nonlinear dependence weights, while
  `(dim.red.method = "NNS.caus")` uses
  [NNS.caus](https://OVVO-Financial.github.io/NNS/reference/NNS.caus.md)
  for causal weights. `(dim.red.method = "all")` averages all methods
  for further feature engineering.

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

- seed:

  Optional integer random seed used for reproducible resampling, fold
  construction, and stochastic fitting steps. If \`NULL\`, the current
  random-number-generator state is used.

## Value

Returns a vector of fitted values for the dependent variable test set
for all models.

- `"NNS.reg.n.best"` returns the optimum `"n.best"` parameter for the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  multivariate regression. `"SSE.reg"` returns the SSE for the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  multivariate regression.

- `"OBJfn.reg"` returns the `obj.fn` for the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  regression.

- `"NNS.dim.red.threshold"` returns the optimum `"threshold"` from the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  dimension reduction regression.

- `"OBJfn.dim.red"` returns the `obj.fn` for the
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  dimension reduction regression.

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

## Note

- Incorporate any objective function from external packages (such as
  `Metrics::mape`) via
  `NNS.stack(..., obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")`

- Like a logistic regression, the `(type = "CLASS")` setting is not
  necessary for target variable of two classes e.g. \[0, 1\]. The
  response variable base category should be 1 for multiple class
  problems.

- Missing data should be handled prior as well using
  [na.omit](https://rdrr.io/r/stats/na.fail.html) or
  [complete.cases](https://rdrr.io/r/stats/complete.cases.html) on the
  full dataset.

If error received:

`"Error in is.data.frame(x) : object 'RP' not found"`

reduce the `CV.size`.

## References

Viole, F. (2016) "Classification Using NNS Clustering Analysis"
[doi:10.2139/ssrn.2864711](https://doi.org/10.2139/ssrn.2864711)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
 ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
 if (FALSE) { # \dontrun{
 NNS.stack(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4], type = "CLASS", 
 balance = TRUE)

 ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
 NNS.stack(iris[ , 1:4], iris[ , 5], type = "CLASS")
 } # }
```
