# NNS Regression

Generates a nonlinear regression based on partial moment quadrant means.

## Usage

``` r
NNS.reg(
  x,
  y,
  factor.2.dummy = TRUE,
  order = NULL,
  dim.red.method = NULL,
  tau = NULL,
  type = NULL,
  point.est = NULL,
  location = "top",
  return.values = TRUE,
  plot = TRUE,
  plot.regions = FALSE,
  residual.plot = TRUE,
  confidence.interval = NULL,
  threshold = 0,
  n.best = NULL,
  smooth = FALSE,
  noise.reduction = "off",
  dist = NULL,
  ncores = NULL,
  point.only = FALSE,
  multivariate.call = FALSE
)
```

## Arguments

- x:

  a vector, matrix or data frame of variables of numeric or factor data
  types.

- y:

  a numeric or factor vector with compatible dimensions to `x`.

- factor.2.dummy:

  logical; `TRUE` (default) Automatically augments variable matrix with
  numerical dummy variables based on the levels of factors.

- order:

  integer; Controls the number of partial moment quadrant means. Users
  are encouraged to try different `(order = ...)` integer settings with
  `(noise.reduction = "off")`. `(order = "max")` will force a limit
  condition perfect fit.

- dim.red.method:

  options: ("cor", "NNS.dep", "NNS.caus", "all", "equal",
  `numeric vector`, NULL) method for determining synthetic X\*
  coefficients (per Dana and Dawes (2004)). Selection of a method
  automatically engages the dimension reduction regression. The default
  is `NULL` for full multivariate regression.
  `(dim.red.method = "NNS.dep")` uses
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  for nonlinear dependence weights, while
  `(dim.red.method = "NNS.caus")` uses
  [NNS.caus](https://OVVO-Financial.github.io/NNS/reference/NNS.caus.md)
  for causal weights. `(dim.red.method = "cor")` uses standard linear
  correlation for weights. `(dim.red.method = "all")` averages all
  methods for further feature engineering. `(dim.red.method = "equal")`
  uses unit weights. Alternatively, user can specify a numeric vector of
  coefficients.

- tau:

  options("ts", NULL); `NULL`(default) To be used in conjunction with
  `(dim.red.method = "NNS.caus")` or `(dim.red.method = "all")`. If the
  regression is using time-series data, set `(tau = "ts")` for more
  accurate causal analysis.

- type:

  `NULL` (default). To perform a classification, set to
  `(type = "CLASS")`. Like a logistic regression, it is not necessary
  for target variable of two classes e.g. \[0, 1\].

- point.est:

  a numeric or factor vector with compatible dimensions to `x`. Returns
  the fitted value `y.hat` for any value of `x`.

- location:

  Sets the legend location within the plot, per the `x` and `y`
  co-ordinates used in base graphics
  [legend](https://rdrr.io/r/graphics/legend.html).

- return.values:

  logical; `TRUE` (default), set to `FALSE` in order to only display a
  regression plot and call values as needed.

- plot:

  logical; `TRUE` (default) To plot regression.

- plot.regions:

  logical; `FALSE` (default). Generates 3d regions associated with each
  regression point for multivariate regressions. Note, adds significant
  time to routine.

- residual.plot:

  logical; `TRUE` (default) To plot `y.hat` and `Y`.

- confidence.interval:

  numeric \[0, 1\]; `NULL` (default) Plots the associated confidence
  interval with the estimate and reports the standard error for each
  individual segment. Also applies the same level for the prediction
  intervals.

- threshold:

  numeric \[0, 1\]; `(threshold = 0)` (default) Sets the threshold for
  dimension reduction of independent variables when `(dim.red.method)`
  is not `NULL`.

- n.best:

  integer; `NULL` (default) Sets the number of nearest regression points
  to use in weighting for multivariate regression at
  `sqrt(# of regressors)`. `(n.best = "all")` will select and weight all
  generated regression points. Analogous to `k` in a
  `k Nearest Neighbors` algorithm. Different values of `n.best` are
  tested using cross-validation in
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md).

- smooth:

  logical; `FALSE` (default) Applies a smoothing spline instead of local
  linear fit to regression points.

- noise.reduction:

  the method of determining regression points options: ("mean",
  "median", "mode", "off"); In low signal:noise
  situations,`(noise.reduction = "mean")` uses means for
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  restricted partitions, `(noise.reduction = "median")` uses medians
  instead of means for
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  restricted partitions, while `(noise.reduction = "mode")` uses modes
  instead of means for
  [NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
  restricted partitions. `(noise.reduction = "off")` uses an overall
  central tendency measure for partitions.

- dist:

  options:(NULL, "NNS", "L1", "L2", "FACTOR") the method of distance
  calculation. `dist = NULL` is the default and selects the native
  blended NNS distance, \\sum(abs(z) + z^2)\\, over range-normalized
  coordinates. `dist = "NNS"` is an explicit alias for the default.
  `dist = "L2"` selects Euclidean distance; `dist = "L1"` selects
  Manhattan distance; `dist = "FACTOR"` uses a frequency.

- ncores:

  integer; value specifying the number of cores to be used in the
  parallelized procedure. If NULL (default), the number of cores to be
  used is equal to the number of cores of the machine - 1.

- point.only:

  Internal argument for abbreviated output.

- multivariate.call:

  Internal argument for multivariate regressions.

## Value

UNIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:

- `"R2"` provides the goodness of fit;

- `"SE"` returns the overall standard error of the estimate between `y`
  and `y.hat`;

- `"Prediction.Accuracy"` returns the correct rounded `"Point.est"` used
  in classifications versus the categorical `y`;

- `"derivative"` for the coefficient of the `x` and its applicable
  range;

- `"Point.est"` for the predicted value generated;

- `"pred.int"` lower and upper prediction intervals for the
  `"Point.est"` returned using the `"confidence.interval"` provided;

- `"regression.points"` provides the points used in the regression
  equation for the given order of partitions;

- `"Fitted.xy"` returns a `data.frame` of `x`, `y`, `y.hat`, `resid`,
  `NNS.ID`, `gradient`;

MULTIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:

- `"R2"` provides the goodness of fit;

- `"equation"` returns the numerator of the synthetic X\* dimension
  reduction equation as a `data.frame` consisting of regressor and its
  coefficient. Denominator is simply the length of all coefficients \>
  0, returned in last row of `equation` `data.frame`.

- `"x.star"` returns the synthetic X\* as a vector;

- `"rhs.partitions"` returns the partition points for each regressor
  `x`;

- `"RPM"` provides the Regression Point Matrix, the points for each `x`
  used in the regression equation for the given order of partitions;

- `"Point.est"` returns the predicted value generated;

- `"pred.int"` lower and upper prediction intervals for the
  `"Point.est"` returned using the `"confidence.interval"` provided;

- `"Fitted.xy"` returns a `data.frame` of `x`,`y`, `y.hat`, `gradient`,
  and `NNS.ID`.

## Note

- Please ensure `point.est` is of compatible dimensions to `x`, error
  message will ensue if not compatible.

- Like a logistic regression, the `(type = "CLASS")` setting is not
  necessary for target variable of two classes e.g. \[0, 1\]. The
  response variable base category should be 1 for classification
  problems.

- For low signal:noise instances, increasing the dimension may yield
  better results using `NNS.stack(cbind(x,x), y, method = 1, ...)`.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Vinod, H. and Viole, F. (2017) "Nonparametric Regression Using Clusters"
[doi:10.1007/s10614-017-9713-5](https://doi.org/10.1007/s10614-017-9713-5)

Vinod, H. and Viole, F. (2018) "Clustering and Curve Fitting by Line
Segments"
[doi:10.20944/preprints201801.0090.v1](https://doi.org/10.20944/preprints201801.0090.v1)

Viole, F. (2020) "Partitional Estimation Using Partial Moments"
[doi:10.2139/ssrn.3592491](https://doi.org/10.2139/ssrn.3592491)

Dana, J., and Dawes, R. M. (2004). The Superiority of Simple
Alternatives to Regression for Social Science Predictions. Journal of
Educational and Behavioral Statistics, 29(3), 317–331.

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.reg(x, y)

## Manual {order} selection
NNS.reg(x, y, order = 2)

## Maximum {order} selection
NNS.reg(x, y, order = "max")

## x-only paritioning (Univariate only)
NNS.reg(x, y, type = "XONLY")

## For Multiple Regression:
x <- cbind(rnorm(100), rnorm(100), rnorm(100)) ; y <- rnorm(100)
NNS.reg(x, y, point.est = c(.25, .5, .75))

## For Multiple Regression based on Synthetic X* (Dimension Reduction):
x <- cbind(rnorm(100), rnorm(100), rnorm(100)) ; y <- rnorm(100)
NNS.reg(x, y, point.est = c(.25, .5, .75), dim.red.method = "cor", ncores = 1)

## IRIS dataset examples:
# Dimension Reduction:
NNS.reg(iris[,1:4], iris[,5], dim.red.method = "cor", order = 5, ncores = 1)

# Dimension Reduction using causal weights:
NNS.reg(iris[,1:4], iris[,5], dim.red.method = "NNS.caus", order = 5, ncores = 1)

# Multiple Regression:
NNS.reg(iris[,1:4], iris[,5], order = 2, noise.reduction = "off")

# Classification:
NNS.reg(iris[,1:4], iris[,5], point.est = iris[1:10, 1:4], type = "CLASS")$Point.est

## To call fitted values:
x <- rnorm(100) ; y <- rnorm(100)
NNS.reg(x, y)$Fitted

## To call partial derivative (univariate regression only):
NNS.reg(x, y)$derivative
} # }
```
