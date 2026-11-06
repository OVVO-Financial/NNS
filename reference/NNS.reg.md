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
  `(noise.reduction = "off")`. `(order = "max")` retains distinct
  predictor locations. A perfect training fit requires compatible
  responses at tied locations and interpolation at those locations;
  duplicate predictors with different responses cannot all be
  interpolated.

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
  for directional association weights. `(dim.red.method = "cor")` uses
  Spearman rank correlation for weights. `(dim.red.method = "all")`
  averages all methods for further feature engineering.
  `(dim.red.method = "equal")` uses unit weights. Alternatively, user
  can specify a numeric vector of coefficients. Named vectors must
  exactly match the encoded predictor columns and are matched by name;
  unnamed vectors follow encoded column order.

- tau:

  options("ts", NULL); `NULL`(default) To be used in conjunction with
  `(dim.red.method = "NNS.caus")` or `(dim.red.method = "all")`. If the
  regression is using time-series data, set `(tau = "ts")` to use the
  time-series lag convention.

- type:

  `NULL` (default). To perform a classification, set to
  `(type = "CLASS")`. Like a logistic regression, it is not necessary
  for target variable of two classes e.g. \[0, 1\]: a numeric response
  with exactly the values 0 and 1, or a factor, character or logical
  response, is classified automatically. Use `type = "REG"` to fit a
  numeric response as regression, including a two-valued response, with
  one or more predictors. `(type = "XONLY")` is the single-predictor
  regression alias and keeps the default x-only partition.

- point.est:

  a numeric or factor vector with compatible dimensions to `x`. Returns
  the fitted value `y.hat` for any value of `x`. Columns are matched to
  `x` by name when the two share the same predictor names, and
  positionally when they share no names at all (as in
  `cbind(test.x_1, test.x_2)` against `cbind(x_1, x_2)`). A CSV UTF-8
  byte-order mark in the first header (for example, `Male` becoming
  `X...Male`) is repaired before matching names. Other partial overlaps,
  including misspelled predictors, are ambiguous and return an error.

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
  interval with the estimate and adds `conf.int.pos` / `conf.int.neg` to
  `Fitted.xy`. The intervals are global: they are built from the type-8
  quantiles of the in-sample residuals (no per-segment standard errors
  are returned). Also applies the same level for the prediction
  intervals.

- threshold:

  numeric \[0, 1\]; `(threshold = 0)` (default) Sets the threshold for
  dimension reduction of independent variables when `(dim.red.method)`
  is not `NULL`. Predictors whose coefficient magnitude is below
  `threshold` are dropped; `NULL` or `NA` means 0. For classification, a
  vector with one cutoff per class (named by class, or in class order)
  applies each cutoff to that class's index, as
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)
  can report in `NNS.dim.red.threshold`.

- n.best:

  integer; `NULL` (default) Sets the number of nearest regression points
  to use in weighting for multivariate regression at
  `floor((1 - D) * sqrt(nrow(x)))` (at least 1, at most the number of
  regression points), where D averages the copula summaries with
  response y and -y so the smoothing choice is unchanged by response
  reflection. `(n.best = "all")` will select and weight all generated
  regression points. Analogous to `k` in a `k Nearest Neighbors`
  algorithm. Different values of `n.best` are tested using
  cross-validation in
  [NNS.stack](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md).
  `NA` selects the default, as `NULL` does.

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

- `"R2"` is squared Pearson correlation between fitted and observed
  responses (Racine–Hastie goodness of fit), or fitted classification
  accuracy;

- `"SE"` returns the overall standard error of the estimate between `y`
  and `y.hat`;

- `"Prediction.Accuracy"` returns the proportion of fitted class
  decisions equal to the observed class;

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

- `"R2"` is squared Pearson correlation between fitted and observed
  responses (Racine–Hastie goodness of fit), or fitted classification
  accuracy;

- `"equation"` returns the numerator of the synthetic X\* dimension
  reduction equation as a `data.frame` consisting of regressor and its
  coefficient. Denominator is the number of non-zero coefficients
  (`abs(coefficient) > 0`), returned in last row of `equation`
  `data.frame`.

- `"x.star"` returns the synthetic X\* as a vector;

- `"rhs.partitions"` returns the partition points for each regressor
  `x`;

- `"RPM"` provides the Regression Point Matrix, the points for each `x`
  used in the regression equation for the given order of partitions;

- `"Point.est"` returns the predicted value generated;

- `"pred.int"` lower and upper prediction intervals for the
  `"Point.est"` returned using the `"confidence.interval"` provided;

- `"Fitted.xy"` returns a `data.frame` of the regressors, `y`, `y.hat`,
  `NNS.ID` and `residuals` for full multivariate regression, or of the
  synthetic `x`, `y`, `y.hat`, `NNS.ID`, `gradient` and `residuals` for
  dimension reduction regression. `confidence.interval` adds
  `conf.int.neg` and `conf.int.pos`.

## Details

For regression, `R2` measures association between fitted and observed
responses: it is \\cor(y, \hat y)^2\\, with the package's finite
constant-vector convention. It is not \\1-SSE/SST\\; a biased or
incorrectly scaled prediction can have `R2 = 1`. Evaluate held-out
prediction errors for predictive accuracy.

The default partition order combines a finite-partition dependence score
with the average of the duplicated-predictor copula score under both
response orientations. This complexity heuristic is invariant to
reflecting the response; the underlying directional dependence and
copula definitions are unchanged. A supplied order bypasses this
heuristic.

Univariate regression linearly interpolates the selected regression
points and extends the first or last segment beyond their range (unless
a smoothing spline is requested). Endpoint slopes can be unstable with
sparse or noisy boundary points. Multivariate distance weighting is a
weighted average of regression-point responses, so it does not
automatically continue a trend outside the training support. Neither
construction supplies a general extrapolation guarantee. With the
default `n.best = NULL`, purely categorical predictors use exact
observed joint cells, and an unobserved joint combination uses its
nearest cells. An explicit `n.best` weights that many nearest cells, as
for numeric predictors. A declared factor level with no training rows is
accepted in `point.est` and is placed by the same distance rule. A
coordinate with exactly two observed values retains both terminal cells,
including a dummy column in a mixed design. Other finite-order numeric
coordinates retain interval cells; `order = "max"` uses exact joint
support.

Classification forms an indicator for each observed class, retains the
empirical class frequencies in each predictor cell, and interpolates or
distance-weights those probability vectors. Class labels never enter
numeric interpolation. `Point.probabilities`, `Fitted.probabilities`,
and `RPM.probabilities` expose the corresponding vectors. Tied class
decisions follow first appearance in the training response. Factor and
character responses retain numeric output codes with `class.levels`;
numeric responses retain their supplied values. For dimension reduction,
one signed NNS index is fitted per class indicator, and the resulting
probabilities are normalized across classes; `equation` is then a list
indexed by class. `confidence.interval` returns `Prediction.set`, the
smallest leading set of classes attaining the requested estimated
probability mass. These are descriptive probability sets, not calibrated
coverage guarantees; nominal classes have no numeric prediction
interval.

## Note

- Please ensure `point.est` is of compatible dimensions to `x`, error
  message will ensue if not compatible.

- Like a logistic regression, the `(type = "CLASS")` setting is not
  necessary for target variable of two classes e.g. \[0, 1\]. For other
  numeric responses set `type = "CLASS"`; class values can be any
  distinct numbers.

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
