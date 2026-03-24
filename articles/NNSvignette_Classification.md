# Getting Started with NNS: Classification

``` r
library(NNS)
library(data.table)
require(knitr)
require(rgl)
```

## Classification

**`NNS.reg`** is a very robust regression technique capable of nonlinear
regressions of continuous variables and classification tasks in machine
learning problems.

We have extended the **`NNS.reg`** applications per the use of an
ensemble method of classification in **`NNS.boost`**. In short,
**`NNS.reg`** is the base learner instead of trees.

***One major advantage `NNS.boost` has over tree based methods is the
ability to seamlessly extrapolate beyond the current range of
observations.***

### Splits vs. Partitions

Popular boosting algorithms take a series of weak learning decision tree
models, and aggregate their outputs. `NNS` is also a decision tree of
sorts, by partitioning each regressor with respect to the dependent
variable. We can directly control the number of “splits” with the
**`NNS.reg(..., order = , ...)`** parameter.

#### NNS Partitions

We can see how `NNS` partitions each regressor by calling the
`$rhs.partitions` output. You will notice that each partition is not an
equal interval, nor of equal length, which differentiates `NNS` from
other bandwidth or tree-based techniques.

Higher dependence between a regressor and the dependent variable will
allow for a larger number of partitions. This is determined internally
with the **`NNS.dep`** measure.

``` r
NNS.reg(iris[,1:4], iris[,5], residual.plot = FALSE, ncores = 1)$rhs.partitions
```

    ##           V1       V2       V3        V4
    ##        <num>    <num>    <num>     <num>
    ##  1: 4.300000 2.000000 1.000000 0.1000000
    ##  2: 4.381250 2.829734 1.050000 0.2000000
    ##  3: 4.577396 3.358016 1.200000 0.3000000
    ##  4: 4.700000 4.400000 1.300000 0.4083333
    ##  5: 4.800000       NA 1.400000 1.0630397
    ##  6: 4.900000       NA 1.500000 1.3742424
    ##  7: 5.000000       NA 1.600000 1.7789810
    ##  8: 5.100000       NA 1.700000 2.2285885
    ##  9: 5.205000       NA 1.900000 2.5000000
    ## 10: 5.400000       NA 3.416305        NA
    ## 11: 5.500000       NA 3.834865        NA
    ## 12: 5.600000       NA 4.000000        NA
    ## 13: 5.700000       NA 4.184722        NA
    ## 14: 5.800000       NA 4.400000        NA
    ## 15: 5.900000       NA 4.500000        NA
    ## 16: 6.000000       NA 4.670803        NA
    ## 17: 6.100000       NA 4.863889        NA
    ## 18: 6.200000       NA 5.000000        NA
    ## 19: 6.300000       NA 5.100000        NA
    ## 20: 6.400000       NA 5.200000        NA
    ## 21: 6.500000       NA 5.337500        NA
    ## 22: 6.600000       NA 5.500000        NA
    ## 23: 6.700000       NA 5.617708        NA
    ## 24: 6.800000       NA 5.849554        NA
    ## 25: 6.900000       NA 6.336875        NA
    ## 26: 7.050000       NA 6.900000        NA
    ## 27: 7.224375       NA       NA        NA
    ## 28: 7.687079       NA       NA        NA
    ## 29: 7.900000       NA       NA        NA
    ##           V1       V2       V3        V4
    ##        <num>    <num>    <num>     <num>

## `NNS.boost()`

Through resampling of the training set and letting each iterated set of
data speak for themselves (while paying extra attention to the residuals
throughout), we can test various regressor combinations in these dynamic
decision trees…only keeping those combinations that add predictive
value. From there we simply aggregate the predictions.

**`NNS.boost`** will automatically search for an accuracy `threshold`
from the training set, reporting iterations remaining and level obtained
in the console. A plot of the frequency of the learning accuracy on the
training set is also provided.

Once a `threshold` is obtained, **`NNS.boost`** will test various
feature combinations against different splits of the training set and
report back the frequency of each regressor used in the final estimate.

Let’s have a look and see how it works. We use 140 random `iris`
observations as our training set with the 10 holdout observations as our
test set. For brevity, we set
`epochs = 10, learner.trials = 10, folds = 1`.

**NOTE: Base category of response variable should be 1, not 0 for
classification problems when using `NNS.boost(..., type = "CLASS")`**.

``` r
test.set = 141:150
 
a = NNS.boost(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              epochs = 10, learner.trials = 10, 
              status = FALSE, balance = TRUE,
              type = "CLASS", folds = 5)

a
$results
 [1] 3 3 3 3 3 3 3 3 3 3

$pred.int
NULL

$feature.weights
 Petal.Width Petal.Length Sepal.Length 
   0.4285714    0.4285714    0.1428571 

$feature.frequency
 Petal.Width Petal.Length Sepal.Length 
           3            3            1 
   
mean( a$results == as.numeric(iris[test.set, 5]) )
[1] 1
```

A perfect classification, using the features weighted per the output
above.

## Cross-Validation Classification Using `NNS.stack()`

The
**[`NNS.stack()`](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)**
routine cross-validates for a given objective function the `n.best`
parameter in the multivariate **`NNS.reg`** function as well as the
`threshold` parameter in the dimension reduction **`NNS.reg`** version.
**`NNS.stack`** can be used for classification via
**`NNS.stack(..., type = "CLASS", ...)`**.

For brevity, we set `folds = 1`.

**NOTE: Base category of response variable should be 1, not 0 for
classification problems when using `NNS.stack(..., type = "CLASS")`**.

``` r
b = NNS.stack(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              type = "CLASS", balance = TRUE,
              ncores = 1, folds = 5)

b
```

``` r
$OBJfn.reg
[1] 0.955787

$NNS.reg.n.best
[1] 1

$probability.threshold
[1] 0.6429167

$OBJfn.dim.red
[1] 0.955787

$NNS.dim.red.threshold
[1] 0.925

$reg
 [1] 3 3 3 3 3 3 3 3 3 3

$reg.pred.int
NULL

$dim.red
 [1] 3 3 3 3 3 3 3 3 3 3

$dim.red.pred.int
NULL

$stack
 [1] 3 3 3 3 3 3 3 3 3 3

$pred.int
NULL
```

``` r
mean( b$stack == as.numeric(iris[test.set, 5]) )
```

``` r
[1] 1
```

### Brief Notes on Other Parameters

- `depth = "max"` will force all observations to be their own partition,
  forcing a perfect fit of the multivariate regression. In essence, this
  is the basis for a `kNN` nearest neighbor type of classification.

- `n.best = 1` will use the single nearest neighbor. When coupled with
  `depth = "max"`, `NNS` will emulate a `kNN = 1` but as the dimensions
  increase the results diverge demonstrating `NNS` is less sensitive to
  the curse of dimensionality than `kNN`.

- `extreme` will use the maximum or minimum `threshold` obtained, and
  may result in errors if that threshold cannot be eclipsed by
  subsequent iterations.

## References

If the user is so motivated, detailed arguments further examples are
provided within the following:

- [Nonlinear Nonparametric Statistics: Using Partial
  Moments](https://ovvo-financial.github.io/NNS/book/)

- [Deriving Nonlinear Correlation Coefficients from Partial
  Moments](https://doi.org/10.2139/ssrn.2148522)

- [Nonparametric Regression Using
  Clusters](https://doi.org/10.1007/s10614-017-9713-5)

- [Clustering and Curve Fitting by Line
  Segments](https://doi.org/10.2139/ssrn.2861339)

- [Classification Using NNS Clustering
  Analysis](https://doi.org/10.2139/ssrn.2864711)
