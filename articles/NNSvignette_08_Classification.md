# Getting Started with NNS: Classification

``` r

library(NNS)
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

For numeric regression, extrapolation follows the fitted NNS
construction and needs boundary validation. Nominal classification
predicts among the observed classes; numeric class codes do not define
an extrapolation scale.

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

    ## [[1]]
    ##  [1] 4.300000 4.612500 5.012000 5.388889 5.704762 6.057895 6.380952 6.766667
    ##  [9] 7.475000 7.900000
    ## 
    ## [[2]]
    ## [1] 2.000000 2.493939 2.924000 3.245238 3.752000 4.400000
    ## 
    ## [[3]]
    ##  [1] 1.000000 1.236364 1.400000 1.500000 1.636364 1.900000 3.000000 3.300000
    ##  [9] 3.575000 3.983333 4.300000 4.581250 4.900000 5.150000 5.536364 5.892308
    ## [17] 6.600000 6.900000
    ## 
    ## [[4]]
    ##  [1] 0.10 0.20 0.30 0.40 0.50 0.60 1.00 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80
    ## [16] 1.90 2.00 2.10 2.20 2.30 2.45 2.50

## `NNS.boost()`

**`NNS.boost`** screens feature subsets. Each learner trial fits
**`NNS.reg`** on a random subset of the regressors and scores it on a
held-out split. Subsets that clear a learned objective threshold are
re-tested on new holdouts in the epoch stage, and each feature is
counted every time a subset containing it passes again. The features are
then replicated in proportion to those counts and a single
**`NNS.stack`** (method 1) is fitted on that design. It is not an
ensemble of voting or averaged learners, and it does not reweight
residuals.

**`NNS.boost`** will automatically search for an accuracy `threshold`
from the training set, reporting iterations remaining and level obtained
in the console. A plot of the frequency of the learning accuracy on the
training set is also provided.

Once a `threshold` is obtained, **`NNS.boost`** will test various
feature combinations against different splits of the training set and
report back the frequency of each regressor used in the final estimate.

Let’s have a look and see how it works. We hold out four observations
from each iris species, using 138 observations for training and 12 for
testing. The split is stratified and reproducible; the last ten rows
alone would contain only virginica. For brevity, we set
`epochs = 10, learner.trials = 10, folds = 1`.

**NOTE:** class labels can be factor levels or any distinct numbers.
Numeric and logical responses are returned on their own values; factor
and character responses are returned as codes, which `class.levels` maps
back to the labels.

``` r

set.seed(123)
test.set <- sort(unlist(lapply(split(seq_len(nrow(iris)), iris$Species),
                               function(rows) sample(rows, 4))))
 
a = NNS.boost(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              epochs = 10, learner.trials = 10, 
              status = FALSE, balance = TRUE,
              type = "CLASS")

a
mean(a$results == as.numeric(iris[test.set, 5]))
```

The reported accuracy is for this held-out split. Repeated stratified
splits are needed to assess variability; a small demonstration does not
establish general superiority or robustness to increasing dimension.

## Cross-Validation Classification Using `NNS.stack()`

The
**[`NNS.stack()`](https://OVVO-Financial.github.io/NNS/reference/NNS.stack.md)**
routine cross-validates, for a given objective function, the neighbour
count `n.best` of the multivariate **`NNS.reg`** fit and the number of
predictors in the dimension reduction **`NNS.reg`** version, then blends
the two. For classification (**`NNS.stack(..., type = "CLASS", ...)`**)
it tunes and blends class-probability vectors. The `OBJfn` values it
returns are the best scores of that search on the out-of-fold
predictions, so they overstate held-out performance; the held-out test
set below is the assessment.

For brevity, we set `folds = 1`.

**NOTE:** class labels can be factor levels or any distinct numbers.
Numeric and logical responses are returned on their own values; factor
and character responses are returned as codes, which `class.levels` maps
back to the labels.

``` r

b = NNS.stack(IVs.train = iris[-test.set, 1:4], 
              DV.train = iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              type = "CLASS", balance = TRUE,
              ncores = 1, folds = 5)

b
```

``` r

mean(b$stack == as.numeric(iris[test.set, 5]))
```

The returned `class.levels` gives the label mapping for encoded class
results. Compare original labels after decoding when a classifier’s API
returns labels.

### Brief Notes on Other Parameters

- `depth = "max"` will force all observations to be their own partition,
  forcing a perfect fit of the multivariate regression. In essence, this
  is the basis for a `kNN` nearest neighbor type of classification.

- `n.best = 1` will use the single nearest neighbor. When coupled with
  `depth = "max"`, `NNS` will emulate a `kNN = 1` for the chosen
  partition and distance configuration. Performance as dimension
  increases requires a matched comparison of sample size, features,
  tuning, distance metric, and held-out loss; this example does not
  establish a general dimensionality advantage.

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
