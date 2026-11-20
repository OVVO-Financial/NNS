# NNS Distance

Internal kernel function for NNS multivariate regression
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
parallel instances.

## Usage

``` r
NNS.distance(rpm, dist.estimate, k = "all", class = NULL)
```

## Arguments

- rpm:

  REGRESSION.POINT.MATRIX from
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)

- dist.estimate:

  Query coordinates. Matching names select RPM features. If no names
  match, coordinates use the first numeric features in order. A partial
  match or a repeated matching name gives an error. Each RPM row must
  share at least one finite coordinate with the query.

- k:

  `n.best` from
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)

- class:

  if classification problem.

## Value

A single numeric prediction at `dist.estimate`, in the units of `y.hat`
for regression or a numeric class code for classification. Equal
distances are ordered by RPM row.
