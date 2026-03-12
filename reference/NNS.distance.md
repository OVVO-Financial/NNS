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

  Vector to generate distances from.

- k:

  `n.best` from
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)

- class:

  if classification problem.

## Value

Returns sum of weighted distances.
