# NNS Partition Map

Creates partitions based on partial-moment quadrant centroids. Numeric
orders use the compiled recursive partitioner. \`order = "max"\` returns
the maximum observable partition representation without passing an
invalid integer to C++.

## Usage

``` r
NNS.part(
  x,
  y,
  Voronoi = FALSE,
  type = NULL,
  order = NULL,
  obs.req = 8,
  min.obs.stop = TRUE,
  noise.reduction = "off"
)
```

## Arguments

- x:

  Numeric vector.

- y:

  Numeric vector of the same length as x.

- Voronoi:

  Logical; draw the partition map.

- type:

  NULL or "XONLY".

- order:

  NULL, a positive integer, or "max".

- obs.req:

  Nonnegative integer minimum-observation stopping control.

- min.obs.stop:

  Logical stopping control.

- noise.reduction:

  One of "mean", "median", "mode", "mode_class", or "off".

## Value

A list containing \`order\`, \`dt\`, and \`regression.points\`.
