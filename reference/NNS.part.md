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

- `order`: the order actually completed. It can be smaller than the
  requested `order` when `obs.req` / `min.obs.stop` stop the recursion
  (e.g. `order = 3` returns 2 for `n = 200` with the default
  `obs.req = 8`).

- `dt`: the observations with their final `quadrant` and
  `prior.quadrant` labels.

- `regression.points`: for a numeric `order`, one (x, y) centroid per
  `prior.quadrant`, i.e. per cell of the previous order: `order = 1`
  gives 4 final quadrants but a single regression point (`"q"`),
  `order = 2` gives up to 16 quadrants and 4 points.

## Details

For a numeric `order`, each cell is split at a centre per coordinate,
set by `noise.reduction`: `"mean"` and `"median"` use the mean and
median, `"mode"` uses the continuous `NNS.mode(multi = FALSE)`,
`"mode_class"` uses
[NNS.gravity](https://OVVO-Financial.github.io/NNS/reference/NNS.gravity.md)
for x and `NNS.mode(discrete = TRUE)` for y, and `"off"` uses
[NNS.gravity](https://OVVO-Financial.github.io/NNS/reference/NNS.gravity.md).
Regression points are the same centres of each `prior.quadrant`.
`order = "max"` makes no splits: every observation is its own regression
point, or with `type = "XONLY"` each distinct x gets the y centre of its
observations. Quadrant labels append one digit per order. With split
centres `cx` and `cy`, the digit for the default joint partition
(`type = NULL`) is `1 + (x <= cx) + 2 * (y <= cy)`: 1 = upper-right, 2 =
upper-left, 3 = lower-right, 4 = lower-left. For `type = "XONLY"` the
digit is 1 for `x <= cx` and 2 for `x > cx`.
