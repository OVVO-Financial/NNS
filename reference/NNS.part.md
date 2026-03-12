# NNS Partition Map

Creates partitions based on partial moment quadrant centroids,
iteratively assigning identifications to observations based on those
quadrants (unsupervised partitional and hierarchical clustering method).
Basis for correlation, dependence
[NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md),
regression
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
routines.

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

  a numeric vector.

- y:

  a numeric vector with compatible dimensions to `x`.

- Voronoi:

  logical; `FALSE` (default) Displays a Voronoi type diagram using
  partial moment quadrants.

- type:

  `NULL` (default) Controls the partitioning basis. Set to
  `(type = "XONLY")` for X-axis based partitioning. Defaults to `NULL`
  for both X and Y-axis partitioning.

- order:

  integer; Number of partial moment quadrants to be generated.
  `(order = "max")` will institute a perfect fit.

- obs.req:

  integer; (8 default) Required observations per cluster where quadrants
  will not be further partitioned if observations are not greater than
  the entered value. Reduces minimum number of necessary observations in
  a quadrant to 1 when `(obs.req = 1)`.

- min.obs.stop:

  logical; `TRUE` (default) Stopping condition where quadrants will not
  be further partitioned if a single cluster contains less than the
  entered value of `obs.req`.

- noise.reduction:

  the method of determining regression points options for the dependent
  variable `y`: ("mean", "median", "mode", "off");
  `(noise.reduction = "mean")` uses means for partitions.
  `(noise.reduction = "median")` uses medians instead of means for
  partitions, while `(noise.reduction = "mode")` uses modes instead of
  means for partitions. Defaults to `(noise.reduction = "off")` where an
  overall central tendency measure is used, which is the default for the
  independent variable `x`.

## Value

Returns:

- `"dt"` a `data.table` of `x` and `y` observations with their partition
  assignment `"quadrant"` in the 3rd column and their prior partition
  assignment `"prior.quadrant"` in the 4th column.

- `"regression.points"` the `data.table` of regression points for that
  given `(order = ...)`.

- `"order"` the `order` of the final partition given `"min.obs.stop"`
  stopping condition.

## Note

`min.obs.stop = FALSE` will not generate regression points due to
unequal partitioning of quadrants from individual cluster observations.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.part(x, y)

## Data.table of observations and partitions
NNS.part(x, y, order = 1)$dt

## Regression points
NNS.part(x, y, order = 1)$regression.points

## Voronoi style plot
NNS.part(x, y, Voronoi = TRUE)

## Examine final counts by quadrant
DT <- NNS.part(x, y)$dt
DT[ , counts := .N, by = quadrant]
DT
} # }
```
