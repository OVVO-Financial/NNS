# NNS SD-based Clustering

Clusters a set of variables by iteratively extracting Stochastic
Dominance (SD)-efficient sets, subject to a minimum cluster size.

## Usage

``` r
NNS.SD.cluster(
  data,
  degree = 1,
  type = "discrete",
  min_cluster = 1,
  dendrogram = FALSE
)
```

## Arguments

- data:

  A numeric matrix or data frame of variables to be clustered.

- degree:

  Numeric options: (1, 2, 3). Degree of stochastic dominance test.

- type:

  Character, either `"discrete"` (default) or `"continuous"`; specifies
  the type of CDF.

- min_cluster:

  Integer. The minimum number of elements required for a valid cluster.

- dendrogram:

  Logical; `FALSE` (default). If `TRUE`, a dendrogram is produced based
  on a simple "distance" measure between clusters.

## Value

A list with the following components:

- `Clusters`: A named list of cluster memberships where each element is
  the set of variable names belonging to that cluster.

- `Dendrogram` (optional): If `dendrogram = TRUE`, an `hclust` object is
  also returned.

## Details

The function applies
[`NNS.SD.efficient.set`](https://OVVO-Financial.github.io/NNS/reference/NNS.SD.efficient.set.md)
iteratively, peeling off the SD-efficient set at each step if it meets
or exceeds `min_cluster` in size, until no more subsets can be extracted
or all variables are exhausted. Variables in each SD-efficient set form
a cluster, with any remaining variables aggregated into the final
cluster if it meets the `min_cluster` threshold.

## References

Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the
Computation of the SD Efficient Set." Journal of Mathematical Finance,
6, 105-126.
[doi:10.4236/jmf.2016.61012](https://doi.org/10.4236/jmf.2016.61012) .

Viole, F. (2017) "A Note on Stochastic Dominance."
[doi:10.2139/ssrn.3002675](https://doi.org/10.2139/ssrn.3002675)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)
A <- cbind(x, y, z)

# Perform SD-based clustering (degree 1), requiring at least 2 elements per cluster
results <- NNS.SD.cluster(data = A, degree = 1, min_cluster = 2)
print(results$Clusters)

# Produce a dendrogram as well
results_with_dendro <- NNS.SD.cluster(data = A, degree = 1, min_cluster = 2, dendrogram = TRUE)
} # }
```
