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

  At degree 1, `"discrete"` (default) uses empirical FSD and
  `"continuous"` uses the degree-one area-ratio (Omega/expectile)
  ordering. The latter does not imply ordinary FSD. SSD and TSD are
  unaffected by `type`.

- min_cluster:

  Integer. Peeling stops once at most `min_cluster` variables remain;
  those form the final cluster, and a final cluster smaller than
  `min_cluster` is merged into the previous one. SD-efficient sets
  peeled before that are kept whatever their size.

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
iteratively, peeling off the SD-efficient set at each step (whatever its
size) while more than `min_cluster` variables remain, or until no more
subsets can be extracted. Variables in each SD-efficient set form a
cluster. The remaining variables form the final cluster, which is merged
into the previous cluster if it has fewer than `min_cluster` members.

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
