# NNS Dependence

Returns the dependence and nonlinear correlation between two variables
based on higher order partial moment matrices measured by frequency or
area.

## Usage

``` r
NNS.dep(x, y = NULL, asym = FALSE, p.value = FALSE, print.map = FALSE)
```

## Arguments

- x:

  a numeric vector, matrix or data frame.

- y:

  `NULL` (default) or a numeric vector with compatible dimensions to
  `x`.

- asym:

  logical; `FALSE` (default) returns the symmetric measures. If `TRUE`,
  returns the asymmetric (directional) measures: the partition is built
  on `x` only and the result measures the dependence of `y` on `x` (`y`
  as a function of `x`); swap the arguments for the other direction.

- p.value:

  logical; `FALSE` (default) If `TRUE` (bivariate only), compares the
  results against 100 random permutations of `y`: the Dependence p-value
  is the upper-tail permutation p-value
  `(1 + #(null >= observed)) / (1 + 100)`, the Correlation p-value is
  the two-sided permutation p-value
  `min(1, 2 * min(1 + #(null <= observed), 1 + #(null >= observed)) / (1 + 100))`,
  and the 95 percent intervals of the permutation null are returned (and
  plotted when `print.map = TRUE`).

- print.map:

  logical; `FALSE` (default) Plots quadrant means, or p-value
  replicates.

## Value

Returns the bi-variate `"Correlation"` and `"Dependence"` or correlation
/ dependence matrix for matrix input.

## Details

Dependence is a finite-resolution summary of the recursive partitions of
[NNS.part](https://OVVO-Financial.github.io/NNS/reference/NNS.part.md),
which split each cell at the gravity of its observations
(`noise.reduction = "off"`), not at its mean. At each cell, degree-0
(frequency) and degree-1 (area) concordance around the cell means are
measured relative to a within-cell permutation reference and rescaled.
Their cell-size-weighted averages are evaluated across partition levels.
The reported Dependence is the strongest level, floored at zero;
Correlation is a signed summary at that selected level and is bounded by
Dependence in absolute value. With `asym = FALSE`, partitions on x, on
y, and on joint quadrants are considered. With `asym = TRUE`, only
partitions on x are considered, so the direction describes the
resolution supplied by the conditioning variable rather than a causal
effect.

The split depth, stopping rules, finite occupancy, and maximization
across levels affect the score. A zero score does not prove
independence, a noiseless nonlinear function need not score one, and
finite independent samples can produce positive scores. The permutation
calculation repeats the same selection rule, providing a reference for
the selected statistic rather than removing every finite-sample effect
from the score itself.

The degree-0 cell reference uses the hypergeometric distribution, with a
normal approximation for cell sizes above 400. The degree-1 reference
enumerates distinct arrangements only when at most 720 are required (all
arrangements for up to six distinct values); larger cells use an
approximation. These are properties of the cell calibration, separate
from the 100 whole-sample permutations used by `p.value = TRUE`.

Whole-sample shuffling assumes exchangeability of observations under the
null. It is not an automatic test for serially dependent or clustered
observations. The returned fields named `95% CIs` are central intervals
of the permutation null, not confidence intervals for a population
dependence parameter. Inclusive-tail counts with the plus-one correction
avoid zero Monte Carlo p-values; ties make the test conservative.

## Note

For asymmetrical `(asym = TRUE)` matrices, entry `[i, j]` is the
directional dependence (\[row variable i\] —\> \[column variable j\]),
i.e. `NNS.dep(x[, i], x[, j], asym = TRUE)`, the dependence of column
`j` on row variable `i`.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.dep(x, y)

## Correlation / Dependence Matrix
x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
B <- cbind(x, y, z)
NNS.dep(B)
} # }
```
