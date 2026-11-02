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
  the smaller of the two one-sided tail frequencies, and the 95 percent
  intervals of the permutation null are returned (and plotted when
  `print.map = TRUE`).

- print.map:

  logical; `FALSE` (default) Plots quadrant means, or p-value
  replicates.

## Value

Returns the bi-variate `"Correlation"` and `"Dependence"` or correlation
/ dependence matrix for matrix input.

## Details

Each partition cell's concordance (degree 0 and degree 1) is measured
against its expectation under the within-cell permutation null and
rescaled so that perfect dependence stays at 1. Independent variables
therefore score near 0 at any sample size, instead of a floor that
shrank only like `n^(-1/4)` (about 0.32 at n = 1000). Cells below chance
can pull the average down; Dependence is floored at 0 after averaging,
and `|Correlation|` never exceeds the Dependence of the same direction.

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
