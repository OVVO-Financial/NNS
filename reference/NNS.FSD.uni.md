# NNS FSD Test uni-directional

Uni-directional test of first degree stochastic dominance using lower
partial moments used in SD Efficient Set routine.

## Usage

``` r
NNS.FSD.uni(x, y, type = "discrete")
```

## Arguments

- x:

  a numeric vector.

- y:

  a numeric vector.

- type:

  `"discrete"` (default) uses empirical FSD. At degree 1, `"continuous"`
  retains the degree-one area-ratio (Omega/expectile) ordering, which
  need not imply FSD or preference for every increasing utility. SSD and
  TSD are unaffected by `type`.

## Value

Returns (1) if X dominates Y in the selected order (ordinary FSD or
Omega/expectile order), else (0).

## Note

The curves are compared on the merged sample points and, for TSD and
continuous area order, also at their closed-form extrema between those
points. Differences smaller than `64 * .Machine$double.eps * S^d` (`S` =
largest absolute value in `x` and `y`, `d` = 1 for SSD and the mean
condition, 2 for TSD; a ratio analogue for continuous area order) are
treated as ties, so that rounding noise does not decide dominance.
Discrete FSD compares exact empirical CDF values.

## References

Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the
Computation of the SD Efficient Set." Journal of Mathematical Finance,
6, 105-126.
[doi:10.4236/jmf.2016.61012](https://doi.org/10.4236/jmf.2016.61012)

Viole, F. (2017) "A Note on Stochastic Dominance."
[doi:10.2139/ssrn.3002675](https://doi.org/10.2139/ssrn.3002675)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.FSD.uni(x, y)
} # }
```
