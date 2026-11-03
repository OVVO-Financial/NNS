# NNS FSD Test

Compares empirical first-order stochastic dominance or degree-one
area-ratio (Omega/expectile) dominance using partial moments.

## Usage

``` r
NNS.FSD(x, y, type = "discrete", plot = TRUE)
```

## Arguments

- x:

  a numeric vector.

- y:

  a numeric vector.

- type:

  `"discrete"` (default) compares empirical CDFs for ordinary
  first-order stochastic dominance. The legacy option `"continuous"`
  compares degree-one ratios \\L_1/(L_1+U_1)\\, equivalently
  Omega/expectile order; it does not carry the ordinary FSD guarantee
  for every increasing utility.

- plot:

  logical; `TRUE` (default) plots the FSD test.

## Value

With `type = "discrete"`, returns `"X FSD Y"`, `"Y FSD X"`, or
`"NO FSD EXISTS"`. With `type = "continuous"`, returns `"X OMEGA Y"`,
`"Y OMEGA X"`, or `"NO OMEGA DOMINANCE EXISTS"`; the numerical
area-order criterion is unchanged.

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
[doi:10.4236/jmf.2016.61012](https://doi.org/10.4236/jmf.2016.61012) .

Viole, F. (2017) "A Note on Stochastic Dominance."
[doi:10.2139/ssrn.3002675](https://doi.org/10.2139/ssrn.3002675) .

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.FSD(x, y)
} # }
```
