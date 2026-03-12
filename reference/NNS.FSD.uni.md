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

  options: ("discrete", "continuous"); `"discrete"` (default) selects
  the type of CDF.

## Value

Returns (1) if `"X FSD Y"`, else (0).

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
