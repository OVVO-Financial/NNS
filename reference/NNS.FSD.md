# NNS FSD Test

Bi-directional test of first degree stochastic dominance using lower
partial moments.

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

  options: ("discrete", "continuous"); `"discrete"` (default) selects
  the type of CDF.

- plot:

  logical; `TRUE` (default) plots the FSD test.

## Value

Returns one of the following FSD results: `"X FSD Y"`, `"Y FSD X"`, or
`"NO FSD EXISTS"`.

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
