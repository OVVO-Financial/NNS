# NNS TSD Test

Bi-directional test of third degree stochastic dominance using lower
partial moments.

## Usage

``` r
NNS.TSD(x, y, plot = TRUE)
```

## Arguments

- x:

  a numeric vector.

- y:

  a numeric vector.

- plot:

  logical; `TRUE` (default) plots the TSD test.

## Value

Returns one of the following TSD results: `"X TSD Y"`, `"Y TSD X"`, or
`"NO TSD EXISTS"`.

## References

Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the
Computation of the SD Efficient Set." Journal of Mathematical Finance,
6, 105-126.
[doi:10.4236/jmf.2016.61012](https://doi.org/10.4236/jmf.2016.61012) .

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.TSD(x, y)
} # }
```
