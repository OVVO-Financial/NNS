# NNS TSD Test uni-directional

Uni-directional test of third degree stochastic dominance using lower
partial moments used in SD Efficient Set routine.

## Usage

``` r
NNS.TSD.uni(x, y)
```

## Arguments

- x:

  a numeric vector.

- y:

  a numeric vector.

## Value

Returns (1) if `"X TSD Y"`, else (0).

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
NNS.TSD.uni(x, y)
} # }
```
