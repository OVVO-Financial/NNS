# NNS SD Efficient Set

Determines the set of stochastic dominant variables for various degrees.

## Usage

``` r
NNS.SD.efficient.set(x, degree, type = "discrete", status = TRUE)
```

## Arguments

- x:

  a numeric matrix or data frame.

- degree:

  numeric options: (1, 2, 3); Degree of stochastic dominance test from
  (1, 2 or 3).

- type:

  options: ("discrete", "continuous"); `"discrete"` (default) selects
  the type of CDF.

- status:

  logical; `TRUE` (default) Prints status update message in console.

## Value

Returns set of stochastic dominant variable names.

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
x <- rnorm(100) ; y<-rnorm(100) ; z<-rnorm(100)
A <- cbind(x, y, z)
NNS.SD.efficient.set(A, 1)
} # }
```
