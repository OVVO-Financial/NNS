# NNS mode

Mode of a distribution, either continuous or discrete.

## Usage

``` r
NNS.mode(x, discrete = FALSE, multi = TRUE)
```

## Arguments

- x:

  vector of data.

- discrete:

  logical; `FALSE` (default) for discrete distributions.

- multi:

  logical; `TRUE` (default) returns multiple mode values.

## Value

Returns a numeric value representing the mode of the distribution.

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)
NNS.mode(x)
} # }
```
