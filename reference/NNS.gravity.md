# NNS gravity

Alternative central tendency measure more robust to outliers.

## Usage

``` r
NNS.gravity(x, discrete = FALSE)
```

## Arguments

- x:

  vector of data.

- discrete:

  logical; `FALSE` (default) for discrete distributions.

## Value

Returns a numeric value representing the central tendency of the
distribution.

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)
NNS.gravity(x)
} # }
```
