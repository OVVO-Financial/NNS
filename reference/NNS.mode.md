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

  logical; `FALSE` (default) for continuous distributions, `TRUE` for
  discrete distributions.

- multi:

  logical; `TRUE` (default) returns multiple mode values. `FALSE`
  returns the highest peak; among equally high peaks, the one nearest
  [NNS.gravity](https://OVVO-Financial.github.io/NNS/reference/NNS.gravity.md)
  of the data (then the lower value), so the result is always one of the
  modes.

## Value

Returns the mode of the distribution: a numeric vector of all modes when
`multi = TRUE` (default), otherwise a single numeric value.

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
