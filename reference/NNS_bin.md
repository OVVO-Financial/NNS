# Fast binning of numeric vector into equidistant bins

Missing values (NA, Inf, NaN) are added at the end of the vector as the
last bin returned if missinglast is set to TRUE

## Usage

``` r
NNS_bin(x, width, origin = 0, missinglast = FALSE)
```

## Arguments

- x:

  A matrix of regressor variables. Must have the same number of rows as
  the length of y.

- width:

  The width of the bins

- origin:

  The starting point for the bins. Any number smaller than origin will
  be disregarded

- missinglast:

  Boolean. Should the missing observations be added as a separate
  element at the end of the returned count vector.

## Value

An list with elements counts (the frequencies), origin (the origin),
width (the width), missing (the number of missings), and
last_bin_is_missing (boolean) telling whether the missinglast is true or
not.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
x <- sample(10, 20, replace = TRUE)
NNS_bin(x, 15)
} # }
```
