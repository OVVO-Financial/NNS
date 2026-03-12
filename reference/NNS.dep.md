# NNS Dependence

Returns the dependence and nonlinear correlation between two variables
based on higher order partial moment matrices measured by frequency or
area.

## Usage

``` r
NNS.dep(x, y = NULL, asym = FALSE, p.value = FALSE, print.map = FALSE)
```

## Arguments

- x:

  a numeric vector, matrix or data frame.

- y:

  `NULL` (default) or a numeric vector with compatible dimensions to
  `x`.

- asym:

  logical; `FALSE` (default) Allows for asymmetrical dependencies.

- p.value:

  logical; `FALSE` (default) Generates 100 independent random
  permutations to test results against and plots 95 percent confidence
  intervals along with all results.

- print.map:

  logical; `FALSE` (default) Plots quadrant means, or p-value
  replicates.

## Value

Returns the bi-variate `"Correlation"` and `"Dependence"` or correlation
/ dependence matrix for matrix input.

## Note

For asymmetrical `(asym = TRUE)` matrices, directional dependence is
returned as (\[column variable\] —\> \[row variable\]).

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.dep(x, y)

## Correlation / Dependence Matrix
x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
B <- cbind(x, y, z)
NNS.dep(B)
} # }
```
