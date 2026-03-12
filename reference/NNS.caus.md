# NNS Causation

Returns the causality from observational data between two variables.

## Usage

``` r
NNS.caus(
  x,
  y = NULL,
  factor.2.dummy = FALSE,
  tau = 0,
  plot = FALSE,
  p.value = FALSE,
  nperm = 100L,
  permute = c("y", "x", "both"),
  seed = NULL,
  conf.int = 0.95
)
```

## Arguments

- x:

  a numeric vector, matrix or data frame.

- y:

  `NULL` (default) or a numeric vector with compatible dimensions to
  `x`.

- factor.2.dummy:

  logical; `FALSE` (default) Automatically augments variable matrix with
  numerical dummy variables based on the levels of factors. Includes
  dependent variable `y`.

- tau:

  options: ("cs", "ts", integer); 0 (default) Number of lagged
  observations to consider (for time series data). Otherwise, set
  `(tau = "cs")` for cross-sectional data. `(tau = "ts")` automatically
  selects the lag of the time series data, while `(tau = [integer])`
  specifies a time series lag.

- plot:

  logical; `FALSE` (default) Plots the raw variables, tau normalized,
  and cross-normalized variables.

- p.value:

  logical; `FALSE` (default) If `TRUE`, runs a permutation test to
  compute empirical p-values for the signed causation from x -\> y.

- nperm:

  integer; number of permutations to use when `p.value = TRUE`. Default
  100.

- permute:

  one of "both", "y", or "x"; which variable(s) to shuffle when
  constructing the null distribution.

- seed:

  optional integer seed for reproducibility of the permutation test.

- conf.int:

  numeric; 0.95 (default) confidence level for the partial-moment based
  interval computed on the permutation null distribution.

## Value

If `p.value=FALSE` returns the original causation vector of length 3
(directional given/received and net), named either "C(x—\>y)" or
"C(y—\>x)" in the third slot. If `p.value=TRUE` returns a list with
components: \* `causation`: the original causation vector as above. \*
`p.value`: a list with empirical two-sided and one-sided p-values
(x_causes_y, y_causes_x), the null distribution, the observed signed
statistic, and metadata (permute, nperm). If `p.value=TRUE` for a
matrix, the function returns a list with components: \* `causality`: the
causality matrix. \* `lower_CI`: matrix of lower confidence bounds
(partial-moment based). \* `upper_CI`: matrix of upper confidence bounds
(partial-moment based). \* `p.value`: matrix of empirical two-sided
p-values.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
## x causes y...
set.seed(123)
x <- rnorm(1000) ; y <- x ^ 2
NNS.caus(x, y, tau = "cs")

## Causal matrix without per factor causation
NNS.caus(iris, tau = 0)

## Causal matrix with per factor causation
NNS.caus(iris, factor.2.dummy = TRUE, tau = 0)
} # }
```
