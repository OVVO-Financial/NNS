# NNS Seasonality Test

Seasonality test based on the coefficient of variation for the variable
and lagged component series. A result of 1 signifies no seasonality
present.

## Usage

``` r
NNS.seas(variable, modulo = NULL, mod.only = TRUE, plot = TRUE)
```

## Arguments

- variable:

  a numeric vector.

- modulo:

  integer(s); NULL (default) Used to find the nearest multiple(s) in the
  reported seasonal period.

- mod.only:

  logical; `TRUE` (default) Limits the number of seasonal periods
  returned to the specified `modulo`.

- plot:

  logical; `TRUE` (default) Returns the plot of all periods exhibiting
  seasonality and the variable level reference.

## Value

Returns a matrix of all periods exhibiting less coefficient of variation
than the variable with `"all.periods"`; and the single period exhibiting
the least coefficient of variation versus the variable with
`"best.period"`; as well as a vector of `"periods"` for easy call into
[NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md).
If no seasonality is detected, `NNS.seas` will return ("No Seasonality
Detected").

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)

## To call strongest period based on coefficient of variation:
NNS.seas(x, plot = FALSE)$best.period

## Using modulos for logical seasonal inference:
NNS.seas(x, modulo = c(2,3,5,7), plot = FALSE)
} # }
```
