# Lower Partial Moment Ratio

This function generates a standardized univariate lower partial moment
of any non‑negative degree for a given target.

## Usage

``` r
LPM.ratio(degree, target, variable)
```

## Arguments

- degree:

  numeric; degree = 0 gives frequency (CDF), degree = 1 gives area.

- target:

  numeric vector; threshold(s). Defaults to mean(variable).

- variable:

  numeric vector or data‑frame column to evaluate.

## Value

Numeric vector of standardized lower partial moments.

## References

Viole, F. & Nawrocki, D. (2013) \*Nonlinear Nonparametric Statistics:
Using Partial Moments\* (ISBN:1490523995)

Viole, F. (2017) Continuous CDFs and ANOVA with NNS.
[doi:10.2139/ssrn.3007373](https://doi.org/10.2139/ssrn.3007373)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
  set.seed(123)
  x <- rnorm(100)
  LPM.ratio(0, mean(x), x)
#> [1] 0.51
if (FALSE) { # \dontrun{
  plot(sort(x), LPM.ratio(0, sort(x), x))
  plot(sort(x), LPM.ratio(1, sort(x), x))
} # }
```
