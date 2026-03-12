# Upper Partial Moment Ratio

This function generates a standardized univariate upper partial moment
of any non‑negative degree for a given target.

## Usage

``` r
UPM.ratio(degree, target, variable)
```

## Arguments

- degree:

  numeric; degree = 0 gives frequency, degree = 1 gives area.

- target:

  numeric vector; threshold(s). Defaults to mean(variable).

- variable:

  numeric vector or data‑frame column to evaluate.

## Value

Numeric vector of standardized upper partial moments.

## References

Viole, F. & Nawrocki, D. (2013) \*Nonlinear Nonparametric Statistics:
Using Partial Moments\* (ISBN:1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples
