# Divergent Partial Moment nD

This function generates the aggregate n‑dimensional divergent partial
moment (n \>= 2) for any degree or target.

## Usage

``` r
DPM_nD(data, target, degree = 0, norm = TRUE)
```

## Arguments

- data:

  A numeric matrix with observations in rows and variables in columns.

- target:

  A numeric vector, length equal to ncol(data).

- degree:

  numeric; degree for upper deviations (0 = frequency, 1 = area).

- norm:

  logical; if `TRUE` (default) normalize to the maximum observed value
  (→ \[0,1\]), otherwise return the raw moment.

## Value

Numeric; the n-dimensional divergent partial moment.

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- matrix(rnorm(200), ncol = 4)
DPM_nD(mat, rep(0, ncol(mat)), degree = 1, norm = FALSE)
} # }
```
