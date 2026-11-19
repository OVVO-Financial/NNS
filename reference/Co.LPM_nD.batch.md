# Batched Co-Lower Partial Moment nD

Internal batched backend for evaluating `Co.LPM_nD` over many targets.

## Usage

``` r
Co.LPM_nD.batch(data, targets, degree = 0, norm = TRUE, strict = FALSE)
```

## Arguments

- data:

  A numeric matrix with observations in rows and variables in columns.

- targets:

  A numeric matrix with target rows and the same number of columns as
  data.

- degree:

  numeric; degree for lower deviations.

- norm:

  logical; normalize result.

- strict:

  logical; count only rows strictly below the target in the degree-0
  share.

## Value

Numeric vector, one value per row of targets.
