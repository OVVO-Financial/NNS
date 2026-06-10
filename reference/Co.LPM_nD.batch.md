# Batched Co-Lower Partial Moment nD

Internal batched backend for evaluating \codeCo.LPM_nD over many
targets.

## Usage

``` r
Co.LPM_nD.batch(data, targets, degree = 0, norm = TRUE)
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

## Value

Numeric vector, one value per row of targets.
