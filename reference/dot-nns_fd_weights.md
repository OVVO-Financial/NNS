# Finite-difference weights for a local polynomial stencil

Construct polynomial finite-difference weights on distinct local nodes.
Invalid or degenerate stencils return NA coefficients.

## Usage

``` r
.nns_fd_weights(nodes, at, degree = 1L)
```

## Arguments

- nodes:

  Numeric vector of distinct stencil coordinates.

- at:

  Numeric evaluation coordinate.

- degree:

  Integer derivative order; defaults to 1.

## Value

Numeric vector of finite-difference weights, or NA coefficients when the
stencil is invalid.
