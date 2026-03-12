# NNS rescale

Rescale a vector using either min-max scaling or risk-neutral
adjustment.

## Usage

``` r
NNS.rescale(x, a, b, method = "minmax", T = NULL, type = "Terminal")
```

## Arguments

- x:

  numeric vector; data to rescale (e.g., terminal prices for
  risk-neutral method).

- a:

  numeric; defines the scaling target: - For `method = "minmax"`: the
  lower limit of the output range (e.g., 5 to scale to \[5, b\]). - For
  `method = "riskneutral"`: the initial price \\ S_0 \\ (must be
  positive, e.g., 100), used to set the target mean.

- b:

  numeric; defines the scaling range or rate: - For `method = "minmax"`:
  the upper limit of the output range (e.g., 10 to scale to \[a,
  10\]). - For `method = "riskneutral"`: the risk-free rate \\ r \\
  (e.g., 0.05), used with \\ T \\ to adjust the mean.

- method:

  character; scaling method: `"minmax"` (default) for min-max scaling,
  or `"riskneutral"` for risk-neutral adjustment.

- T:

  numeric; time to maturity in years (required for
  `method = "riskneutral"`, ignored otherwise; e.g., 1). Default is
  NULL.

- type:

  character; for `method = "riskneutral"`: `"Terminal"` (default) or
  `"Discounted"` (mean = \\ S_0 \\).

## Value

Returns a rescaled distribution: - For `"minmax"`: values scaled
linearly to the range `[a, b]`. - For `"riskneutral"`: values scaled
multiplicatively to a risk-neutral mean (\\ S_0 e^(rT) \\ if
`type = "Terminal"`, or \\ S_0 \\ if `type = "Discounted"`).

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
# Min-max scaling: a = lower limit, b = upper limit
x <- rnorm(100)
NNS.rescale(x, a = 5, b = 10, method = "minmax")  # Scales to [5, 10]

# Risk-neutral scaling (Terminal): a = S_0, b = r  # Mean approx 105.13
prices <- 100 * exp(cumsum(rnorm(100, 0.001, 0.02)))
NNS.rescale(prices, a = 100, b = 0.05, method = "riskneutral", T = 1, type = "Terminal")

# Risk-neutral scaling (Discounted): a = S_0, b = r  # Mean approx 100
NNS.rescale(prices, a = 100, b = 0.05, method = "riskneutral", T = 1, type = "Discounted")
} # }
```
