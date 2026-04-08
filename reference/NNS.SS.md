# NNS Stochastic Superiority

Computes stochastic superiority between two numeric vectors as the
empirical probability that an observation from `x` exceeds an
observation from `y`, with optional tie adjustment and optional
confidence intervals via maximum entropy bootstrap.

## Usage

``` r
NNS.SS(
  x,
  y,
  confidence.interval = FALSE,
  reps = 999,
  ci = 0.95,
  rho = 1
)
```

## Arguments

- x:

  a numeric vector.

- y:

  a numeric vector.

- confidence.interval:

  logical; `FALSE` (default) returns only the empirical stochastic
  superiority measures. Set to `TRUE` to compute bootstrap confidence
  intervals for `p_star`.

- reps:

  numeric; number of maximum entropy bootstrap replicates used when
  `confidence.interval = TRUE`. Default is `999`.

- ci:

  numeric in \\(0, 1)\\; confidence level used for the bootstrap
  interval when `confidence.interval = TRUE`. Default is `0.95`.

- rho:

  numeric; dependence target passed to
  [`NNS.meboot`](https://OVVO-Financial.github.io/NNS/reference/NNS.meboot.md).
  Default is `1`.

## Value

If `confidence.interval = FALSE`, returns a list containing:

- `p_gt`:

  empirical probability that `x > y`.

- `p_tie`:

  empirical probability that `x = y`.

- `p_star`:

  tie-adjusted stochastic superiority probability.

If `confidence.interval = TRUE`, returns a list containing:

- `p_gt`:

  empirical probability that `x > y`.

- `p_tie`:

  empirical probability that `x = y`.

- `p_star`:

  tie-adjusted stochastic superiority probability.

- `lower`:

  lower confidence bound for `p_star`.

- `upper`:

  upper confidence bound for `p_star`.

- `ci`:

  confidence level used.

- `reps`:

  number of bootstrap replicates used.

- `boot_vals`:

  bootstrap replicate values of `p_star`.

## Details

`NNS.SS` returns: \$\$P(X \> Y),\$\$ the tie probability \$\$P(X =
Y),\$\$ and the tie-adjusted stochastic superiority measure \$\$P^\* =
P(X \> Y) + \frac{1}{2} P(X = Y).\$\$

When `confidence.interval = TRUE`, confidence bounds for `P^*` are
computed from
[`NNS.meboot`](https://OVVO-Financial.github.io/NNS/reference/NNS.meboot.md)
bootstrap replicates using
[`LPM.VaR`](https://OVVO-Financial.github.io/NNS/reference/LPM.VaR.md)
and
[`UPM.VaR`](https://OVVO-Financial.github.io/NNS/reference/UPM.VaR.md)
with `degree = 0`.

Missing values are removed from both `x` and `y` using
[`stats::na.omit`](https://rdrr.io/r/stats/na.fail.html). The empirical
estimates are computed via a fast sorted comparison routine rather than
explicit pairwise expansion of all `x`-`y` combinations.

For continuous data, `p_tie` will typically be zero, so `p_star` and
`p_gt` will be identical up to numerical precision. For discrete data,
`p_star` provides the standard tie-adjusted superiority measure.

When `confidence.interval = TRUE`, the interval is constructed from the
empirical bootstrap distribution of `p_star`, where \\\alpha = 1 - ci\\.
The lower bound is obtained from
[`LPM.VaR`](https://OVVO-Financial.github.io/NNS/reference/LPM.VaR.md)
evaluated at \\\alpha / 2\\, and the upper bound is obtained from
[`UPM.VaR`](https://OVVO-Financial.github.io/NNS/reference/UPM.VaR.md)
evaluated at \\\alpha / 2\\, both with `degree = 0`.

## Note

This function measures stochastic superiority as a pairwise exceedance
probability. This is distinct from first-, second-, or third-degree
stochastic dominance; see
[`NNS.FSD`](https://OVVO-Financial.github.io/NNS/reference/NNS.FSD.md),
[`NNS.SSD`](https://OVVO-Financial.github.io/NNS/reference/NNS.SSD.md),
and
[`NNS.TSD`](https://OVVO-Financial.github.io/NNS/reference/NNS.TSD.md)
for dominance testing.

## References

- Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank
  Correlations in Maximum Entropy Bootstrap and Improved Monte Carlo
  Simulations.
  [doi:10.2139/ssrn.3621614](https://doi.org/10.2139/ssrn.3621614)

- Viole, F. and Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics:
  Using Partial Moments*. ISBN: 1490523995, 2nd edition:
  <https://ovvo-financial.github.io/NNS/book/>.

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(200, mean = 0.4, sd = 1)
y <- rnorm(200, mean = 0.0, sd = 1)

# Empirical stochastic superiority
NNS.SS(x, y)

# With confidence intervals
NNS.SS(x, y, confidence.interval = TRUE, reps = 999, ci = 0.95)

# Discrete example with ties
x <- sample(1:5, 100, replace = TRUE)
y <- sample(1:5, 100, replace = TRUE)
NNS.SS(x, y)
} # }
```
