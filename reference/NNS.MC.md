# NNS Monte Carlo Sampling

Monte Carlo sampling from the maximum entropy bootstrap routine
[NNS.meboot](https://OVVO-Financial.github.io/NNS/reference/NNS.meboot.md),
ensuring the replicates are sampled from the full \[-1,1\] correlation
space.

## Usage

``` r
NNS.MC(
  x,
  reps = 30,
  lower_rho = -1,
  upper_rho = 1,
  by = 0.01,
  exp = 1,
  type = "spearman",
  drift = TRUE,
  target_drift = NULL,
  target_drift_scale = NULL,
  xmin = NULL,
  xmax = NULL,
  ...
)
```

## Arguments

- x:

  vector of data.

- reps:

  numeric; number of replicates to generate, `30` default.

- lower_rho:

  numeric `[-1,1]`; `.01` default will set the `from` argument in
  `seq(from, to, by)`.

- upper_rho:

  numeric `[-1,1]`; `.01` default will set the `to` argument in
  `seq(from, to, by)`.

- by:

  numeric; `.01` default will set the `by` argument in
  `seq(-1, 1, step)`.

- exp:

  numeric; `1` default will exponentially weight maximum rho value if
  `exp > 1`. Shrinks values towards `upper_rho`.

- type:

  options("spearman", "pearson", "NNScor", "NNSdep");
  `type = "spearman"`(default) dependence metric desired.

- drift:

  logical; `drift = TRUE` (default) preserves the drift of the original
  series.

- target_drift:

  numerical; `target_drift = NULL` (default) Specifies the desired drift
  when `drift = TRUE`, i.e. a risk-free rate of return.

- target_drift_scale:

  numerical; instead of calculating a `target_drift`, provide a scalar
  to the existing drift when `drift = TRUE`.

- xmin:

  numeric; the lower limit for the left tail.

- xmax:

  numeric; the upper limit for the right tail.

- ...:

  possible additional arguments to be passed to
  [NNS.meboot](https://OVVO-Financial.github.io/NNS/reference/NNS.meboot.md).

## Value

- ensemble average observation over all replicates as a vector.

- replicates maximum entropy bootstrap replicates as a list for each
  `rho`.

## References

Vinod, H.D. and Viole, F. (2020) Arbitrary Spearman's Rank Correlations
in Maximum Entropy Bootstrap and Improved Monte Carlo Simulations.
[doi:10.2139/ssrn.3621614](https://doi.org/10.2139/ssrn.3621614)

## Examples

``` r
if (FALSE) { # \dontrun{
# To generate a set of MC sampled time-series to AirPassengers
MC_samples <- NNS.MC(AirPassengers, reps = 10, lower_rho = -1, upper_rho = 1, by = .5, xmin = 0)
} # }
```
