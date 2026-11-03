# NNS Causation

Returns the directional causation between two variables from
observational data, as an index in \[0, 1\] for each direction (0.5 = no
net direction).

## Usage

``` r
NNS.caus(
  x,
  y = NULL,
  factor.2.dummy = FALSE,
  tau = 0,
  plot = FALSE,
  p.value = FALSE,
  nperm = 100L,
  permute = c("auto", "y", "x", "both", "shift"),
  seed = NULL,
  conf.int = 0.95
)
```

## Arguments

- x:

  a numeric vector, matrix or data frame.

- y:

  `NULL` (default) or a numeric vector with compatible dimensions to
  `x`.

- factor.2.dummy:

  logical; `FALSE` (default) Automatically augments variable matrix with
  numerical dummy variables based on the levels of factors. Includes
  dependent variable `y`.

- tau:

  options: ("cs", "ts", integer); 0 (default) Lag between cause and
  effect. `(tau = [integer])` compares the lagged pairs
  `(x[t - tau], y[t])` for x —\> y and `(y[t - tau], x[t])` for y —\> x.
  `(tau = 0)` or `(tau = "cs")` compares the contemporaneous pairs
  (cross-sectional data). `(tau = "ts")` uses one lag for both
  directions, the mean of the first
  [NNS.seas](https://OVVO-Financial.github.io/NNS/reference/NNS.seas.md)
  period at most `sqrt(length(x))` of each series.

- plot:

  logical; `FALSE` (default) Plots both lagged pairs on the copula scale
  with their indices.

- p.value:

  logical; `FALSE` (default) If `TRUE`, runs a permutation test of the
  signed statistic `C(x--->y) - C(y--->x)`.

- nperm:

  integer; number of permutations to use when `p.value = TRUE`. Default
  100.

- permute:

  one of "auto" (default), "y", "x", "both" or "shift"; how the null
  distribution is built. `"shift"` circularly shifts `y` by random
  distinct offsets at least `max(tau + 1, n / 10)` steps from alignment,
  so each series keeps its own dynamics and only their alignment is
  broken: a nonparametric counterpart to conditioning on the effect's
  own past. When at most `nperm` offsets are admissible, each is used
  once and the reported `nperm` is their number. `"y"`, `"x"` and
  `"both"` shuffle the values, an independence null for cross-sections.
  `"auto"` uses `"shift"` when the lag is positive (including
  `tau = "ts"`) and `"y"` otherwise. The lag of `tau = "ts"` is resolved
  once on the observed series and kept for every permutation. For a
  matrix `x` each choice applies to every pair, with column `i` as `x`
  and column `j` as `y`; `"auto"` shuffles every column for
  cross-sections.

- seed:

  optional integer seed for reproducibility of the permutation test.

- conf.int:

  numeric; 0.95 (default) confidence level for the partial-moment based
  interval computed on the permutation null distribution.

## Value

If `p.value = FALSE` returns a named vector
`c("C(x--->y)", "C(y--->x)")` of directional causation indices in \[0,
1\] that sum to 1: above 0.5 is evidence that the first direction leads,
0.5 is no net direction. If `p.value = TRUE` returns a list with
components: \* `causation`: the vector above. \* `p.value`: a list with
empirical two-sided and one-sided p-values (x_causes_y, y_causes_x) of
the signed statistic `C(x--->y) - C(y--->x)`, its null distribution and
observed value, the partial-moment interval of the null, and metadata
(permute, nperm). For a matrix `x`, entry `[i, j]` is
`C(column i ---> column j)`, with a 0 diagonal. If `p.value = TRUE` for
a matrix, the function returns a list with components: \* `causality`:
the causation matrix. \* `lower_CI`, `upper_CI`: matrices with the
partial-moment interval of each index under the permutation null. \*
`p.value_matrix`: matrix of empirical two-sided p-values of
`C(i ---> j) - C(j ---> i)`.

## Details

Following book chapter 14, each direction is evaluated on explicitly
lagged pairs mapped to the copula scale with per-column ranks (the
degree-0 lower partial moment), so the result is invariant to shifts,
rescaling and any increasing transform of either variable (reversing a
variable's order changes the copula, and so the index). For x —\> y the
raw score is

\$\$\tilde{C}\_{x \to y} = \frac{1}{2}\left\[P\_{x \to y} +
\max(\rho\_{x \to y} - \rho\_{y \to x}, 0)\right\],\$\$

where \\P\_{x \to y}\\ is the share of the deviation of `y[t]` resolved
by conditioning on the position of `x[t - tau]`, the deviation being the
degree-2 lower plus upper partial moment about the mean with the
\\m/(m-1)\\ factor, which is exactly unbiased under permutation for any
distribution, and \\\rho\_{x \to y}\\ is the debiased directional
dependence of `y[t]` on `x[t - tau]` from
[NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)'s
cells (the lagged co-partial moments of the pair), averaged over the
levels of the partition on `x[t - tau]`.
[NNS.dep](https://OVVO-Financial.github.io/NNS/reference/NNS.dep.md)
reports the strongest level; the average has no selection step, so
independent series are centred at 0.5 whatever their margins. The index
is

\$\$C\_{x \to y} = \frac{\tilde{C}\_{x \to y} + \kappa}{\tilde{C}\_{x
\to y} + \tilde{C}\_{y \to x} + 2\kappa}, \quad C\_{y \to x} = 1 - C\_{x
\to y},\$\$

with a neutral pseudo-count \\\kappa = 1/\sqrt{n - \tau}\\, the order of
the raw scores' sampling variation, so that weak or absent directional
evidence stays near 0.5. An invertible relationship (e.g.
`y = 2 * x + 1`) has no cross-sectional direction and returns 0.5; the
lag is what identifies direction in that case. The index measures
directional lagged association: it does not condition on other
variables, so a common driver that reaches `y` later than `x` makes x
—\> y look causal.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
## x causes y (y is a function of x, x is not a function of y)
set.seed(123)
x <- rnorm(1000) ; y <- x ^ 2
NNS.caus(x, y, tau = "cs")

## Lagged causation: y follows x by one period
x <- rnorm(500) ; y <- c(0, 0.8 * x[-500]) + 0.3 * rnorm(500)
NNS.caus(x, y, tau = 1)

## Causal matrix: entry [i, j] is C(column i ---> column j)
NNS.caus(iris, tau = 0)

## Causal matrix with per factor causation
NNS.caus(iris, factor.2.dummy = TRUE, tau = 0)
} # }
```
