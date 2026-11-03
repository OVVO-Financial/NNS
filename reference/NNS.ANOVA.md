# NNS ANOVA: Nonparametric Analysis of Variance

Compares control and treatment groups through benchmark-relative partial
moments. The chosen mode compares means, medians, or a distribution
profile with central and tail benchmarks, without selecting a parametric
family.

## Usage

``` r
NNS.ANOVA(
  control,
  treatment,
  means.only = FALSE,
  medians = FALSE,
  confidence.interval = 0.95,
  tails = "Both",
  pairwise = FALSE,
  plot = TRUE,
  robust = FALSE
)
```

## Arguments

- control:

  Numeric vector of control group observations

- treatment:

  Numeric vector of treatment group observations

- means.only:

  Logical; `FALSE` (default) uses full distribution analysis. Set `TRUE`
  for the central comparison only (means, or medians when
  `medians = TRUE`)

- medians:

  Logical; `FALSE` (default) uses means. Set `TRUE` for median-based
  analysis

- confidence.interval:

  Numeric \[0,1\]; confidence level for effect size bounds (e.g., 0.95)

- tails:

  Character; specifies CI tail(s): "both", "left", or "right"

- pairwise:

  logical; `FALSE` (default) Returns pairwise certainty tests when set
  to `pairwise = TRUE`.

- plot:

  Logical; `TRUE` (default) generates distribution plot

- robust:

  logical; `FALSE` (default) Two-sample (`control`, `treatment`) path
  only: draws 100 independent bootstrap resamples of each group (with
  replacement, from its own full sample), and returns / plots the robust
  central tendency
  ([NNS.gravity](https://OVVO-Financial.github.io/NNS/reference/NNS.gravity.md))
  and the `confidence.interval` bounds of the resampled `Certainty`
  values. It is ignored when `control` is a matrix or list (multi-group,
  including `pairwise = TRUE`).

## Value

Returns a list containing:

- `Control`: Mean/median of control group

- `Treatment`: Mean/median of treatment group

- `Grand_Statistic`: Mean/median of all observations pooled

- `Control_CDF`: CDF value at grand statistic (control)

- `Treatment_CDF`: CDF value at grand statistic (treatment)

- `Certainty`: Agreement index for the chosen central or
  distribution-profile comparison, with the finite-sample adjustment;
  not a calibrated probability.

- `Effect_Size_LB`: Lower bound of treatment effect (if
  confidence.interval requested)

- `Effect_Size_UB`: Upper bound of treatment effect (if
  confidence.interval requested)

- `Confidence_Level`: Confidence level used (if confidence.interval
  requested)

## Details

The key output, `Certainty`, is a descriptive agreement index in \\\[0,
1\]\\ for the chosen target:

- `means.only = TRUE, medians = FALSE`: compares group means using
  degree-one deviation-area ratios at the pooled mean.

- `means.only = TRUE, medians = TRUE`: compares medians using the
  distance of 0.5 from each interval \\\[F(t-), F(t)\]\\ at the pooled
  median \\t\\. This interval contains 0.5 exactly when \\t\\ is a
  median, including tied and discrete data.

- `means.only = FALSE`: adds the 12.5% and 25% lower and upper
  deviation-area benchmarks to the chosen central comparison. This
  examines distributional shape through a finite profile; agreement at
  these benchmarks does not prove equality of laws.

At a constant group's only value both deviation areas vanish. For
agreement comparisons this satisfies the expectile balance at every
level, so its area discrepancy is zero there. Returned CDF values retain
the public partial-moment step convention; other benchmarks and
nonconstant samples use the ordinary ratio discrepancy.

The raw score is 1 at maximal benchmark agreement. The returned score is
multiplied by \\((n_1 + n_2 - 2)/(n_1 + n_2))^2\\, so the finite-sample
maximum is below 1. `Certainty` is not a calibrated probability, a
p-value, or statistical power, and has no universal significance cutoff.
Equal population means need not produce the same score distribution when
sample sizes or spreads differ. `robust = TRUE` describes sampling
variability of the score; it does not generate a null reference
distribution.

Random relabelings of the pooled observations provide a reference under
identical distributions only when observations are exchangeable. They
are not a general unequal-variance equal-means test. Any inferential
calibration must match the chosen target and the sampling design;
dependent observations require a suitable resampling design.

Effect-size bounds retain the conservative difference of the groups'
separate bootstrap endpoints. They are not quantiles of the bootstrap
difference and need not have exactly the requested nominal coverage. The
independent-group bootstrap assumes independent observations within and
between groups.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995, 2nd edition:
<https://ovvo-financial.github.io/NNS/book/>)

Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"
[doi:10.2139/ssrn.3007373](https://doi.org/10.2139/ssrn.3007373)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
 if (FALSE) { # \dontrun{
### Binary analysis and effect size
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
NNS.ANOVA(control = x, treatment = y)

### Two variable analysis with no control variable
A <- cbind(x, y)
NNS.ANOVA(A)

### Medians test
NNS.ANOVA(A, means.only = TRUE, medians = TRUE)

### Multiple variable analysis with no control variable
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
A <- cbind(x, y, z)
NNS.ANOVA(A)

### Different length vectors used in a list
x <- rnorm(30) ; y <- rnorm(40) ; z <- rnorm(50)
A <- list(x, y, z)
NNS.ANOVA(A)
} # }
```
