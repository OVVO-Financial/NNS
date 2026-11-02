# NNS ANOVA: Nonparametric Analysis of Variance

Performs a distribution-free ANOVA using partial-moment statistics to
assess differences between control and treatment groups. Depending on
the setting of `means.only`, the procedure tests either differences in
central tendency (means or medians) or differences across the full
empirical distributions.

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
  for mean-only comparison

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

- `Control_Statistic`: Mean/median of control group

- `Treatment_Statistic`: Mean/median of treatment group

- `Grand_Statistic`: Mean/median of all observations pooled

- `Control_CDF`: CDF value at grand statistic (control)

- `Treatment_CDF`: CDF value at grand statistic (treatment)

- `Certainty`: Similarity index in \[0, 1\] (1 = the groups are the
  *same*; means-only or full distribution depending on `means.only`);
  not a calibrated probability.

- `Effect_Size_LB`: Lower bound of treatment effect (if
  confidence.interval requested)

- `Effect_Size_UB`: Upper bound of treatment effect (if
  confidence.interval requested)

- `Confidence_Level`: Confidence level used (if confidence.interval
  requested)

## Details

The key output is the `Certainty` metric, a descriptive similarity index
in \\\[0, 1\]\\ (1 = the groups are the same) with respect to the chosen
comparison mode:

- If `means.only = TRUE`: similarity of the group *means* (or medians,
  if `medians = TRUE`), measured by how far each group's partial-moment
  CDF at the grand statistic is from 0.5.

- If `means.only = FALSE`: similarity of the two *entire distributions*,
  adding the same comparison at the 12.5% and 25% tail quantiles.

`Certainty` is *not* a calibrated probability and not the inverse of a
p-value: for two samples from the same distribution it is not uniform on
\\\[0, 1\]\\ and rises toward 1 with sample size, it is shrunk by the
population adjustment \\((n_1 + n_2 - 2) / (n_1 + n_2))^2\\ for small
samples, and it should not be compared with p-value thresholds such as
0.05. `robust = TRUE` bounds the sampling variability of `Certainty`
itself; it is not a reference distribution under no difference. For
that, recompute `Certainty` on random relabelings of the pooled
observations.

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
