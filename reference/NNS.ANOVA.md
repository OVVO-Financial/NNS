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

  logical; `FALSE` (default) Generates 100 independent random
  permutations to test results, and returns / plots 95 percent
  confidence intervals along with robust central tendency of all results
  for pairwise analysis only.

## Value

Returns a list containing:

- `Control_Statistic`: Mean/median of control group

- `Treatment_Statistic`: Mean/median of treatment group

- `Grand_Statistic`: Grand mean/median

- `Control_CDF`: CDF value at grand statistic (control)

- `Treatment_CDF`: CDF value at grand statistic (treatment)

- `Certainty`: Probability that the groups are the *same* (means-only or
  full distribution depending on `means.only`).

- `Effect_Size_LB`: Lower bound of treatment effect (if
  confidence.interval requested)

- `Effect_Size_UB`: Upper bound of treatment effect (if
  confidence.interval requested)

- `Confidence_Level`: Confidence level used (if confidence.interval
  requested)

## Details

The key output is the `Certainty` metric, a calibrated probability in
\\\[0, 1\]\\ representing the likelihood that the groups being compared
are the \*same\* with respect to the chosen comparison mode:

- If `means.only = TRUE`: `Certainty` is the probability that the group
  *means* (or medians, if `medians = TRUE`) are the same.

- If `means.only = FALSE`: `Certainty` is the probability that the two
  *entire distributions* are the same.

This makes `Certainty` the conceptual inverse of a classical p-value. A
\*low\* Certainty (e.g., \< 0.10) indicates strong evidence of
difference, while a \*high\* Certainty (e.g., \> 0.90) indicates strong
evidence of similarity.

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
