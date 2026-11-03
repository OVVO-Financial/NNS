# NNS gravity

Central tendency averaging the median, the mean, the midpoint of the
quartiles and a histogram mode; less sensitive to outliers than the
mean.

## Usage

``` r
NNS.gravity(x, discrete = FALSE)
```

## Arguments

- x:

  vector of data.

- discrete:

  logical; `FALSE` (default) for continuous distributions, `TRUE` for
  discrete distributions.

## Value

Returns a numeric value representing the central tendency of the
distribution.

## Details

\$\$(Q_2 + m + \bar{x} + (Q_1 + Q_3)/2) / 4,\$\$ with the quartiles
\\Q_1, Q_2, Q_3\\ of `quantile(x, type = 7)` and, as the mode \\m\\, the
mean of the observations in the tallest bin of a histogram with bin
width IQR / sqrt(n) and its two neighbouring bins (all observations when
several bins are equally tall). Up to three values give the median. The
mean term remains unbounded under one-point contamination: a single 1e6
among 99 zeros gives 2500, against 10000 for the mean. The other
components can also change, so the one-quarter factor is an example, not
a universal sensitivity bound. Its asymptotic replacement breakdown
point is zero; it is not a high-breakdown estimator.

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100)
NNS.gravity(x)
} # }
```
