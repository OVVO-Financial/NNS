# NNS

![](articles/images/NNS_hex_sticker.png)

NNS (Nonlinear Nonparametric Statistics) leverages partial moments – the
fundamental [elements of
variance](https://github.com/OVVO-Financial/NNS/blob/NNS-Beta-Version/examples/Partial%20Moments%20Equivalences.md)
that [asymptotically approximate the area of
f(x)](https://ovvo-financial.github.io/NNS/book/numerical-integration-via-partial-moments.html)
– to provide a robust foundation for nonlinear analysis while
maintaining linear equivalences. Designed for real-world data that
violates symmetry, linearity, or distributional assumptions.

NNS delivers a comprehensive suite of advanced statistical techniques,
including: - Numerical Integration & Numerical Differentiation -
Partitional & Hierarchical Clustering - Nonlinear Correlation &
Dependence - Causal Analysis - Nonlinear Regression & Classification -
ANOVA - Seasonality & Autoregressive Modeling - Normalization -
Stochastic Superiority / Dominance - Advanced Monte Carlo Sampling

Companion R-package and datasets to: \#### Viole, F. and Nawrocki, D.
(2013) “*Nonlinear Nonparametric Statistics: Using Partial Moments*”
(ISBN: 1490523995)

2nd edition available here: <https://ovvo-financial.github.io/NNS/book/>

#### For a direct quantitative finance implementation of NNS, see [OVVO Labs](https://www.ovvolabs.com)

## Current Version

Current
[![NNS](https://img.shields.io/badge/NNS--blue.svg)](https://cran.r-project.org/package=NNS)
CRAN version is
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/NNS)](https://www.r-pkg.org/badges/version/NNS)

## Installation

[![NNS](https://img.shields.io/badge/NNS--blue.svg)](https://cran.r-project.org/package=NNS)
requires [![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/).
See <https://cran.r-project.org/> or
[![installr](https://img.shields.io/badge/installr-0.18.0-blue.svg)](https://cran.r-project.org/package=installr)
for upgrading to the latest R release.

``` r

library(remotes); remotes::install_github('OVVO-Financial/NNS', ref = "NNS-Beta-Version")
```

or via CRAN

``` r

install.packages('NNS')
```

## Examples

The nine numbered files under
[`vignettes/`](https://OVVO-Financial.github.io/NNS/vignettes/) are the
canonical NNS example curriculum and the source of truth for companion
language ports. Python follows the same 01–09 sequence while using
Python-native syntax and containers.

See
[`examples/index.md`](https://OVVO-Financial.github.io/NNS/examples/index.md)
for the canonical cross-language mapping plus applied studies in
statistics, regression, machine learning, forecasting, and econometrics.

## Citation

    @Manual{,
        title = {NNS: Nonlinear Nonparametric Statistics},
        author = {Fred Viole},
        year = {2016},
        note = {R package version 13.1},
        url = {https://CRAN.R-project.org/package=NNS},
      }

## Thank you for your interest in NNS!

![](https://cranlogs.r-pkg.org/badges/NNS)![](https://cranlogs.r-pkg.org/badges/grand-total/NNS)
