[![packageversion](https://img.shields.io/badge/NNS%20version-10.5-blue.svg?style=flat-square)](https://github.com/OVVO-Financial/NNS/commits/NNS-Beta-Version)   [![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)


# NNS
Nonlinear nonparametric statistics using partial moments.  Partial moments are the [elements of variance](https://www.linkedin.com/pulse/elements-variance-fred-viole) and [asymptotically approximate the area of f(x)](https://www.ssrn.com/abstract=2186471).  These robust statistics provide the basis for nonlinear analysis while retaining linear equivalences.

NNS offers: 
  - Numerical Integration & Numerical Differentiation
  - Partitional & Hierarchial Clustering
  - Nonlinear Correlation & Dependence
  - Causal Analysis
  - Nonlinear Regression & Classification
  - ANOVA
  - Seasonality & Autoregressive Modeling
  - Normalization 
  - Stochastic Dominance
  - Advanced Monte Carlo Sampling


Companion R-package and datasets to: 
#### Viole, F. and Nawrocki, D. (2013) ["*Nonlinear Nonparametric Statistics: Using Partial Moments*"]( https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp)


#### For a quantitative finance implementation of NNS, see [OVVO Labs](https://www.ovvolabs.com)



## Current Version
[![NNS](https://img.shields.io/badge/NNS%3E%3D-0.5.5-blue.svg)](https://cran.r-project.org/package=NNS) is built on [![data.table](https://img.shields.io/badge/data.table%3E%3D-1.10.4-6666ff.svg)](https://cran.r-project.org/package=data.table) and [![doParallel](https://img.shields.io/badge/doParallel%3E%3D-1.0.14-6666ff.svg)](https://cran.r-project.org/package=doParallel) architecture 
and [![NNS](https://img.shields.io/badge/NNS%3E%3D-0.9.0-blue.svg)](https://cran.r-project.org/package=NNS) is built on [![Rcpp](https://img.shields.io/badge/Rcpp%3E%3D-1.0.8.3-6666ff.svg)](https://cran.r-project.org/package=Rcpp) with notable performance enhancements.

*Current [![NNS](https://img.shields.io/badge/NNS--blue.svg)](https://cran.r-project.org/package=NNS) CRAN version is  [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/NNS)](https://cran.r-project.org/package=NNS)

## Installation
[![NNS](https://img.shields.io/badge/NNS--blue.svg)](https://cran.r-project.org/package=NNS) requires [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/).  See https://cran.r-project.org/ or [![installr](https://img.shields.io/badge/installr-0.18.0-blue.svg)](https://cran.r-project.org/package=installr) for upgrading to latest R release.

```r
library(remotes); remotes::install_github('OVVO-Financial/NNS', ref = "NNS-Beta-Version")
```
or via CRAN
```r
install.packages('NNS')
```

## Examples
Please see https://github.com/OVVO-Financial/NNS/blob/NNS-Beta-Version/examples/index.md for basic partial moments equivalences and hands-on statistics, machine learning and econometrics examples.


## Citation
```
@Manual{,
    title = {NNS: Nonlinear Nonparametric Statistics},
    author = {Fred Viole},
    year = {2016},
    note = {R package version 10.5},
    url = {https://CRAN.R-project.org/package=NNS},
  }
```

## Thank you for your interest in NNS!
![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/NNS)
