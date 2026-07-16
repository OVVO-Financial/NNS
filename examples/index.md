<img src="https://github.com/OVVO-Financial/NNS/raw/NNS-Beta-Version/vignettes/images/NNS_hex_sticker.png" width="150" style="border: none; outline: none; margin: 0; padding: 0; display: block;"/>

# NNS examples

The numbered package vignettes are the **canonical NNS example curriculum**.
They define the statistical narrative, section order, datasets, and intended
interpretation for companion implementations. The Python package follows this
sequence and treats R as the source of truth.

## Canonical package vignettes

| # | Topic | R vignette | Python companion |
|---|---|---|---|
| 01 | Overview | [Overview](https://ovvo-financial.github.io/NNS/articles/NNSvignette_01_Overview.html) | [`01_overview.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/01_overview.py) |
| 02 | Partial Moments | [Partial Moments](https://ovvo-financial.github.io/NNS/articles/NNSvignette_02_Partial_Moments.html) | [`02_partial_moments.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/02_partial_moments.py) |
| 03 | Correlation and Dependence | [Correlation and Dependence](https://ovvo-financial.github.io/NNS/articles/NNSvignette_03_Correlation_and_Dependence.html) | [`03_correlation_and_dependence.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/03_correlation_and_dependence.py) |
| 04 | Normalization and Rescaling | [Normalization and Rescaling](https://ovvo-financial.github.io/NNS/articles/NNSvignette_04_Normalization_and_Rescaling.html) | [`04_normalization_and_rescaling.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/04_normalization_and_rescaling.py) |
| 05 | Sampling and Simulation | [Sampling and Simulation](https://ovvo-financial.github.io/NNS/articles/NNSvignette_05_Sampling.html) | [`05_sampling_and_simulation.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/05_sampling_and_simulation.py) |
| 06 | Comparing Distributions | [Comparing Distributions](https://ovvo-financial.github.io/NNS/articles/NNSvignette_06_Comparing_Distributions.html) | [`06_comparing_distributions.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/06_comparing_distributions.py) |
| 07 | Clustering and Regression | [Clustering and Regression](https://ovvo-financial.github.io/NNS/articles/NNSvignette_07_Clustering_and_Regression.html) | [`07_clustering_and_regression.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/07_clustering_and_regression.py) |
| 08 | Classification | [Classification](https://ovvo-financial.github.io/NNS/articles/NNSvignette_08_Classification.html) | [`08_classification.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/08_classification.py) |
| 09 | Forecasting | [Forecasting](https://ovvo-financial.github.io/NNS/articles/NNSvignette_09_Forecasting.html) | [`09_forecasting.py`](https://github.com/OVVO-Financial/NNS-python/blob/main/examples/vignettes/09_forecasting.py) |

Canonical changes should be made in the R vignette first. A companion port may
use language-appropriate syntax and containers, but should preserve the same
statistical demonstration and interpretation.

## Applied studies and extended examples

The following material applies NNS to specific problems, comparisons, and
research questions. These examples are useful applications, but they do not
replace the numbered package curriculum.

### Basic statistics

1. [Partial Moment Equivalences](Partial%20Moments%20Equivalences.md)
2. [Bayes' Theorem](Bayes'%20Theorem%20From%20Partial%20Moments.pdf)
3. [CDFs and ANOVA](Continuous_CDFs_and_ANOVA_with_NNS.pdf)
4. [Bias and Confidence Intervals](https://ovvo-financial.github.io/NNS/examples/Bias_and_CI.html)
5. [Partial Moments Estimation Error](https://github.com/OVVO-Financial/Finance/blob/main/Data/Estimation_Error_Replication.md)

### Regression

1. [Overview](https://ssrn.com/abstract=3389938)
2. [Curve Fitting](https://ovvo-financial.github.io/NNS/examples/Curve_Fitting.html)
3. [Nonparametric Regression Using Clusters](http://rdcu.be/tz0J)
4. [Clustering and Curve Fitting By Line Segments](https://ssrn.com/abstract=2861339)
5. [Regression Residuals](https://ovvo-financial.github.io/NNS/examples/Regression_Residuals.html)
6. [Multiple Imputation](NNS_MI_vs_MICE.md)
7. [Logistic Regression Binary Classification](https://ovvo-financial.github.io/NNS/examples/Logistic_Comparison.html)
8. [Boston Housing](https://ovvo-financial.github.io/NNS/examples/Boston_Housing.html)

### Machine learning

1. [Partitional Estimation Using Partial Moments](https://ssrn.com/abstract=3592491)
2. [NNS Regression in Machine Learning](Machine_Learning.pdf)
3. [Classification Using NNS Clustering Analysis](https://ssrn.com/abstract=2864711)
4. [NNS vs. xgboost](https://ovvo-financial.github.io/NNS/examples/xgboost_example.html)
5. [Time-Series Classification](https://ovvo-financial.github.io/NNS/examples/Time_Series_Classification.html)
6. [Time-Series Classification II](https://ovvo-financial.github.io/NNS/examples/Time_Series_Classification_Expanded.html)
7. [Spiral Matching Example](Spiral%20Matching%20Example.pdf)
8. [MNIST](NNS%20vs%20KNN%20MNIST%20dataset.pdf)

### Time-series forecasting

1. [Overview](https://ssrn.com/abstract=3382300)
2. [NNS vs. KERAS](https://ovvo-financial.github.io/NNS/examples/Sunspots_example.html)
3. [NNS vs. prophet](https://ovvo-financial.github.io/NNS/examples/prophet_NNS_comparison.html)
4. [Tides](https://ovvo-financial.github.io/NNS/examples/tides.html)
5. [NNS vs. N-HiTS](NNS.ARMA%20vs%20N-Hits.md)
6. [NNS Time-Series Prediction Interval Benchmark](nns_arma_conformal_benchmark_report.md)

### Econometrics

1. [Econometrics Critiques and Solutions](https://ovvo-financial.github.io/NNS/examples/7_Econometric_Reasons.html)
2. [VAR Alternative](https://ovvo-financial.github.io/NNS/examples/VAR_example.html)
3. [Nowcasting](https://ssrn.com/abstract=3589816)
4. [Causal Analysis](https://ovvo-financial.github.io/NNS/examples/PWT.html)
5. [Federal Reserve Causal Analysis](https://ovvo-financial.github.io/NNS/examples/Causal_Inference_Amongst_Macroeconomic_Variables_Using_NNS.html)
6. [Causal Inference](Causal_Inference_with_NNS_stack.pdf)

## References

The applied examples are demonstrations rather than exhaustive proofs. See the
[book](https://ovvo-financial.github.io/NNS/book/) and the
[research papers](https://papers.ssrn.com/sol3/cf_dev/AbsByAuth.cfm?per_id=1421356)
for the underlying arguments.
