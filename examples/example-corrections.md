# Example corrections and reproducibility

These corrections accompany the November 2026 implementation and book review.
The numbered R vignettes remain the canonical curriculum. Positive-degree
partial moments retain their directional severity meaning; event probabilities,
uncertainty, and predictive performance require their corresponding definitions
and evidence.

## Corrections to applied examples

| Topic | Correct interpretation or revised procedure |
|---|---|
| Multivariate fitted contrasts | Describe conditional associations and heterogeneity. Resample sampling units and refit the model to include fitting uncertainty; do not compare a fixed-fit contrast bootstrap to an OLS coefficient interval as if they had the same target. The editable [workflow](Multivariate%20NNS%20Inference.md) now refits. |
| Time-series interval comparison | Tune on past data and score an explicit future horizon. Match forecast origins and available lags. The invalid historical ranking is withdrawn; the [corrected report](nns_arma_conformal_benchmark_report.md) links executable source. Three interval summary points do not identify a native predictive density. |
| MNIST and machine-learning claims | The supplied PDF's NNS score is one reported run; it does not include a matched kNN fit establishing a 75% error reduction. Those comparative and general dimensionality claims are withdrawn pending common splits, preprocessing, tuning, metrics, and repeated runs. |
| Five-way ranking | Normalized degree-one UPMs are shares of above-target deviation severity. They are not probabilities of winning a race. Sampling winners using those same shares merely reproduces the supplied weights; it does not validate their event-frequency interpretation. |
| Partitional estimation | For `rnorm(n) + runif(n)`, mean noise is 0.5. The conditional mean is `f(x) + 0.5`, or use `runif(n) - 0.5` to target `f(x)`. Compare all fitted methods on a common held-out set; random-forest OOB error and NNS training error are different assessments. |
| Dependence permutation output | Shuffled-statistic quantiles are a null reference distribution under exchangeability, not a confidence interval for population dependence. Specify dependence-preserving resampling for serial or clustered data. |
| Copula sampling | Degree-zero marginal transforms or an explicitly randomized probability transform at atoms target frequency. Positive-degree inverses change the marginal law; their directional geometry does not imply preservation of the original marginal distribution. |
| Iris classification | Rows 141:150 are all virginica. Vignette08 now uses a seeded stratified split, removes historical perfect-accuracy output, and makes class decoding explicit. |
| Derivative summaries | Quantiles across fitted derivatives describe heterogeneity over evaluation locations. Sampling uncertainty requires resampling and refitting with a design appropriate to the observations. Smoothness and admissible local steps still matter. |
| Normalization | Positive constant multiplication of columns preserves Pearson correlation; a negative multiplier reverses its sign. Changing portfolio weights is a separate operation. Check scaling-factor signs before interpreting transformed associations. |
| Forecast period selection | Periods selected using a test window make that window validation. Vignette09 now selects periods within the training prefix, reserves the final observations for testing, and computes current outputs rather than printing stale values. |
| Noisy differentiation | Adding frozen real noise to a complex-valued evaluation leaves the imaginary part unchanged. Tiny complex-step error in that construction is not evidence of generic noise robustness. A reproducible noise generator, complex-input contract, and repeated evaluation protocol are required. |

## Published PDF source status

The twelve published PDFs are retained as historical documents. Matching original
editable sources were not present in the audited repository. The new
[revised-companions.Rmd](revised-companions.Rmd) is an explicitly revised,
executable companion, **not a reconstruction of the original authors' source**.
It provides corrected demonstrations and clear limits for each PDF below.

| Historical PDF | Revised companion section |
|---|---|
| `Bayes' Theorem From Partial Moments.pdf` | Conditional event probabilities |
| `Continuous_CDFs_and_ANOVA_with_NNS.pdf` | Frequency, severity, and comparison targets |
| `Nonlinear_Correlation_and_Dependence_Using_NNS.pdf` | Directional dependence and null calibration |
| `Normalization.pdf` | Positive scaling and association |
| `5_way_rank_probabilities_using_R.pdf` | Severity shares and winner frequencies |
| `Partitional_Estimation_Using_Partial_Moments.pdf` | A correctly centered regression target |
| `Machine_Learning.pdf` | Held-out classification |
| `NNS vs KNN MNIST dataset.pdf` | A reproducible MNIST comparison contract |
| `Spiral Matching Example.pdf` | Geometric matching |
| `Nonlinear_FWL.pdf` | Residualization and conditional targets |
| `Causal_Inference_with_NNS_stack.pdf` | Predictive contrasts and identification |
| `NNS ARMA.pdf` | Training, validation, and test chronology |

The [rendered companion](https://ovvo-financial.github.io/NNS/examples/revised-companions.html) is rebuilt by the documentation workflow. Build it locally with:

```sh
Rscript examples/build-revised-companions.R
```

The build records `sessionInfo()` and the NNS version. Original source recovery
remains a provenance task: when those sources are available, add them with their
assets and build instructions rather than representing revised material as an
exact source restoration. Historical numerical tables are not carried forward
as current validation results.
