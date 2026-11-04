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
| Independent 3-D copula output | Page 9 of `Nonlinear_Correlation_and_Dependence_Using_NNS.pdf` prints `0.09571785` for `set.seed(123); x <- rnorm(1000); y <- rnorm(1000); z <- rnorm(1000); NNS.copula(cbind(x,y,z))`. The reviewed implementation returns `0.197845241162`; the revised companion recomputes it. This descriptive finite-sample score is not a p-value or an independence threshold. The population event mass `2 * (1/2)^3` at independent marginal medians is a different quantity. |
| Copula sampling | Degree-zero marginal transforms or an explicitly randomized probability transform at atoms target frequency. Positive-degree inverses change the marginal law; their directional geometry does not imply preservation of the original marginal distribution. |
| Iris classification | Rows 141:150 are all virginica. Vignette08 now uses a seeded stratified split, removes historical perfect-accuracy output, and makes class decoding explicit. |
| Derivative summaries | Quantiles across fitted derivatives describe heterogeneity over evaluation locations. Sampling uncertainty requires resampling and refitting with a design appropriate to the observations. Smoothness and admissible local steps still matter. |
| Normalization | Positive constant multiplication of columns preserves Pearson correlation; a negative multiplier reverses its sign. Changing portfolio weights is a separate operation. Check scaling-factor signs before interpreting transformed associations. |
| Forecast period selection | Periods selected using a test window make that window validation. Vignette09 now selects periods within the training prefix, reserves the final observations for testing, and computes current outputs rather than printing stale values. |
| Noisy differentiation | Adding frozen real noise to a complex-valued evaluation leaves the imaginary part unchanged. Tiny complex-step error in that construction is not evidence of generic noise robustness. A reproducible noise generator, complex-input contract, and repeated evaluation protocol are required. |

## Historical HTML output regeneration

The six archived HTML pages listed in issue #72 have been rebuilt with the
current NNS from R Markdown sources in [`sources/`](sources/), or retired. Each
source records its data provenance and its training, validation and test
windows, prints `sessionInfo()`, and rebuilds with
`Rscript examples/sources/build.R <name>`. Scores from the archived pages are
not carried forward.

- [x] [Macroeconomic causal analysis](Causal_Inference_Amongst_Macroeconomic_Variables_Using_NNS.html): `NNS.nowcast()` is no longer in NNS. The same 29 FRED series are aligned monthly with `NNS.VAR(h = 0)`, which makes no forecast and so runs no validation step. `NNS.caus` is reported as a one-month directional lagged association between transformed (detrended) series, with circular-shift p-values and a Holm adjustment. The archived conclusion that Fed assets cause the listed variables is withdrawn.
- [x] [Archived overview](NNSvignette_Overview.html): retired. The page forwards to the canonical numbered vignettes.
- [x] [Sunspots](Sunspots_example.html): each slice selects its periods with `NNS.seas` and `NNS.ARMA.optim` on its own training window, before its 120-month evaluation horizon. The KERAS LSTM scores are quoted from the cited article and labelled as not rerun.
- [x] [VAR example](VAR_example.html): the full pipeline (VAR, `NNS.ARMA`, auto ARIMA, `NNS.boost` and `NNS.stack`, random forest, `NNS.VAR`) is rerun on FRED data for the archived 1954 Q4 to 2007 Q4 window. Internal-validation scores and the 12-quarter test scores are reported in separate, labelled tables. The Kolmogorov-Smirnov comparisons are withdrawn.
- [x] [Prophet comparison](prophet_NNS_comparison.html): both methods use the same data, forecast origins and horizons, and neither is tuned on a test window. The Peyton Manning periods come from the training data, the lynx forecast has the test window's length, and the S&P 500 section uses FRED's `SP500` series at three fixed origins. A last-value forecast is scored for reference, and the archived tallies are withdrawn.
- [x] [Tides](tides.html): rerun on the Win-Vector data file, checked by MD5, with the training and test windows printed. The forecast origin precedes the test window, and the archived R-squared values are not carried forward.

## Published PDF source status

The twelve published PDFs are retained as historical documents. Their original
editable sources were not recovered. [`sources/`](sources/) now holds an R
Markdown source for each PDF, reconstructed in November 2026 from the published
text (OVVO Labs' Markdown conversions of the PDFs). Each source runs against the
current NNS and notes, at the top and in "Note (2026)" boxes, where its output
or interpretation differs from the published PDF. Rebuild them with
`Rscript examples/sources/build.R` (all sources) or name one or more, for
example `Rscript examples/sources/build.R Normalization`. PDFs need a LaTeX
installation; without one, or with `--html`, they render to HTML. Re-rendering
replaces a published PDF, so the PDFs below have not been regenerated.

The [revised-companions.Rmd](revised-companions.Rmd) remains an explicitly
revised, executable companion. It provides corrected demonstrations and clear
limits for each PDF below.

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
