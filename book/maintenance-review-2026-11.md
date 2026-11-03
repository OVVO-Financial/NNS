# NNS maintenance review — November 2026

This record reconciles the 20 issues and 18 pull requests reviewed together on
3 November 2026. The baseline was `1075a338d932db13a8210e8ed23fb94f70a3f21f`
on `NNS-Beta-Version`. The implementation preserves the book's partial-moment
framework, directional measures, gravity reducer, recursive partitions, and
two-stage multivariate forecasting. Corrections distinguish an algorithm's
defined target from stronger inference or performance claims.

## Pull requests

The following original PRs were reviewed and merged individually:

| PR | Repair |
| --- | --- |
| #132 | Mixed-degree co-partial moments |
| #125 | Scale-invariant integer-degree VaR root stopping |
| #141 | Two-sided inclusive correlation permutation tails |
| #134 | Paired-length validation before native parallel work |
| #123 | Target hazard/cumulative-hazard evaluation |
| #126 | Stochastic-dominance cluster ordering |
| #127 | Residual-cluster handling |
| #147 | Repeated small-sample derivative bandwidths |
| #135 | Training-only ARMA support selection |
| #146 | Monte Carlo exponent grids within the requested bounds |
| #148 | Superiority intervals with constant or singleton groups |
| #149 | Bootstrap handling of exactly zero residuals |
| #133 | Scalar/bulk distance exponential-rate agreement |
| #143 | Degenerate copula targets |

The completion branch integrates #130, #138, #139, and #140 together. Their
interactions require joint verification: gravity changes affect forecast
fixtures; nonlinear shrinkage affects optimizer selection; returned objectives
must score the selected forecast; and a single shrunk lag must follow the same
recursive rule as a multivariate lag combination. Original PR heads remain in
the completion commit's ancestry.

## Issue dispositions

“Resolved” below describes the implemented or explicitly documented decision.
It does not assert universal statistical calibration or a new benchmark win.

| Issue | Disposition and evidence in the repository |
| --- | --- |
| #65 | **Partly resolved; open.** Partial derivatives use the actual selected predictor coordinates, observed evaluation tuples, bounded local probability windows, and explicit one-sided/nonuniform differences. Binary predictors return finite contrasts; constant coordinates have undefined derivatives. A one-column data frame of evaluation points is accepted like a numeric vector. The estimates are still attenuated toward zero, and the attenuation does not shrink with the sample size: for noiseless `y = x1 * x2` with `x ~ U(-1, 1)^2` (`set.seed(1)`), `dy/dx1` at (0.5, 0.5) is 0.24 at n = 500 and 0.11 at n = 2000 against 0.5, and the mixed partial is 0.31 and 0.22 against 1. The issue remains open for the derivative estimator. See `R/dy_d_wrt.R` and local-support/Stack geometry tests. |
| #72 | **Open for six historical HTML outputs.** The exact book/vignette text corrections, current executable vignette sources, regenerated Chapter 25 sine figure, and independent 3-D copula PDF erratum are addressed. The site workflow builds the package reference, full book HTML, and revised companions, but copies six archived HTML analyses without rerunning them. Their staged copies carry notices, and `examples/example-corrections.md` lists the bounded regeneration checklist. A site build is not a full historical recomputation. Original PDF-source provenance remains #118. |
| #114 | **Resolved.** Mode detection retains discrete data values and edge peaks; regression's continuous mode reducer no longer rounds observations to class codes. See the mode regression tests. |
| #115 | **Resolved.** General one-sided bounds use Cantelli or directional Markov inequalities. The halved Chebyshev expression is limited to its symmetry assumptions. |
| #116 | **Resolved.** Degree-zero frequencies, positive-degree deviation-area shares, and the induced degree-one CDF/expectiles are distinguished. The mean anchor is retained. Type-7 interpolation and future-coverage limits are explicit. |
| #117 | **Resolved.** Gravity uses corrected quartiles and observation-centered peaks while retaining its original location blend. The unbounded mean component and lack of a positive breakdown point are documented. |
| #118 | **Open: original editable sources not recovered.** The 12 historical PDFs are inventoried in `examples/example-corrections.md`. `examples/revised-companions.Rmd` supplies new editable, executable companions and a build script; it is explicitly not represented as the original source. Recovering the authors' original Rmd/TeX/Word files remains necessary for faithful regeneration of those PDFs. |
| #119 | **Resolved.** Chapters 1–14 correct the paired copula surface, moment targets, independence products, feasible moment examples, atom/PIT qualifications, hinge recovery, and directional-causation interpretation. |
| #120 | **Resolved.** Chapters 15–28 distinguish classification indicators from numeric labels, joint RPM cells from one-dimensional bounds, assumptions sufficient for cell-mean consistency from bounded defaults, and the actual regression/forecast algorithms from unsupported general guarantees. |
| #121 | **Resolved with explicit estimands.** Mean, median, and distribution-profile agreement remain distinct. Median jumps and constant-group area balance are handled explicitly. Certainty is descriptive; a pooled permutation reference is not automatically a heteroscedastic mean test. Separate-endpoint bootstrap intervals retain their conservative construction. |
| #122 | **Resolved by accurate labeling.** Continuous degree-one dominance retains its numerical criterion and uses `OMEGA` labels. Discrete FSD and numerical dominance/set/cluster interfaces retain their meanings. |
| #124 | **Resolved.** Nominal classification combines class-indicator probabilities before selecting labels. Tests cover label permutations, numeric class values, factors, Stack, and Boost. Prediction sets are descriptive, without a calibration guarantee. |
| #128 | **Resolved with a documented boundary policy.** Exact categorical cells, factor-level invariance, reflection-invariant complexity selection, and Stack's actual multivariate geometry are repaired. Squared correlation remains the documented R-squared statistic. Existing endpoint extrapolation is retained and its sensitivity to terminal spacing/noise is explicit; no new local-linear estimator is silently substituted. |
| #129 | **Partly implemented; open for selection policy.** `dispersion = "absolute"` provides translation-invariant seasonal scoring while `"auto"` preserves historical behavior. The default scan, explicit `modulo`, predictive-lag target, and permutation assumptions are documented. Fundamental-period selection, a different minimum-cycle rule, and a default significance-filter policy require an explicit design and validation decision; no uncalibrated divisor heuristic was made the default. |
| #131 | **Resolved.** VAR validation reconstructs prefix completion, first-stage forecasts, predicted contemporaneous/later-lag features, and relevance selection at a chronological origin. Only observed holdout targets are scored. All forecasting fits explicitly use numeric regression, including two-valued responses. Final fitting uses full available history. `validation = "observed"` remains an explicit compatibility option. |
| #136 | **Resolved as a documented scenario/inference distinction.** Residual expansion has its own switch. Forced-mean dispersion is identified as an iid reference; zero rank correlation does not establish independence, and interior mixing need not preserve margins or seasonal structure. No universal serial-dependence sampling calibration is claimed. |
| #137 | **Resolved through corrected sources and explicit historical errata.** Examples use refitted sampling uncertainty, reserved forecast outcomes, indicator classification, correctly described permutation references, and matched statistical targets. Unsupported rankings and numerical tables are withdrawn. Revised PDF companions identify limitations and do not invent replacement benchmark results. The `NNS.norm` demonstrations in vignettes 01 and 04 use positive-mean columns, the `NNS.boost` and `NNS.stack` help examples and vignette 08 use a stratified iris test set (rows 141:150 alone are all virginica), and vignette 08 describes `NNS.boost` as feature screening followed by one stacked fit. |
| #142 | **Resolved: no shift by default.** `bias.correction = "none"` is the default. Neither acceptance rule made the shift carry over to held-out forecasts: on rolling-origin forecasts of nine series (six expanding origins, h = 12; AirPassengers, UKgas, LakeHuron, a random walk, sunspot.month from 1974, nottem, USAccDeaths, co2 and Nile) the uncorrected forecasts had the lower mean MASE on eight, against both `"split"` (estimate on the first half of the validation window, accept on the second) and the former rule (estimate and accept on the same window). `"split"` remains an option; the same-window rule is removed. The returned raw/corrected forecasts and validation metadata expose the policy. Whole-model performance still requires independent outer evaluation. |
| #144 | **Resolved by scope and null clarification, with #141/#143 numerical fixes.** Directional lag association is not structural causal identification or own-history conditioning. Finite partition resolution, fixed-target copula blind spots, circular-shift assumptions, and exact/Monte Carlo permutation budgets are documented. |
| #145 | **Resolved by the supported-interface decision.** Multivariate `type = "hazard"` rejects the undefined estimator instead of reporting CDF/survival as density/survival. Cumulative hazard retains its stated transform. Derivative documentation identifies noisy second differences and the observation-weighted meaning of `"overall"`; it does not silently switch to a different weighting target. |

## Compatibility decisions

- `NNS.FSD(type = "continuous")` returns `X OMEGA Y`, `Y OMEGA X`, or
  `NO OMEGA DOMINANCE EXISTS`. Code comparing the old continuous strings must
  update.
- Classification exposes probabilities, descriptive prediction sets, and
  `class.levels` for factor/character output codes. Nominal Stack models report
  a selected dimension count and class-specific equations; a scalar coefficient
  threshold is not meaningful for that model and is `NA`.
- `NNS.reg` and `NNS.stack` accept `type = "REG"` for explicit numeric regression, including two-valued responses.
- `NNS.ARMA.optim` defaults to `bias.correction = "none"`, so it no longer
  adds a bias shift; `"split"` remains available, and the former same-window
  rule is removed. `NNS.VAR`
  defaults to `validation = "pipeline"`, whose validation uses no later data;
  on rolling-origin forecasts of ten panels it matched `"observed"` (mean
  scaled ensemble error 1.433 for both), so it is not an accuracy claim. Their
  compatibility settings are explicit.
- `NNS.seas` keeps `dispersion = "auto"`; `"absolute"` is opt-in. Bootstrap
  residual expansion also retains its prior default.
- Geometry, gravity, nominal-model, and validation fixes can change predictions.
  A single shrunk ARMA lag now incurs the recursive fitting cost already present
  for multiple lags.

## Verification and reproducibility

The completion PR is the authoritative record of final CI results. Its checks
cover the full installed-package tests on Windows, Linux, and macOS, the native
sanitizer workflow, and the combined pkgdown/book/example build.

Focused local checks cover formula invariants, categorical-label permutations,
response reflection, scalar/all-k/single-k/bulk distance parity, degenerate
ANOVA groups, derivative support boundaries, prefix-only forecast validation,
returned-model objectives, and serial/fork/PSOCK agreement. The corrected
inference example executes 199 model refits. A forecast-origin sentinel check
changes future outcomes and confirms that the earlier forecast is unchanged.

The historical N-Hits, MNIST, and multi-seed conformal rankings were not rerun
under a matched protocol. They are not claimed as validated results. External
data examples and historical PDF figures are distinguished from newly executed
examples and regenerated figures. The original PDF files remain archived;
their corrections and new editable companions are linked from the example index.

## Review follow-up

Changes made after the review of the completion branch:

- `NNS.stack` method 1 keeps the same predictor coordinates in every fold,
  ranked by the dimension-reduction coefficients, instead of the columns with
  a nonzero coefficient in that fold. A predictor that was constant in one
  chronological training fold gave that fold a univariate design, and the
  stack stopped because no neighbour count could be tuned there.
- A factor level that is declared but has no training rows is a valid test
  value in `NNS.reg`, `NNS.stack` and `NNS.boost`, as it was before the
  completion branch; an undeclared character value is still rejected.
  `NNS.stack` also accepts a single test row with a factor predictor.
- Nominal classification computes each fold's cells, class shares and sorted
  distances once for every neighbour count, skips fitted classes in
  point-only calls, and shares the rank weights across points. Outputs are
  unchanged. For a three-class problem with two predictors and `folds = 2`,
  `NNS.stack` takes 1.3 s and `NNS.boost` (3 learner trials, 2 epochs) 1.3 s
  at n = 1000, against about 20 s each before; at n = 5000 they take 5.8 s
  and 8.0 s.
- Categorical `NNS.reg` honours an explicit `n.best`; the exact-cell rule
  applies with the default `n.best = NULL`.
- `dy.d_` accepts a one-column data frame of evaluation points, and the bulk
  distance helpers match an unnamed test matrix by position.
- The classification `NNS.stack` passes `tau = "ts"` to its dimension
  reduction when `ts.test` is set, as the regression path does. `Uni.caus`
  scores aligned pairs for both values (#96), so the coefficients are
  unchanged.
- Documentation: degree-zero `LPM.VaR` is the type-7 quantile and need not
  satisfy `F0(t) = p`; degree-one smoothness is not robustness (chapters 15
  and 17); regression points compress and average only where cells hold
  several observations (chapter 22); `NNS.boost` epochs re-test the surviving
  subsets on a balanced schedule and `dist` has no DTW option (chapter 24);
  the chapter 23 workflow summary; the full multi-level `NNS.dep` definition
  (chapter 10); partitions split at gravity, not the mean (`NNS.dep` help,
  chapter 14); and `OBJfn.*` values are selection-optimized scores.

Decisions:

- The multivariate estimator is not changed to weight regression points by
  their case counts. Balancing duplicates rows, and in fine-cell multivariate
  fits duplicate rows share one regression point, so balancing has little
  effect there; `NNS.boost` fits its final stack with `balance = FALSE`.
- `NNS.copula` returns 1 when every observation is concordant. That is perfect
  concordance, the measure's maximum, and the correct value; short windows
  reach it more often only because a few observations can all be concordant
  (#144).
- The same-window bias rule is removed rather than kept as a compatibility
  option.
