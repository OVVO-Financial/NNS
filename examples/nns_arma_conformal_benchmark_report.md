# NNS time-series prediction interval benchmark: corrected protocol

The earlier numerical ranking on this page is withdrawn. Its NNS call included
scored outcomes in the optimizer's internal validation data, while the other
methods used observed one-step lags and NNS forecast a longer chunk from one
origin. Those scores cannot support a comparative performance or coverage claim.
The original material remains available in Git history at commit
`1075a338d932db13a8210e8ed23fb94f70a3f21f`.

## Reproducible evaluation

The executable source is [nns_arma_interval_benchmark.R](nns_arma_interval_benchmark.R).
It retains the simulated nonlinear heteroskedastic process and the existing
comparison methods, with these corrections:

1. All methods issue one-step forecasts from the same chronological origins.
   `MAX_H` is one and any different value is rejected.
2. NNS receives only the observations available at the forecast origin and an
   explicit `h = 1`. Seasonal candidates are selected from the earlier internal
   training prefix. The remaining observed prefix is internal validation; the
   next outcome is scored only after forecasting.
3. Baseline estimation, calibration, and scored test periods remain separate.
   Online calibration methods may use a test outcome only for later forecasts.
4. Native NNS interval endpoints are scored for coverage, width, and interval
   score. A lower endpoint, point forecast, and upper endpoint do not specify a
   predictive distribution. CRPS and logscore are therefore `NA` for NNS, as for
   other interval-only methods. A degree-two partial-moment transform of three
   supplied points would define an additional modelling choice, not a native
   predictive density or a validated probability calibration.
5. Each run records the selected periods, method, weights, shift, forecast
   origin, and horizon, together with raw per-seed scores. Keep package versions,
   `sessionInfo()`, seeds, and optional dependency availability with results.

```r
source("examples/nns_arma_interval_benchmark.R")
result <- run_all()
sessionInfo()
```

The full five-seed, 3,500-point experiment now entails thousands of honest
one-step optimizer fits and can take substantial time. The corrected source is
provided for that experiment; this page does not substitute a small smoke run
for a full benchmark or carry forward the invalid historical numbers.

## What can be concluded

The comparison estimates performance for this process, seed set, horizon,
training policy, and loss. Report uncertainty across series or origins and
inspect coverage by volatility regime. The true conditional oracle is a
reference with privileged information, not an available forecasting method.
Conformal finite-sample guarantees require the applicable exchangeability or
other calibration assumptions; this serially dependent, changing process does
not satisfy ordinary iid split-conformal assumptions automatically. NNS residual
bands also require empirical out-of-sample validation.

The baseline uses `glmnet` ridge when installed and `lm` otherwise; these are
separate configurations and must be identified in results. Model refitting
policies are likewise part of the compared procedure: NNS updates its fit at
origins while the original lagged baseline is fitted on the initial training
window. A claim about estimator quality alone requires an additional matched
refitting-policy comparison.

See [example corrections](example-corrections.md) for the related distinction
between frequency calibration and directional severity thresholds.
