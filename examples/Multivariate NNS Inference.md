# Multivariate NNS fitted contrasts and resampling uncertainty

This example uses `mtcars` to compare fitted miles per gallon at automatic and
manual transmission settings, holding the other predictors fixed. NNS preserves
the nonlinear conditional response surface and allows the contrast to vary
across covariate profiles. These are **model-based conditional associations**.
Interpreting them causally additionally requires a suitable design,
exchangeability conditional on the controls, overlap, and a defensible outcome
model; this small observational dataset does not establish those conditions.

## Fit and describe the response surface

```r
library(NNS)
X <- mtcars[c("am", "wt", "hp", "disp")]
y <- mtcars$mpg

fit <- NNS.reg(X, y, ncores = 1)
data.frame(
  squared_fitted_correlation = fit$R2,
  training_MAE = mean(abs(fit$Fitted.xy$residuals)),
  training_RMSE = sqrt(mean(fit$Fitted.xy$residuals^2))
)

conditional_contrasts <- function(train_X, train_y, reference_X) {
  X0 <- X1 <- reference_X
  X0$am <- 0
  X1$am <- 1
  prediction <- NNS.reg(train_X, train_y, point.est = rbind(X0, X1),
                        ncores = 1)$Point.est
  n <- nrow(reference_X)
  as.numeric(prediction[n + seq_len(n)] - prediction[seq_len(n)])
}

delta <- conditional_contrasts(X, y, X)
summary(delta)
plot(delta, xlab = "Observed covariate profile", ylab = "Fitted manual - automatic mpg")
abline(h = 0, lty = 2)
```

`fit$R2` is squared observed/fitted correlation. It is not automatically the
residual-sum-of-squares definition used by `summary(lm(...))`. Training error
describes the fitted sample; predictive comparisons require a common held-out
sample and the same loss. The spread of `delta` describes **heterogeneity across
profiles**, not sampling uncertainty of the fitted surface.

## Refit when estimating model uncertainty

Resampling the already computed `delta` vector holds the fitted NNS surface
fixed. Such resampling can describe variation of summaries of that fixed
profile distribution, but excludes uncertainty from fitting the model.
The following pairs bootstrap refits NNS in every replicate and evaluates each
refit at the same reference profiles. Its target is the average or median fitted
contrast over this fixed empirical reference population.

```r
set.seed(20261103)
B <- 199L                       # increase for a substantive analysis
bootstrap_delta <- vapply(seq_len(B), function(b) {
  rows <- sample.int(nrow(X), replace = TRUE)
  conditional_contrasts(X[rows, , drop = FALSE], y[rows], X)
}, numeric(nrow(X)))

alpha <- 0.05
bootstrap_summary <- rbind(
  mean = colMeans(bootstrap_delta),
  median = apply(bootstrap_delta, 2, median)
)
intervals <- t(apply(bootstrap_summary, 1, quantile,
                     probs = c(alpha / 2, 1 - alpha / 2), type = 7))
data.frame(estimate = c(mean(delta), median(delta)), intervals,
           row.names = c("Mean fitted contrast", "Median fitted contrast"))
```

The resampling unit and calibration must match the study design. Independent
pairs resampling assumes independent sampling units; clustered or time-dependent
data need an appropriate design. These percentile intervals are approximate and
do not guarantee nominal coverage for small samples or an adaptive estimator.
If predictor selection or tuning is part of the analysis, repeat it inside
every resample. To target a changing population of profiles, resample or
otherwise model that reference population as well.

## Compare like targets

```r
linear_fit <- lm(mpg ~ am + wt + hp + disp, data = mtcars)
linear_contrast <- unname(coef(linear_fit)["am"])
data.frame(model = c("NNS", "Additive linear model"),
           mean_fitted_contrast = c(mean(delta), linear_contrast))
```

The additive linear model constrains the transmission contrast to be constant;
the NNS fit permits heterogeneity. A linear coefficient interval and an interval
obtained by resampling fixed NNS contrasts measure different uncertainty and
cannot establish that one method detects an effect the other misses. Compare
common contrasts with full refitting and common held-out predictive losses.
The sign and size of this example's results should be recomputed for the package
version in use; no historical fit or interval is asserted as a current result.

See [the example corrections and reproducibility notes](example-corrections.md)
and the book's regression chapter for the estimator's intended partial-moment
construction.
