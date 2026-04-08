# Multivariate NNS Inference Workflow

## Case Study: Controlled Effect of Transmission (`am`) on `mpg`

This example illustrates a multivariate inference workflow using the NNS framework to estimate the controlled effect of a binary regressor. The response variable is miles per gallon (`mpg`), and the predictor set consists of transmission type (`am`), weight (`wt`), horsepower (`hp`), and displacement (`disp`). This section also includes a direct comparison to a traditional linear-inference workflow.

```{r setup, message=FALSE, warning=FALSE}
library(NNS)
library(effectsize)
```

```{r data}
data(mtcars)

y <- mtcars$mpg

X <- data.frame(
  am   = mtcars$am,
  wt   = mtcars$wt,
  hp   = mtcars$hp,
  disp = mtcars$disp
)
```

## 1. Directional Relevance Screen

The first step is to establish nonlinear dependence between each predictor and the response prior to fitting the multivariate model.

```{r dependence-screen}
dep_summary <- data.frame(
  variable   = colnames(X),
  dependence = sapply(X, function(x) NNS.dep(x, y)$Dependence)
)

dep_summary <- dep_summary[order(-dep_summary$dependence), ]

print("--- Nonlinear Dependence Screen ---")
dep_summary
```

The nonlinear dependence screen ranks the predictors in descending order of directional relevance to `mpg`. In this specification, `disp` is the strongest predictor, followed by `hp`, `wt`, and then `am`.

## 2. Multivariate NNS Regression Fit

The multivariate regression step estimates the conditional mean surface without imposing parametric linearity.

```{r nns-fit}
fit_full <- NNS.reg(x = X, y = y)

cat("\nModel R2: ", fit_full$R2, "\n")
cat("MAE:  ", mean(abs(fit_full$Fitted.xy$residuals)), "\n")
cat("RMSE: ", sqrt(mean(fit_full$Fitted.xy$residuals^2)), "\n")
```

The model explains approximately 89.8 percent of the variation in `mpg`, with low in-sample absolute and quadratic error. This indicates that the fitted multivariate NNS surface provides a strong nonlinear representation of the data-generating structure.

```{r nns-fit-plot, fig.height=5, fig.width=6}
plot(
  fit_full$Fitted.xy$y,
  fit_full$Fitted.xy$y.hat,
  xlab = "Observed mpg",
  ylab = "Fitted mpg",
  main = "NNS Multivariate Fit: Observed vs Fitted",
  pch = 16,
  col = "steelblue"
)
abline(0, 1, lty = 2, col = "red")
```

## 3. Benchmark Counterfactual at Median Controls

To isolate the effect of transmission, two counterfactual points are constructed: one automatic and one manual, while holding `wt`, `hp`, and `disp` fixed at their median values.

```{r benchmark-counterfactual}
x_ref <- data.frame(
  am   = c(0, 1),
  wt   = median(X$wt),
  hp   = median(X$hp),
  disp = median(X$disp)
)

preds_ref <- NNS.reg(x = X, y = y, point.est = x_ref)$Point.est

controlled_am_effect_at_medians <- preds_ref[2] - preds_ref[1]

cat("\n--- Counterfactual at Medians ---")
cat("\nPredicted mpg (Automatic):", preds_ref[1])
cat("\nPredicted mpg (Manual):   ", preds_ref[2])
cat("\nControlled Effect:        ", controlled_am_effect_at_medians, "\n")
```

At the benchmark covariate profile, the fitted manual-transmission effect is 1.8 miles per gallon relative to automatic transmission.

## 4. Heterogeneous Controlled Effects

Rather than relying on a single benchmark effect, the transmission variable can be flipped for every observed covariate profile in the sample. This produces a heterogeneous controlled-effect distribution.

```{r heterogeneous-effects}
X0 <- X
X0$am <- 0

X1 <- X
X1$am <- 1

pred0 <- NNS.reg(x = X, y = y, point.est = X0)$Point.est
pred1 <- NNS.reg(x = X, y = y, point.est = X1)$Point.est

delta_am <- pred1 - pred0

cat("\n--- Distribution of Controlled am Effects ---")
summary(delta_am)
```

This distribution shows that the controlled transmission effect is heterogeneous across the observed support of the covariates. The average effect is positive, but the full effect profile ranges from negative values to large positive gains.

## 5. Uncertainty Quantification with MEBoot

To quantify uncertainty in the heterogeneous controlled-effect distribution, a maximum-entropy bootstrap is applied to `delta_am`. Confidence intervals are then constructed for the mean and median controlled effects.

```{r meboot}
set.seed(123)
B <- 1000
alpha <- 0.05

boot_delta <- NNS.meboot(x = delta_am, reps = B, rho = 1, drift = FALSE)
boot_delta_mat <- boot_delta["replicates", ]$replicates

boot_mean_effect   <- apply(boot_delta_mat, 2, mean)
boot_median_effect <- apply(boot_delta_mat, 2, median)

mean_ci <- c(
  LPM.VaR(alpha / 2, 0, boot_mean_effect),
  UPM.VaR(alpha / 2, 0, boot_mean_effect)
)

median_ci <- c(
  LPM.VaR(alpha / 2, 0, boot_median_effect),
  UPM.VaR(alpha / 2, 0, boot_median_effect)
)
```

## 6. Results Visualization and Summary

Pointwise interval estimates can be generated for the heterogeneous effect series, and summary intervals can be reported for the mean and median controlled effects.

```{r interval-plot, fig.height=5, fig.width=7}
lower_effect <- apply(boot_delta_mat, 1, function(z) LPM.VaR(alpha / 2, 0, z))
upper_effect <- apply(boot_delta_mat, 1, function(z) UPM.VaR(alpha / 2, 0, z))

plot(
  delta_am,
  ylim = range(c(lower_effect, upper_effect)),
  pch = 16,
  xlab = "Observation Index",
  ylab = "mpg Benefit (Manual - Auto)",
  main = "Heterogeneous am Effect with 95% MEBoot CIs"
)

arrows(
  seq_along(delta_am), lower_effect,
  seq_along(delta_am), upper_effect,
  angle = 90, code = 3, length = 0.03, col = "grey"
)

abline(h = 0, lty = 2)
```

```{r summary-table}
summary_table <- data.frame(
  Statistic = c("Mean Controlled Effect", "Median Controlled Effect"),
  Estimate  = c(mean(delta_am), median(delta_am)),
  Lower_95  = c(mean_ci[1], median_ci[1]),
  Upper_95  = c(mean_ci[2], median_ci[2])
)

cat("\n--- Final Inference Summary ---\n")
summary_table
```

## 7. Traditional Statistical Inference Workflow

For comparison, the same question can be approached using a conventional linear workflow. This highlights the difference between a homogeneous coefficient-based framework and the multivariate NNS conditional-effect framework.

### 7.1 Pearson Correlation Screen

```{r pearson-screen}
cor_summary <- data.frame(
  variable    = colnames(X),
  correlation = sapply(X, function(x) cor(x, y))
)

print("--- Pearson Correlation Screen ---")
cor_summary[order(-abs(cor_summary$correlation)), ]
```

The Pearson screen shows the strongest linear associations with `mpg`. Weight and displacement dominate in absolute magnitude, followed by horsepower and then transmission.

### 7.2 Multivariate OLS Regression

```{r ols-fit}
fit_ols <- lm(mpg ~ am + wt + hp + disp, data = mtcars)

print("--- OLS Model Summary ---")
summary(fit_ols)
```

The OLS fit explains less of the variation in `mpg` than the NNS fit. More importantly, the coefficient on `am` is positive but not conventionally statistically significant at the 5 percent level.

### 7.3 Standardized Effect Sizes

```{r eta-squared}
print("--- Standardized Effect Sizes (Partial Eta-Squared) ---")
eta_squared(fit_ols, partial = TRUE)
```

The partial effect-size view suggests that `am` and `wt` both carry substantial partial explanatory content, even though `am` is not declared statistically significant by the standard t-test.

### 7.4 Controlled Effect and Confidence Interval

```{r ols-coefficient}
am_coeff <- coef(fit_ols)["am"]
am_ci <- confint(fit_ols)["am", ]

cat("\n--- Controlled Effect of Transmission ---")
cat("\nCoefficient (Effect size in mpg):", am_coeff)
cat("\n95% Confidence Interval:        ", am_ci[1], "to", am_ci[2], "\n")
```

The OLS model estimates a positive average transmission effect of about 2.16 mpg, but its confidence interval crosses zero.

### 7.5 ANOVA F-Test for Transmission

```{r anova-test}
fit_reduced <- lm(mpg ~ wt + hp + disp, data = mtcars)

print("--- ANOVA F-test for Transmission Effect ---")
anova(fit_reduced, fit_ols)
```

The nested-model ANOVA reaches the same inferential conclusion as the coefficient test: under the linear model, the added transmission term is not statistically significant at conventional thresholds.

## 8. Discussion: NNS versus Traditional Inference

The contrast between the NNS workflow and the traditional OLS workflow is instructive.

### 8.1 Predictor Hierarchy

The Pearson screen shows that `wt` and `disp` have the strongest linear association with `mpg`, whereas the NNS directional screen ranks `disp`, `hp`, `wt`, and `am`. These are related but not identical views. Pearson correlation measures straight-line association, while `NNS.dep` is built to capture more general nonlinear dependence.

### 8.2 The Linear Significance Trap

In the OLS model, the p-value on `am` is 0.14405. Under conventional inference, that leads to the conclusion that transmission type does not have a statistically significant effect on `mpg` once the controls are included. The NNS workflow reaches a different conclusion because it does not force the controlled effect of `am` into a single global linear coefficient.

This matters when the transmission effect is not constant across the support of `wt`, `hp`, and `disp`. If the effect is positive for some car profiles and weaker or negative for others, a single global coefficient can become unstable and its standard error can expand. The NNS framework instead estimates the conditional surface and then studies the induced heterogeneous effect distribution directly.

### 8.3 Effect Size versus Statistical Significance

The OLS output shows an estimated `am` effect of 2.159271 mpg, which is directionally consistent with the NNS mean controlled effect of 2.526736 mpg and the benchmark controlled effect of 1.8 mpg. The difference is interpretive:

- OLS treats the effect as a homogeneous coefficient with one confidence interval.
- NNS treats the effect as a heterogeneous conditional contrast that varies by covariate profile.

The OLS confidence interval for `am` crosses zero, whereas the NNS summary intervals for the mean and median controlled effects are both strictly positive.

### 8.4 Model Fit

The NNS model achieves:

- `R2 = 0.8978988`
- `MAE = 1.011309`
- `RMSE = 1.963105`

The OLS model achieves:

- `R2 = 0.8402`
- residual standard error `= 2.581`

So even at the level of fit, the multivariate NNS model provides a stronger representation of the observed response surface.

### 8.5 Controlled Effect Interpretation

The NNS benchmark counterfactual implies that, at the median values of `wt`, `hp`, and `disp`, manual transmission is associated with a fitted increase of 1.8 mpg. Across the observed covariate support, the effect distribution is heterogeneous:

- mean effect = `2.526736`
- median effect = `1.800000`
- minimum effect = `-3.400`
- maximum effect = `8.900`

This reveals something the OLS coefficient cannot: the controlled effect of transmission is not constant across all car profiles.

## 9. Comparative Summary

| Feature | Multivariate NNS Results | Traditional OLS Results |
|---|---:|---:|
| Model Fit | `R2 = 0.8978988` | `R2 = 0.8402` |
| Error Structure | `MAE = 1.011309`, `RMSE = 1.963105` | Residual SE = `2.581` |
| Transmission Effect | `+1.8` at medians | `+2.159271` average coefficient |
| Heterogeneity | Explicitly modeled | Not modeled |
| Mean Effect Interval | `[2.322152, 2.731320]` | Coefficient CI crosses zero |
| Median Effect Interval | `[1.559981, 3.444243]` | Not available as an effect-distribution summary |
| Inference | Conditional, nonlinear, heterogeneous | Linear, homogeneous, coefficient-based |

## 10. Conclusion

This example shows that the effect of transmission on fuel efficiency can be analyzed within a fully multivariate, nonlinear, nonparametric framework. The multivariate `NNS.reg` fit indicates strong predictive structure in the data, while the counterfactual contrast isolates the conditional effect of `am`. The resulting controlled-effect distribution shows that the transmission effect is positive on average, but heterogeneous across the observed support of the controls. The MEBoot intervals reinforce that the mean and median controlled effects remain positive.

The traditional OLS workflow gives a positive point estimate for transmission, but because it forces the effect into a single linear coefficient, it fails to characterize the heterogeneity revealed by the NNS analysis and does not declare the effect statistically significant at conventional levels.

Methodologically, the key distinction is that the NNS workflow answers a conditional and heterogeneous effect question, while the traditional workflow answers a homogeneous coefficient question. Although `mtcars` is only a toy dataset, the inferential structure of the NNS workflow is intended to generalize: nonlinear relevance screening, multivariate conditional mean estimation, fitted counterfactual contrasts for binary regressors, heterogeneous effect analysis, and bootstrap interval estimation can all be carried forward to richer empirical settings.
