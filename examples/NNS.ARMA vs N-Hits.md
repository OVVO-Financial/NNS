# 🚀 NNS.ARMA.optim vs. n-HiTS: A Surprising Forecasting Result

When it comes to time-series forecasting, deep learning models like **N-BEATS** and **n-HiTS** often dominate the conversation. But what happens when we put them head-to-head against a nonparametric statistical approach?

I ran a comparison using hourly traffic volume data  
**(original article here: https://www.datasciencewithmarco.com/blog/all-about-n-hits-the-latest-breakthrough-in-time-series-forecasting)**

## 📊 Results (MAE)

| Model                    | MAE     |
|--------------------------|---------|
| Baseline                 | 249.0   |
| N-BEATS                  | 292.0   |
| N-BEATS + covariates     | 288.0   |
| n-HiTS                   | 266.0   |
| **NNS.ARMA.optim()**     | **236.17** |

**Yes — the nonlinear nonparametric `NNS.ARMA.optim()` outperformed all of them**, including the deep learning-based n-HiTS.

## ✅ Exact Reproducible R Code

```r
library(NNS)

# Read Data
daily_traffic <- read.csv("https://raw.githubusercontent.com/marcopeix/time-series-analysis/refs/heads/master/data/daily_traffic.csv")

# Create train / test sets (last 120 observations = test set)
train_set <- head(daily_traffic$traffic_volume, length(daily_traffic$traffic_volume) - 120)
test_set  <- tail(daily_traffic$traffic_volume, 120)

# Determine seasonal periods (modulo = 24 because data is hourly)
periods <- NNS.seas(train_set, modulo = 24)$periods

# Optimize seasonal periods + ARMA parameters using MAE as objective
nns_estimates <- NNS.ARMA.optim(train_set, 
                                h = 120, 
                                seasonal.factor = periods,
                                obj.fn = expression(Metrics::mae(actual, predicted)), 
                                objective = "min",
                                plot = TRUE, 
                                negative.values = FALSE)

# Final MAE on test set
Metrics::mae(nns_estimates$results, test_set)

# Plot actual vs. forecast
plot(test_set, 
     col = "blue", type = "l", lwd = 2,
     main = "NNS.ARMA.optim() Forecast", 
     ylab = "traffic_volume", xlab = "Index")
lines(nns_estimates$results, col = "red", lwd = 2)
legend("topleft", legend = c("Actual", "NNS.ARMA.optim() Forecast"), 
       col = c("blue", "red"), lwd = 2, bty = "n")
```

<img src="/examples/nhits_2.png"  style="border: none; outline: none; margin: 0; padding: 0; display: block;"/>

```
Console output from the run (exact values from the screenshot):
textMetrics::mae(nns_estimates$results, test_set)
[1] 236.17
```


✅ Takeaway
Sometimes simplicity + interpretability beats complexity.
Before jumping into the latest neural architecture, it’s worth asking:
Can a nonparametric approach solve the problem faster, with fewer resources, and better performance?

NNS on CRAN: [Install the latest version](https://cran.r-project.org/package=NNS)
