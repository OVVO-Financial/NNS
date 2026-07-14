from pathlib import Path

path = Path("R/Regression.R")
text = path.read_text()
old = '''.nns_reg_r2 <- function(actual, predicted) {
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)
  if (sst == 0) return(if (sse == 0) 1 else 0)
  1 - sse / sst
}
'''
new = '''.nns_reg_r2 <- function(actual, predicted) {
  # Racine-Hastie nonparametric goodness-of-fit: squared sample correlation
  # between observed and fitted values. Unlike predictive 1 - SSE / SST, this
  # within-sample descriptive measure is bounded in [0, 1].
  actual.centered <- actual - mean(actual)
  predicted.centered <- predicted - mean(predicted)
  denominator <- sum(actual.centered^2) * sum(predicted.centered^2)

  if (denominator == 0) {
    return(if (isTRUE(all.equal(actual, predicted))) 1 else 0)
  }

  r2 <- sum(actual.centered * predicted.centered)^2 / denominator
  # Protect the mathematical bounds against floating-point roundoff.
  min(1, max(0, r2))
}
'''
count = text.count(old)
if count != 1:
    raise RuntimeError(f"Expected exactly one legacy R2 helper, found {count}")
path.write_text(text.replace(old, new, 1))
