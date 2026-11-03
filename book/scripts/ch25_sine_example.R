# Run from the repository root with the current NNS version installed:
# Rscript book/scripts/ch25_sine_example.R [output.png]
# An optional PNG path exports the optimizer's current forecast plot.
library(NNS)
RcppParallel::setThreadOptions(numThreads = 1)
options(mc.cores = 1)

set.seed(123)
raw_series <- sin(seq_len(480) / 8) + rnorm(480, sd = .35)
validation_origin <- floor(.8 * length(raw_series))
training <- raw_series[seq_len(validation_origin)]
z <- (raw_series - mean(training)) / sd(training)

# Fit scaling and scan periods at the same earlier validation origin.
seasonal_period <- NNS.seas(z[seq_len(validation_origin)], plot = FALSE)
print(head(seasonal_period$all.periods))

output <- commandArgs(trailingOnly = TRUE)
if (length(output) > 1L) stop("Usage: Rscript book/scripts/ch25_sine_example.R [output.png]")
local({
  if (length(output)) {
    png(output[1L], width = 1600, height = 1000, res = 160)
    on.exit(dev.off(), add = TRUE)
  }
  fit <- NNS.ARMA.optim(z, h = 48, training.set = validation_origin,
                       seasonal.factor = seasonal_period$periods,
                       plot = length(output) == 1L, ncores = 1,
                       print.trace = FALSE)
  print(fit)
})

sessionInfo()
