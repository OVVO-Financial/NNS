# Benchmark harness for corrected reference vs native stack paths.
# Run from package root after devtools::load_all().
library(NNS)
run_case <- function(native) {
  options(NNS.native.stack = native, NNS.native.mreg = native, NNS.native.univariate = native)
  set.seed(123); x <- data.frame(a = rnorm(150), b = rnorm(150), c = rnorm(150)); y <- x$a + sin(x$b) - x$c^2
  t <- system.time(fit <- NNS.stack(x, y, IVs.test = x[1:25,], method = 1, folds = 3, ncores = 1, status = FALSE))
  data.frame(version = if (native) "native" else "reference", elapsed = unname(t["elapsed"]), nbest = fit$NNS.reg.n.best, checksum = sum(as.numeric(fit$stack)))
}
res <- rbind(run_case(FALSE), run_case(TRUE))
write.csv(res, "benchmarks/stack_native_results.csv", row.names = FALSE)
print(res)
