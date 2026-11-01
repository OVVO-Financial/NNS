x <- as.numeric(AirPassengers)
t <- seq_along(x)
res <- resid(lm(x ~ t))

meboot_reps <- function(...) {
  set.seed(1)
  NNS.meboot(x, reps = 50, rho = 1, ...)["replicates", ]$replicates
}
mean_slope <- function(R) mean(apply(R, 2, function(r) coef(lm(r ~ t))[2]))
mean_detrended_sd <- function(R) mean(apply(R, 2, function(r) sd(resid(lm(r ~ t)))))

test_that("drift = FALSE replicates are a flat baseline plus residual-scale noise", {
  R <- meboot_reps(drift = FALSE)
  expect_lt(abs(mean_slope(R)), 0.1)
  expect_lt(mean_detrended_sd(R), 1.1 * sd(res))
  expect_lt(mean(apply(R, 2, sd)), 1.1 * sd(res))
})

test_that("target_drift and target_drift_scale set the replicate slope", {
  expect_equal(mean_slope(meboot_reps(target_drift = 1)), 1, tolerance = 0.05)
  expect_equal(mean_slope(meboot_reps(target_drift_scale = 0.5)), 0.5 * coef(lm(x ~ t))[[2]], tolerance = 0.05)
  expect_lt(mean_detrended_sd(meboot_reps(target_drift = 1)), 1.1 * sd(res))
})

test_that("scl.adjustment kappa uses the residual mean and sd", {
  set.seed(1)
  b <- NNS.meboot(x, reps = 50, rho = 1, scl.adjustment = TRUE)
  zz <- c(b["xmin", ]$xmin, b["z", ]$z, b["xmax", ]$xmax)
  uv <- (sum((b["desintxb", ]$desintxb - mean(res))^2) + sum(diff(zz^2) / 12)) / length(x)
  expect_equal(b["kappa", ]$kappa, sd(res) / sqrt(uv) - 1, tolerance = 1e-10)
  expect_equal(mean(apply(b["replicates", ]$replicates, 2, sd)),
               mean(apply(meboot_reps(), 2, sd)), tolerance = 0.01)
})
