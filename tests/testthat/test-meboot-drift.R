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

test_that("force.clt spreads changed-drift replicate means on the residual scale", {
  set.seed(7)
  tt <- 1:100
  y <- 1000 + 10 * tt + rnorm(100)
  fit <- lm(y ~ tt)
  clt_reps <- function(...) {
    set.seed(1)
    NNS.meboot(y, reps = 10, rho = 1, ...)["replicates", ]$replicates
  }

  # The spread is that of a mean of n = 100 observations (#101).
  flat <- colMeans(clt_reps(drift = FALSE))
  expect_equal(mean(flat), mean(y), tolerance = 1e-6)
  expect_equal(sd(flat), sd(resid(fit)) / sqrt(100), tolerance = 1e-5)

  sloped <- colMeans(clt_reps(target_drift = 5))
  ref <- 5 * tt + resid(fit)
  expect_equal(mean(sloped), mean(y), tolerance = 1e-6)
  expect_equal(sd(sloped), sd(ref) / sqrt(100), tolerance = 1e-5)
})

test_that("drift = FALSE keeps mean(x) as the replicate level", {
  for (clt in c(TRUE, FALSE)) {
    R <- meboot_reps(drift = FALSE, force.clt = clt)
    expect_equal(mean(R), mean(x), tolerance = 1e-3)
    expect_gt(min(R), 0)
  }
})

test_that("a target_drift equal to the original slope up to rounding keeps the default output", {
  ap <- as.numeric(AirPassengers); tt <- seq_along(ap)
  d_lm <- unname(coef(lm(ap ~ tt))[2])   # differs from fast_lm's slope by ~1e-15
  reps <- function(...) { set.seed(1); NNS.meboot(ap, reps = 20, rho = 1, ...)["replicates", ][[1]] }
  base <- reps()
  expect_equal(reps(target_drift = d_lm), base)
  expect_equal(reps(target_drift_scale = 1), base)
  # A real change of drift still takes the changed-drift path.
  expect_false(isTRUE(all.equal(reps(target_drift = 1.01 * d_lm), base)))
  # A target strictly inside the tolerance reproduces the default exactly,
  # which pins the snap to the original slope.
  d <- NNS:::fast_lm(tt, ap)$coef[2]
  tol <- sqrt(.Machine$double.eps) * max(abs(d), sd(ap) / length(ap))
  expect_identical(reps(target_drift = d + tol / 2), base)
})

test_that("the slope tolerance allows for rounding at a large series level", {
  t <- 1:20; set.seed(7); x <- 1e6 + .0001 * t + .01 * rnorm(20)
  d_lm <- unname(coef(lm(x ~ t))[2])   # same OLS fit, ~1e-11 from fast_lm's slope
  reps <- function(...) { set.seed(1); NNS.meboot(x, reps = 5, rho = .5, ...)["replicates", ][[1]] }
  expect_identical(reps(target_drift = d_lm), reps())
  # A real change of drift still takes the changed-drift path.
  expect_false(isTRUE(all.equal(reps(target_drift = .1), reps())))
})

test_that("force.clt spreads replicate means like means of n observations (#101)", {
  set.seed(1); x <- cumsum(rnorm(100))
  spread <- function(r) { set.seed(2); sd(colMeans(NNS.meboot(x, reps = r, rho = 1)["replicates", ]$replicates)) }
  s20 <- spread(20); s80 <- spread(80)
  # Before: 0.841 at 20 replicates and 0.420 at 80 (sd(x) / sqrt(reps)).
  expect_lt(abs(s80 / s20 - 1), 0.2)
  expect_equal(s80, sd(x) / sqrt(length(x)), tolerance = 0.2)
})
