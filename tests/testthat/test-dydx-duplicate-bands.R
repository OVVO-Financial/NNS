band_means <- function(x, y, p) {
  h_s <- round(exp(seq(log(2), log(floor(sqrt(length(x)))), length.out = 5)))
  h_step <- gravity(diff(sort(unique(x)))) * h_s
  lo <- pmax(min(x), p - h_step)
  hi <- pmin(max(x), p + h_step)
  est <- matrix(unlist(NNS.reg(x, y, plot = FALSE, point.est = c(lo, rep(p, 5), hi),
                               point.only = TRUE, ncores = 1, smooth = TRUE)$Point.est), ncol = 3)
  run_1 <- hi - p; run_2 <- p - lo
  rise_1 <- est[, 3] - est[, 2]; rise_2 <- est[, 2] - est[, 1]
  c(mean((rise_1 + rise_2) / (run_1 + run_2)),
    mean((rise_1 / run_1 - rise_2 / run_2) / ((run_1 + run_2) / 2)))
}

test_that("dy.dx works when small samples repeat a bandwidth", {
  for (n in c(10, 12)) {
    x <- seq(0, 3, length.out = n)
    d <- dy.dx(x, x^2, eval.point = 1.5)
    expect_true(all(is.finite(c(d$first.derivative, d$second.derivative))))
    expect_equal(d$first.derivative, 3, tolerance = 0.05, info = paste("n =", n))
  }
})

test_that("dy.dx averages over all five bandwidth positions", {
  for (n in c(20, 30, 50, 10)) {
    x <- seq(0, 3, length.out = n)
    d <- dy.dx(x, x^2, eval.point = 1.5)
    expect_equal(c(d$first.derivative, d$second.derivative), band_means(x, x^2, 1.5),
                 tolerance = 1e-12, info = paste("n =", n))
  }
})
