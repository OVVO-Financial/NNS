test_that("NNS.distance weights ranks by dexp(rank, rate = 1/k)", {
  set.seed(3)
  n <- 300
  X <- cbind(a = runif(n), b = runif(n))
  y <- sin(4 * X[, 1]) + X[, 2]^2 + rnorm(n, sd = .1)
  rpm <- NNS.reg(X, y, plot = FALSE, order = 2)$RPM
  pt <- c(a = .37, b = .61)

  mins <- pmin(vapply(rpm[, 1:2], min, 0), pt)
  maxs <- pmax(vapply(rpm[, 1:2], max, 0), pt)
  Z <- sweep(sweep(as.matrix(rpm[, 1:2]), 2, pt), 2, maxs - mins, "/")
  d <- rowSums(abs(Z) + Z^2)
  norm <- function(v) if (all(is.finite(v)) && sum(v) > 0) v / sum(v) else rep(0, length(v))
  reference <- function(k) {
    ids <- order(d)[seq_len(k)]
    s <- d[ids]; r <- seq_len(k)
    w <- norm(rep(1 / k, k) + norm(dt(s, k)) + norm(1 / s) +
      norm(dexp(r, rate = 1 / k)) + norm(abs(rev(dlnorm(r, 0, sd(r), log = TRUE)))) +
      norm(r^-2) + norm(dnorm(s, 0, sd(s))) + norm(exp(-s / (2 * var(s)))))
    sum(w * rpm$y.hat[ids])
  }

  for (k in c(3, 5, 10)) {
    expect_equal(NNS.distance(rpm, pt, k = k), reference(min(k, nrow(rpm))), tolerance = 1e-10)
  }
})
