reps_of <- function(x, reps = 5, ...) {
  set.seed(1)
  NNS.meboot(x, reps = reps, ...)["replicates", ]$replicates
}
expect_replicate_matrix <- function(actual, expected) {
  expect_s3_class(actual, "NNS.matrix")
  expect_identical(dim(actual), dim(expected))
  expect_equal(as.numeric(actual), as.numeric(expected))
}
slopes <- function(R) apply(R, 2, function(r) coef(lm(r ~ seq_along(r)))[[2]])
max_resid <- function(R) max(abs(apply(R, 2, function(r) resid(lm(r ~ seq_along(r))))))

test_that("an exactly linear series bootstraps like a nearly linear one at any rho", {
  x <- as.numeric(1:10)
  set.seed(5)
  near <- x + 1e-9 * rnorm(10)
  q <- qnorm(1:5 / 6)
  clt_means <- 5.5 + (q - mean(q)) / sd(q) * sd(x) / sqrt(10)
  for (type in c("spearman", "pearson", "nnsdep")) {
    for (rho in c(-1, 0, 0.5, 1)) {
      R <- reps_of(x, rho = rho, type = type)
      expect_lt(max_resid(R), 1e-6)
      expect_true(all(slopes(R) > 1 - 1e-6 & slopes(R) < 1.05 + 1e-6))
      expect_equal(unname(sort(colMeans(R))), clt_means, tolerance = 1e-6)
      expect_equal(unname(sort(colMeans(R))),
                   unname(sort(colMeans(reps_of(near, rho = rho, type = type)))), tolerance = 1e-6)
    }
  }
  expect_true(all(is.finite(reps_of(x, rho = 1, scl.adjustment = TRUE))))
  expect_replicate_matrix(reps_of(x, rho = 1, reps = 1), matrix(x, 10, 1))
  expect_equal(unname(colMeans(reps_of(x, rho = 1, drift = FALSE))), rep(5.5, 5), tolerance = 1e-6)

  r <- reps_of(ts(x, frequency = 4), rho = 0.5)
  expect_equal(attr(r, "tsp"), tsp(ts(x, frequency = 4)))
  expect_lt(max_resid(r), 1e-6)
})

test_that("a constant series bootstraps to constant columns", {
  expect_replicate_matrix(reps_of(rep(3, 10), rho = 0.5), matrix(3, 10, 5))
  expect_replicate_matrix(reps_of(rep(3, 10), rho = 0.5, type = "nnsdep"), matrix(3, 10, 5))
  expect_replicate_matrix(reps_of(rep(3, 10), rho = 0.5, scl.adjustment = TRUE), matrix(3, 10, 5))
  expect_replicate_matrix(reps_of(rep(3, 10), target_drift = 0), matrix(3, 10, 5))
})

test_that("NNS.SS and NNS.MC work on an exactly linear series", {
  set.seed(1)
  x <- rnorm(10, 5)
  set.seed(2)
  r <- NNS.SS(x, 1:10, confidence.interval = TRUE, reps = 49)
  expect_true(all(is.finite(r$boot_vals)))
  expect_true(r$lower <= r$p_star && r$p_star <= r$upper)

  set.seed(1)
  mc <- NNS.MC(1:20, reps = 3, by = 0.5)
  expect_length(mc$replicates, 5)
  expect_true(all(is.finite(unlist(mc$replicates))))
})
