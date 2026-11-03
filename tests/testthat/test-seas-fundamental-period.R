# NNS.seas reports the fundamental period (#129): a multiple of the true period
# scores as well in expectation, and the lowest score often fell on a multiple.

test_that("the jackknife variance matches a direct delete-one computation", {
  set.seed(5)
  x <- 10 + rexp(57)
  direct <- function(x, m, use_cv) {
    comp <- split(x, (seq_along(x) - 1L) %% m)
    gmd <- function(z) { k <- length(z); sum(abs(outer(z, z, "-"))) / (k * (k - 1)) }
    v <- vapply(comp, function(z) {
      k <- length(z)
      d <- vapply(seq_len(k), function(i) {
        g <- gmd(z[-i]); if (use_cv) g / abs(mean(z[-i])) else g
      }, numeric(1))
      (k - 1) / k * sum((d - mean(d))^2)
    }, numeric(1))
    sum((lengths(comp) / length(x))^2 * v)
  }
  for (m in c(3L, 7L)) for (cv in c(TRUE, FALSE))
    expect_equal(NNS:::.nns_seas_score_var(x, m, cv), direct(x, m, cv))
})

test_that("a monthly cycle is reported at 12, not at a multiple", {
  t <- 1:540
  p <- vapply(1:6, function(s) {
    set.seed(s); x <- 10 + sin(2 * pi * t / 12) + rnorm(540, sd = .3)
    out <- NNS.seas(x, plot = FALSE)
    c(best = out$best.period, first = out$all.periods$Period[1L], periods = out$periods[1L],
      lowest = with(out$all.periods, Period[which.min(Coefficient.of.Variation)]))
  }, numeric(4))
  expect_true(all(p["best", ] == 12))
  expect_true(all(p["first", ] == 12) && all(p["periods", ] == 12))
  expect_true(all(p["lowest", ] %% 12 == 0 & p["lowest", ] > 12))   # the lowest score sat on a multiple
  expect_equal(NNS.seas(nottem, plot = FALSE)$best.period, 12)
})

test_that("a longer true period is not reduced to a divisor", {
  t <- 1:540
  for (s in 1:3) {
    set.seed(s)
    expect_equal(NNS.seas(10 + sin(2 * pi * t / 24) + rnorm(540, sd = .3), plot = FALSE)$best.period, 24)
    set.seed(s)
    mix <- 10 + sin(2 * pi * t / 12) + .5 * sin(2 * pi * t / 24) + rnorm(540, sd = .3)
    expect_equal(NNS.seas(mix, plot = FALSE)$best.period, 24)
  }
})

test_that("the table keeps every scored period and its score", {
  set.seed(2); x <- 10 + sin(2 * pi * (1:540) / 12) + rnorm(540, sd = .3)
  out <- NNS.seas(x, plot = FALSE)
  raw <- NNS:::NNS_seas_cpp(variable = x - coef(lm(x ~ seq_along(x)))[2] * (seq_along(x) - mean(seq_along(x))),
                            modulo = NULL, mod_only = TRUE, max_period = NULL)$all.periods
  expect_setequal(out$all.periods$Period, raw$Period)
  expect_equal(out$all.periods[order(out$all.periods$Period), "Coefficient.of.Variation"],
               raw[order(raw$Period), "Coefficient.of.Variation"])
  expect_equal(out$all.periods$Period[-1L], raw$Period[raw$Period != out$best.period])
})

test_that("the permutation filter applies to the fundamental period", {
  set.seed(1); x <- 10 + sin(2 * pi * (1:540) / 12) + rnorm(540, sd = .3)
  set.seed(3); out <- NNS.seas(x, plot = FALSE, p.value = TRUE, nperm = 30)
  expect_equal(out$best.period, 12)
  set.seed(3); w <- rnorm(300, 50, 5)
  set.seed(1); expect_equal(NNS.seas(w, plot = FALSE, p.value = TRUE, nperm = 30)$best.period, 1L)
})
