a <- c(-1.2, 0.3, 0.8, 2.5, 1.1, -0.4)
b <- c(0.1, 1.4, 2.2, 3.0, 3.9, 4.6, 0.5, 2.7, 1.9, -0.3)

test_that("list and two-sample paths share the size-weighted grand mean", {
  two <- NNS.ANOVA(a, b, plot = FALSE, confidence.interval = NULL)
  expect_equal(two$Grand_Statistic, mean(c(a, b)))
  expect_equal(two$Certainty, 0.02320132, tolerance = 1e-6)
  expect_equal(NNS.ANOVA(list(a, b), plot = FALSE), c(Certainty = two$Certainty))
})

x <- 1:9
y <- c(4, 7, 15, 20)

test_that("grand median is the median of the pooled observations", {
  r <- NNS.ANOVA(x, y, medians = TRUE, means.only = TRUE, plot = FALSE,
                 confidence.interval = NULL)
  expect_identical(r$Grand_Statistic, 6)
})

test_that("multi-group medians evaluate at the pooled median", {
  expected <- NNS:::NNS.ANOVA.bin(x, y, medians = TRUE, means.only = TRUE,
                                  mean.of.means = 6, plot = FALSE)$Certainty
  expect_equal(expected, 0.1789941, tolerance = 1e-6)
  expect_equal(NNS.ANOVA(list(x, y), medians = TRUE, means.only = TRUE, plot = FALSE),
               c(Certainty = expected))
})

test_that("a constant group and a variable group with equal means use that mean", {
  for (k in c(0.1, 0.5, 2.7, 100.1)) {
    x <- rep(k, 2)
    y <- k + c(-1/16, 1/16, 0, 0)
    r <- NNS.ANOVA(x, y, means.only = TRUE, confidence.interval = NULL, plot = FALSE)
    expect_identical(r$Grand_Statistic, mean(y))
    expect_equal(r$Certainty, 4/9)
  }
  x <- rep(.Machine$double.xmax, 3)
  expect_equal(NNS.ANOVA(x, x, means.only = TRUE, confidence.interval = NULL, plot = FALSE)$Certainty, 4/9)
})
