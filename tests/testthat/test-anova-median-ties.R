test_that("median agreement uses the entire CDF jump at a tied median", {
  x <- c(rep(0, 80), rep(1, 20))
  fit <- NNS.ANOVA(x, x, means.only = TRUE, medians = TRUE,
                   confidence.interval = NULL, plot = FALSE)
  maximum <- ((length(x) * 2 - 2) / (length(x) * 2))^2
  expect_equal(fit$Grand_Statistic, 0)
  expect_equal(fit$Control_CDF, 0.8)
  expect_equal(fit$Treatment_CDF, 0.8)
  expect_equal(fit$Certainty, maximum)

  # A common median can have different jump heights in the two groups.
  y <- c(rep(0, 60), rep(2, 40))
  expect_equal(NNS.ANOVA(x, y, means.only = TRUE, medians = TRUE,
                         confidence.interval = NULL, plot = FALSE)$Certainty,
               maximum)
  # The full distribution profile still detects differing tail structure.
  expect_lt(NNS.ANOVA(x, y, means.only = FALSE, medians = TRUE,
                      confidence.interval = NULL, plot = FALSE)$Certainty,
            maximum)
})

test_that("median ties are handled in multi-group and pairwise paths", {
  x <- c(rep(0, 80), rep(1, 20))
  maximum <- (198 / 200)^2
  expect_equal(unname(NNS.ANOVA(list(x, x, x), means.only = TRUE,
                                medians = TRUE, plot = FALSE)), maximum)
  pairwise <- NNS.ANOVA(cbind(x, x, x), means.only = TRUE,
                        medians = TRUE, pairwise = TRUE, plot = FALSE)
  expect_equal(as.numeric(pairwise[lower.tri(pairwise)]), rep(maximum, 3))
})

test_that("median distance agrees with the original score away from atoms", {
  x <- c(-4, -2, 0, 2, 4)
  y <- c(1, 3, 5, 7, 9)
  target <- median(c(x, y))
  expect_false(target %in% c(x, y))
  distance <- max(abs(mean(x <= target) - 0.5), abs(mean(y <= target) - 0.5))
  expected <- (0.5 - distance)^2 / 0.25 * (8 / 10)^2
  expect_equal(NNS.ANOVA(x, y, means.only = TRUE, medians = TRUE,
                         confidence.interval = NULL, plot = FALSE)$Certainty,
               expected)
  expect_equal(NNS.ANOVA(rep(0, 10), rep(1, 10), means.only = TRUE,
                         medians = TRUE, confidence.interval = NULL,
                         plot = FALSE)$Certainty, 0)
})

test_that("mean agreement retains its degree-one target", {
  x <- c(-2, -1, 1, 2)
  expected <- (6 / 8)^2
  expect_equal(NNS.ANOVA(x, 4 * x, means.only = TRUE,
                         confidence.interval = NULL, plot = FALSE)$Certainty,
               expected)
})

test_that("constant groups use their zero-area expectile balance", {
  x <- rep(0.25, 10)
  maximum <- (18 / 20)^2
  for (central.only in c(TRUE, FALSE)) {
    for (use.medians in c(TRUE, FALSE)) {
      fit <- NNS.ANOVA(x, x, means.only = central.only,
                       medians = use.medians, confidence.interval = NULL,
                       plot = FALSE)
      expect_equal(fit$Certainty, maximum)
      # Agreement does not change the public CDF step convention.
      expect_equal(fit$Control_CDF, 1)
      expect_equal(fit$Treatment_CDF, 1)
      expect_equal(NNS.ANOVA(x, x + 1, means.only = central.only,
                             medians = use.medians,
                             confidence.interval = NULL,
                             plot = FALSE)$Certainty, 0)
    }
  }

  expect_equal(unname(NNS.ANOVA(list(x, x, x), plot = FALSE)), maximum)
  pairwise <- NNS.ANOVA(cbind(x, x, x + 1), pairwise = TRUE, plot = FALSE)
  expect_equal(pairwise[1, 2], maximum)
  expect_equal(pairwise[1, 3], 0)

  # The weighted formula must not move a common decimal constant by one ULP.
  unequal <- NNS.ANOVA(rep(0.1, 2), rep(0.1, 4),
                       confidence.interval = NULL, plot = FALSE)
  expect_identical(unequal$Grand_Statistic, 0.1)
  expect_equal(unequal$Certainty, (4 / 6)^2)
})

test_that("a constant and a variable group can agree in mean but not profile", {
  x <- rep(0, 10)
  y <- rep(c(-1, 1), 5)
  maximum <- (18 / 20)^2
  means <- NNS.ANOVA(x, y, means.only = TRUE,
                     confidence.interval = NULL, plot = FALSE)
  profile <- NNS.ANOVA(x, y, confidence.interval = NULL, plot = FALSE)
  expect_equal(means$Certainty, maximum)
  expect_equal(means$Control_CDF, 1)
  expect_equal(means$Treatment_CDF, 0.5)
  expect_gt(profile$Certainty, 0)
  expect_lt(profile$Certainty, means$Certainty)
})
