test_that("learner threshold maps probabilities through LPM.VaR", {
  scores <- c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

  default_max <- .nns_boost_threshold(NULL, "max", FALSE, scores)
  expect_identical(default_max$probability, 0.80)
  expect_equal(default_max$cutoff, as.numeric(LPM.VaR(0.80, 1, scores)))

  default_min <- .nns_boost_threshold(NULL, "min", FALSE, scores)
  expect_identical(default_min$probability, 0.20)
  expect_equal(default_min$cutoff, as.numeric(LPM.VaR(0.20, 1, scores)))

  extreme_max <- .nns_boost_threshold(NULL, "max", TRUE, scores)
  expect_identical(extreme_max$probability, 1)
  expect_equal(extreme_max$cutoff, as.numeric(LPM.VaR(1, 1, scores)))

  extreme_min <- .nns_boost_threshold(NULL, "min", TRUE, scores)
  expect_identical(extreme_min$probability, 0)
  expect_equal(extreme_min$cutoff, as.numeric(LPM.VaR(0, 1, scores)))
})

test_that("a supplied threshold is used as the LPM.VaR probability", {
  scores <- c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

  supplied <- .nns_boost_threshold(0.65, "max", FALSE, scores)
  expect_identical(supplied$probability, 0.65)
  expect_equal(supplied$cutoff, as.numeric(LPM.VaR(0.65, 1, scores)))

  # The probability is not treated as a literal objective-score cutoff.
  expect_false(isTRUE(all.equal(supplied$cutoff, 0.65)))

  # extreme overrides a supplied probability.
  overridden <- .nns_boost_threshold(0.65, "max", TRUE, scores)
  expect_identical(overridden$probability, 1)
})

test_that("threshold selection does not use fivenum or empirical quantiles", {
  boost_source <- paste(deparse(body(NNS.boost)), collapse = "\n")
  helper_source <- paste(deparse(body(.nns_boost_threshold)), collapse = "\n")

  for (src in c(boost_source, helper_source)) {
    expect_false(grepl("fivenum", src, fixed = TRUE))
    expect_false(grepl("quantile(", src, fixed = TRUE))
    expect_false(grepl("stats::quantile", src, fixed = TRUE))
  }

  # The cutoff comes from LPM.VaR, and epochs reuse the same cutoff variable.
  expect_true(grepl("LPM.VaR", helper_source, fixed = TRUE))
  expect_true(grepl("learner.threshold", boost_source, fixed = TRUE))
  expect_false(grepl("Learner Accuracy Threshold", boost_source, fixed = TRUE))
})
