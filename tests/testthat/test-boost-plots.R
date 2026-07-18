.boost_plot_fixture <- function() {
  set.seed(42)
  list(
    learner.results = runif(40, 0.5, 1),
    threshold.probability = 0.80,
    learner.threshold = 0.9,
    feature.frequency = c(a = 5L, b = 3L, c = 1L)
  )
}

test_that("diagnostics draws both panels in one two-panel layout", {
  fx <- .boost_plot_fixture()
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .nns_boost_plot_learner_distribution = function(learner.results,
                                                    threshold.probability,
                                                    learner.threshold) {
      seen$dist_layout <- graphics::par("mfrow")
      seen$dist_threshold <- learner.threshold
      invisible(NULL)
    },
    .nns_boost_plot_feature_frequency = function(feature.frequency) {
      seen$freq_layout <- graphics::par("mfrow")
      invisible(NULL)
    },
    .package = "NNS"
  )

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  .nns_boost_plot_diagnostics(
    learner.results = fx$learner.results,
    threshold.probability = fx$threshold.probability,
    learner.threshold = fx$learner.threshold,
    feature.frequency = fx$feature.frequency
  )

  # Both helpers were invoked, under the same undisturbed two-panel layout.
  expect_identical(seen$dist_layout, c(2L, 1L))
  expect_identical(seen$freq_layout, c(2L, 1L))

  # The threshold line receives learner.threshold, not the probability.
  expect_identical(seen$dist_threshold, fx$learner.threshold)

  # The original graphics parameters are restored after both panels.
  expect_identical(graphics::par("mfrow"), c(1L, 1L))
})

test_that("the learner-distribution panel draws the cutoff line and label", {
  src <- paste(deparse(body(.nns_boost_plot_learner_distribution)),
               collapse = "\n")
  expect_true(grepl("abline", src, fixed = TRUE))
  expect_true(grepl("v = learner.threshold", src, fixed = TRUE))
  expect_true(grepl("LPM.VaR(p = %.2f)", src, fixed = TRUE))
  expect_true(grepl("hist", src, fixed = TRUE))

  freq_src <- paste(deparse(body(.nns_boost_plot_feature_frequency)),
                    collapse = "\n")
  expect_true(grepl("barplot", freq_src, fixed = TRUE))
  expect_true(grepl("horiz = TRUE", freq_src, fixed = TRUE))
})

test_that("NNS.boost invokes both plot helpers only when feature.importance = TRUE", {
  calls <- new.env(parent = emptyenv())
  calls$dist <- 0L
  calls$freq <- 0L

  local_mocked_bindings(
    .nns_boost_plot_learner_distribution = function(...) {
      calls$dist <- calls$dist + 1L
      invisible(NULL)
    },
    .nns_boost_plot_feature_frequency = function(...) {
      calls$freq <- calls$freq + 1L
      invisible(NULL)
    },
    .package = "NNS"
  )

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  run_boost <- function(feature.importance) {
    NNS.boost(
      iris[1:60, 1:4],
      iris[1:60, 5],
      IVs.test = iris[61:70, 1:4],
      epochs = 3,
      learner.trials = 5,
      type = "CLASS",
      feature.importance = feature.importance,
      status = FALSE,
      seed = 123L
    )
  }

  invisible(run_boost(feature.importance = TRUE))
  expect_identical(calls$dist, 1L)
  expect_identical(calls$freq, 1L)

  invisible(run_boost(feature.importance = FALSE))
  expect_identical(calls$dist, 1L)
  expect_identical(calls$freq, 1L)
})

test_that("NNS.boost completes its diagnostics on a PDF device without errors", {
  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)

  expect_no_error(
    invisible(NNS.boost(
      iris[1:60, 1:4],
      iris[1:60, 5],
      IVs.test = iris[61:70, 1:4],
      epochs = 3,
      learner.trials = 5,
      type = "CLASS",
      feature.importance = TRUE,
      status = FALSE,
      seed = 123L
    ))
  )

  grDevices::dev.off()
  expect_true(file.exists(pdf_file))
  expect_gt(file.size(pdf_file), 0)
})
