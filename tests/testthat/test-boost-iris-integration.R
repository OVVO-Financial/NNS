# Full iris integration run mirroring the documented NNS.boost example. Heavy:
# 100 learner trials (exhaustive: 15 subsets), 100 epochs, final NNS.stack.
test_that("NNS.boost iris run has probability threshold, epochs, and a replicated multivariate final stack", {
  skip_on_cran()

  stack_args <- new.env(parent = emptyenv())
  suppressMessages(trace(
    NNS.stack,
    where = asNamespace("NNS"),
    print = FALSE,
    tracer = bquote({
      assign("captured", list(
        method = method,
        folds = folds,
        stack = stack,
        n_predictors = ncol(IVs.train),
        predictor_names = names(IVs.train)
      ), envir = .(stack_args))
    })
  ))
  on.exit(suppressMessages(untrace(NNS.stack, where = asNamespace("NNS"))),
          add = TRUE)

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)

  set.seed(123)
  msgs <- capture_messages(
    a <- NNS.boost(
      iris[1:140, 1:4],
      iris[1:140, 5],
      IVs.test = iris[141:150, 1:4],
      epochs = 100,
      learner.trials = 100,
      type = "CLASS",
      depth = NULL,
      balance = TRUE,
      feature.importance = TRUE,
      status = TRUE
    )
  )

  grDevices::dev.off()

  # 1-3: the console states the probability and the distinct LPM.VaR cutoff,
  # and never prints the obsolete status string.
  expect_true(any(grepl("Learner threshold probability = 0.80", msgs,
                        fixed = TRUE)))
  cutoff_line <- grep("objective cutoff = ", msgs, fixed = TRUE, value = TRUE)
  expect_length(cutoff_line, 1L)
  cutoff <- as.numeric(sub(".*objective cutoff = ([0-9.]+).*", "\\1",
                           cutoff_line))
  expect_true(is.finite(cutoff))
  expect_false(isTRUE(all.equal(cutoff, 0.80)))
  expect_false(any(grepl("Learner Accuracy Threshold", msgs, fixed = TRUE)))

  # 4-5: both diagnostic panels were drawn (nonempty device output).
  expect_true(file.exists(pdf_file))
  expect_gt(file.size(pdf_file), 0)

  # 6: iris has 15 unique feature subsets < 100 trials (exhaustive), yet
  # epochs still ran.
  expect_true(any(grepl("% of epochs", msgs, fixed = TRUE)))

  # 7-9: results and feature outputs. Classification results keep the
  # historical numeric coding (integer codes, base category 1) and
  # $class.levels recovers the label for each code.
  expect_length(a$results, 10L)
  expect_true(is.numeric(a$results))
  expect_true(all(a$results %in% seq_along(levels(iris$Species))))
  expect_identical(a$class.levels, levels(iris$Species))
  expect_gte(mean(a$results == as.numeric(iris[141:150, 5])), 0.8)
  expect_gte(mean(a$class.levels[a$results] ==
                    as.character(iris[141:150, 5])), 0.8)
  expect_equal(sum(a$feature.weights), 1)
  expect_true(all(names(a$feature.frequency) %in% names(iris)[1:4]))
  expect_true(length(a$feature.frequency) >= 1L)

  # 10: the final call used replicated multivariate predictors through
  # NNS.stack(method = 1, folds = 5, stack = FALSE).
  captured <- get("captured", envir = stack_args)
  expect_equal(as.numeric(captured$method), 1)
  expect_equal(as.numeric(captured$folds), 5)
  expect_false(captured$stack)
  expect_gt(captured$n_predictors, 1L)
  base_names <- sub("[.][0-9]+$", "", captured$predictor_names)
  expect_true(all(base_names %in% names(iris)[1:4]))
  expect_false(any(grepl("xstar", captured$predictor_names, fixed = TRUE)))

  # Replication counts follow the historical scaling rule.
  expected_replication <- pmax(1L, as.integer(round(
    as.numeric(a$feature.frequency) / min(as.numeric(a$feature.frequency))
  )))
  expect_identical(captured$n_predictors, sum(expected_replication))
})

test_that("NNS.boost passes a supplied threshold to LPM.VaR as a probability", {
  skip_on_cran()

  set.seed(123)
  msgs <- capture_messages(
    b <- NNS.boost(
      iris[1:140, 1:4],
      iris[1:140, 5],
      IVs.test = iris[141:150, 1:4],
      epochs = 100,
      learner.trials = 100,
      type = "CLASS",
      depth = NULL,
      balance = TRUE,
      threshold = 0.65,
      feature.importance = FALSE,
      status = TRUE
    )
  )

  expect_true(any(grepl("Learner threshold probability = 0.65", msgs,
                        fixed = TRUE)))

  # The supplied probability produces a distinct objective cutoff: 0.65 is
  # not applied as a literal accuracy cutoff.
  cutoff_line <- grep("objective cutoff = ", msgs, fixed = TRUE, value = TRUE)
  expect_length(cutoff_line, 1L)
  cutoff <- as.numeric(sub(".*objective cutoff = ([0-9.]+).*", "\\1",
                           cutoff_line))
  expect_true(is.finite(cutoff))
  expect_false(isTRUE(all.equal(cutoff, 0.65)))

  expect_length(b$results, 10L)
})
