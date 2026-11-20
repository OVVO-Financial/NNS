test_that("boost does not decode numeric binary stack results twice", {
  x <- iris[1:100, 1:2]
  query <- x[c(1, 100), ]
  for (values in list(c(-5, 8), c(-.1, .2))) {
    y <- rep(values, each = 50)
    seen <- numeric()
    record <- function(actual) seen <<- c(seen, actual)
    loss <- as.expression(bquote({ .(record)(actual); sum((predicted - actual)^2) }))
    fit <- NNS.boost(x, y, query, obj.fn = loss, learner.trials = 3,
      epochs = 0, threshold = 1, seed = 7, folds = 2, status = FALSE,
      feature.importance = FALSE)
    expect_setequal(unique(seen), values)
    expect_equal(fit$results, values)
  }
})

test_that("boost keeps factor label decisions when class levels are permuted", {
  run <- function(y) NNS.boost(iris[, 1:3], y, iris[c(1, 55, 110), 1:3],
    learner.trials = 3, epochs = 0, balance = TRUE, threshold = 1, type = "CLASS",
    seed = 13, folds = 2, status = FALSE, feature.importance = FALSE)
  a <- run(iris$Species)
  b <- run(factor(iris$Species, levels = rev(levels(iris$Species))))
  expect_equal(a$class.levels[a$results], b$class.levels[b$results])
  expect_equal(a$feature.weights, b$feature.weights)
})

test_that("boost returns double nominal codes and preserves logical labels", {
  x <- iris[1:100, 1:2]
  fit <- function(y) NNS.boost(x, y, x[c(1, 100), ], epochs = 0, folds = 2,
    learner.trials = 3, status = FALSE, feature.importance = FALSE, seed = 9)
  y <- factor(iris$Species[1:100], levels = c("versicolor", "setosa"))
  nominal <- fit(y)
  expect_type(nominal$results, "double")
  expect_identical(nominal$class.levels, levels(y))
  expect_identical(fit(y == levels(y)[2])$results,
                   nominal$class.levels[nominal$results] == levels(y)[2])
})

test_that("explicit boost classification scores every fit on the original scale", {
  x <- iris[1:100, 1:2]
  y <- rep(c(-5, 8), each = 50)
  seen <- numeric()
  record <- function(predicted, actual) seen <<- c(seen, predicted, actual)
  loss <- as.expression(bquote({ .(record)(predicted, actual); mean(predicted == actual) }))
  fit <- NNS.boost(x, y, x[c(1, 100), ], type = "CLASS", obj.fn = loss,
    objective = "max", threshold = 0, learner.trials = 3, epochs = 2, folds = 2, seed = 9,
    feature.importance = FALSE, status = FALSE)
  expect_setequal(seen, y)
  expect_equal(fit$results, y[c(1, 100)])
})

test_that("boost keeps ordered response levels in a one-column frame", {
  x <- iris[1:100, 1:2]
  y <- ordered(iris$Species[1:100], levels = rev(levels(iris$Species)))
  fit <- function(y) NNS.boost(x, y, x[c(1, 100), ], epochs = 0, folds = 2,
    learner.trials = 3, status = FALSE, feature.importance = FALSE, seed = 9)
  control <- fit(y)
  expect_identical(fit(data.frame(y)), control)
})
