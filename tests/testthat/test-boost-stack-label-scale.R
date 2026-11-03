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
