# NNS.stack classification weights its two models as NNS 13.0 did: by their
# squared objective scores, or inverse squares for a minimized objective.

test_that("the weights follow NNS 13.0's rule", {
  expect_equal(.nns_stack_objective_weights(.9, .6, "max"),
               c(reg = .81, dim.red = .36) / 1.17)
  expect_equal(.nns_stack_objective_weights(2, 1, "min"), c(reg = .2, dim.red = .8))
  # A zero loss counts as 1e-10, so that model takes essentially all the weight.
  expect_equal(.nns_stack_objective_weights(0, 1, "min"), c(reg = 1, dim.red = 0))
  expect_equal(.nns_stack_objective_weights(0, 0, "max"), c(reg = .5, dim.red = .5))
  expect_equal(.nns_stack_objective_weights(NA, .8, "max"), c(reg = 0, dim.red = 1))
})

test_that("nearly equal objective scores give a nearly equal blend", {
  # Out-of-fold accuracies of 130 and 131 of 138 rows. Choosing the blend
  # weight by its own out-of-fold accuracy gave dim.red all the weight, and the
  # stack then misclassified the sixth test row, which reg classifies correctly.
  set.seed(123)
  test.set <- sort(unlist(lapply(split(seq_len(nrow(iris)), iris$Species),
                                 function(rows) sample(rows, 4))))
  b <- NNS.stack(IVs.train = iris[-test.set, 1:4], DV.train = iris[-test.set, 5],
                 IVs.test = iris[test.set, 1:4], type = "CLASS", balance = TRUE,
                 ncores = 1, folds = 5, status = FALSE)
  expect_equal(c(b$OBJfn.reg, b$OBJfn.dim.red), c(130, 131) / 138)
  expect_equal(b$weights, c(reg = 130^2, dim.red = 131^2) / (130^2 + 131^2))
  expect_equal(unname(b$probabilities),
               unname(b$weights[["reg"]] * b$reg.probabilities +
                        b$weights[["dim.red"]] * b$dim.red.probabilities))
  expect_equal(b$stack, as.numeric(iris[test.set, 5]))
})

test_that("a minimized objective weights the models by inverse squares", {
  set.seed(5)
  test.set <- sample(nrow(iris), 30)
  b <- NNS.stack(IVs.train = iris[-test.set, 1:4], DV.train = iris[-test.set, 5],
                 IVs.test = iris[test.set, 1:4], type = "CLASS",
                 obj.fn = expression(mean(predicted != actual)), objective = "min",
                 ncores = 1, folds = 5, status = FALSE)
  s <- pmax(c(b$OBJfn.reg, b$OBJfn.dim.red), 1e-10)
  expect_equal(unname(b$weights), (1 / s^2) / sum(1 / s^2))
})
