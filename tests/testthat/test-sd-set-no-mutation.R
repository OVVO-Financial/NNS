test_that("NNS.SD.efficient.set leaves the caller's matrix unnamed", {
  set.seed(1)
  m <- matrix(rnorm(300), 100, 3)
  m2 <- m
  named <- m + 0
  colnames(named) <- c("X_1", "X_2", "X_3")

  out <- NNS.SD.efficient.set(m, 1, status = FALSE)

  expect_null(colnames(m))
  expect_null(colnames(m2))
  expect_identical(out, NNS.SD.efficient.set(named, 1, status = FALSE))
})
