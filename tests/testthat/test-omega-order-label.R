test_that("continuous dominance is identified as area order, not ordinary FSD", {
  x <- c(0, 1, 4, 5)
  y <- c(0, 2, 2.5, 3)
  expect_equal(NNS.FSD(x, y, type = "continuous", plot = FALSE), "X OMEGA Y")
  expect_equal(NNS.FSD(y, x, type = "continuous", plot = FALSE), "Y OMEGA X")
  expect_equal(NNS.FSD(x, y, type = "discrete", plot = FALSE), "NO FSD EXISTS")
  expect_equal(NNS.FSD.uni(x, y, type = "continuous"), 1)
  expect_equal(NNS.FSD.uni(x, y, type = "discrete"), 0)

  # The original CDFs cross and an increasing utility can prefer y.
  expect_gt(mean(x <= 1), mean(y <= 1))
  expect_lt(mean(x <= 3), mean(y <= 3))
  expect_lt(mean(x > 1), mean(y > 1))
  expect_equal(NNS.FSD(x, x, type = "continuous", plot = FALSE),
               "NO OMEGA DOMINANCE EXISTS")
})
