test_that("native partition complete-case handling keeps infinities", {
  part <- NNS.part(
    x = c(Inf, Inf, NA_real_, NaN),
    y = c(1, 2, 3, 4),
    order = 1,
    obs.req = 0,
    noise.reduction = "median"
  )

  expect_true(any(is.infinite(part$regression.points$x)))
  expect_false(any(is.na(part$regression.points$x)))
})
