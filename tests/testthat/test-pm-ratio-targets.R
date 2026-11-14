x <- c(1, 2, 3, 10)

test_that("LPM.ratio and UPM.ratio use integer targets as given", {
  expect_equal(LPM.ratio(0, 2L, x), 0.5)
  expect_equal(LPM.ratio(0, 1:3, x), c(0.25, 0.5, 0.75))
  expect_equal(UPM.ratio(0, 1:3, x), c(0.75, 0.5, 0.25))
  expect_equal(LPM.ratio(1, 1:3, x), LPM.ratio(1, c(1, 2, 3), x))
})

test_that("Co.LPM, Co.UPM, D.LPM and D.UPM use integer targets as given", {
  expect_equal(Co.LPM(1, x, x, 2L, 2L), 0.25)
  expect_equal(Co.UPM(1, x, x, 2L, 2L), Co.UPM(1, x, x, 2, 2))
  expect_equal(D.LPM(1, 1, x, rev(x), 0L, 0L), 0)
  expect_equal(D.UPM(1, 1, x, rev(x), 3L, 2L), D.UPM(1, 1, x, rev(x), 3, 2))
})

test_that("PM.matrix uses an integer target as given", {
  A <- cbind(x, rev(x))
  expect_equal(PM.matrix(1, 1, c(0L, 0L), A, FALSE)$cupm[1, 1], mean(x^2))
  expect_equal(PM.matrix(1, 1, c(2L, 3L), A, FALSE), PM.matrix(1, 1, c(2, 3), A, FALSE))
})

test_that("NNS.CDF of an integer vector is the ECDF", {
  expect_equal(NNS.CDF(1:10, plot = FALSE)$Function$CDF, (1:10) / 10)
})

test_that("ratio and co-moment wrappers reject non-finite input", {
  msg <- "must be finite"
  expect_error(LPM.ratio(0, 2, c(1, 2, NA)), msg)
  expect_error(UPM.ratio(0, Inf, 1:3), msg)
  expect_error(Co.LPM(0, c(1, 2, NaN), 1:3, 2, 2), msg)
  expect_error(Co.UPM(0, 1:3, 1:3, 2, -Inf), msg)
  expect_error(D.LPM(0, 0, 1:3, c(1, NA, 3), 2, 2), msg)
  expect_error(D.UPM(0, 0, 1:3, 1:3, NA_real_, 2), msg)
  expect_error(PM.matrix(1, 1, c(2, 2), cbind(1:3, c(1, 2, NA)), FALSE), msg)
})

test_that("NNS.CDF drops missing observations", {
  y <- c(1, 2, NA, 4)
  expect_equal(NNS.CDF(y, plot = FALSE)$Function$CDF, c(1, 2, 3) / 3)
  expect_equal(NNS.CDF(y, target = 2, plot = FALSE)$target.value, 2 / 3)
})

test_that("LPM counts ties at degree 0 with excess_ret = TRUE", {
  expect_equal(LPM(0, 2, 1:3, excess_ret = TRUE), 2 / 3)
  expect_equal(LPM(0, 2, 1:3, TRUE) + UPM(0, 2, 1:3, TRUE), 1)
})

test_that("LPM.ratio and UPM.ratio fill only the 0/0 case at degree > 0", {
  expect_equal(LPM.ratio(1, 3, rep(3, 4)), 1)
  expect_equal(UPM.ratio(2, 3, rep(3, 4)), 0)
  expect_equal(LPM.ratio(1, rep(3, 40), 3), rep(1, 40))
  expect_equal(NNS.CDF(rep(2, 4), 1, plot = FALSE)$Function$CDF, rep(1, 4))
  na_or_error <- function(expr) tryCatch(all(is.na(expr)), error = function(e) TRUE)
  expect_true(na_or_error(LPM.ratio(NA_real_, 2, 1:3)))
  expect_true(na_or_error(UPM.ratio(NA_real_, 2, 1:3)))
})

test_that("degree > 0 CDF of a large constant sample is identically 1", {
  # Exercise the large-input path without a machine-specific elapsed-time cap.
  z <- NNS.CDF(rep(2, 1e5), 1, plot = FALSE)
  expect_true(all(z$Function$CDF == 1))
})
