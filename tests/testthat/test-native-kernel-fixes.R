direct_lpm <- function(d, t, x) vapply(t, function(tt) LPM(d, tt, x), numeric(1))
direct_upm <- function(d, t, x) vapply(t, function(tt) UPM(d, tt, x), numeric(1))

test_that("vectorised LPM/UPM (> 32 targets) match the direct scan at every degree", {
  set.seed(7)
  price <- 400 * cumprod(1 + rnorm(1500, 4e-4, 0.012))
  ret <- rt(1000, df = 3) * 0.01
  for (x in list(price, ret)) {
    t <- c(seq(min(x) - 0.1 * diff(range(x)), max(x) + 0.1 * diff(range(x)), length.out = 80),
           sort(x)[seq(1, length(x), by = 25)])
    for (d in c(0, 1, 2, 3, 4, 8, 12)) {
      l <- LPM(d, t, x); u <- UPM(d, t, x)
      expect_true(all(l >= 0) && all(u >= 0))
      expect_equal(l, direct_lpm(d, t, x), tolerance = 1e-10)
      expect_equal(u, direct_upm(d, t, x), tolerance = 1e-10)
      r <- LPM.ratio(d, t, x)
      expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
    }
  }
  # Degree 2 just below the maximum used to come out negative.
  x <- rnorm(1e4, 100, 1)
  t <- sort(x)[(1e4 - 99):1e4] - 1e-9
  expect_true(all(UPM(2, t, x) >= 0))
})

test_that("negative and NaN degrees are rejected", {
  expect_error(LPM(-1, 0, c(-2, 1)), "non-negative")
  expect_error(UPM(NaN, 0, c(-2, 1)), "non-negative")
  expect_error(LPM.ratio(-0.5, 0, c(-2, 1)), "non-negative")
  expect_error(Co.LPM(-1, c(-2, 1), c(-2, 1), 0, 0), "non-negative")
  expect_error(D.UPM(1, -2, c(-2, 1), c(-2, 1), 0, 0), "non-negative")
  expect_error(PM.matrix(-1, 1, NULL, matrix(rnorm(20), 10), TRUE), "non-negative")
  expect_error(Co.LPM_nD(matrix(rnorm(20), 10), c(0, 0), -1), "non-negative")
  expect_equal(LPM(0.5, 0, c(-4, 1)), 1)
})

test_that("discrete NNS.mode handles values beyond the int range", {
  expect_equal(NNS.mode(c(3e9, 3e9, 3e9, 1, 2), discrete = TRUE), 3e9)
  expect_equal(NNS.mode(c(-5e9, -5e9, -5e9, 1, 2), discrete = TRUE), -5e9)
  expect_equal(NNS.mode(c(1, 2, 2, 3, 3, 4, 4, 5), discrete = TRUE, multi = TRUE), c(2, 3, 4))
})

test_that("NNS.gravity / NNS.mode stay fast with a tight IQR and an outlier", {
  set.seed(1)
  x <- c(5 + 1e-7 * rnorm(1000), 6)
  elapsed <- system.time({ g <- NNS.gravity(x); m <- NNS.mode(x, multi = FALSE) })[["elapsed"]]
  expect_lt(elapsed, 2)
  expect_equal(g, 5, tolerance = 1e-3)
  expect_equal(m, 5, tolerance = 1e-5)
  # Tiny histograms (<= 3 bins) are well defined; tied modes resolve to the
  # one nearest the gravity (0.375), not their average.
  expect_equal(NNS.mode(c(0, 0, 1, 1), multi = FALSE), 0)
  expect_equal(NNS.gravity(c(0, 0, 1, 1)), 0.375)
})

test_that("vector generators do not read out of bounds or divide by zero", {
  g <- NNS:::generate.vectors(c(1, 2, 3), 7L)
  expect_length(g$Component.series[[1]], 0)
  expect_equal(NNS:::generate.vectors(as.numeric(1:10), 3L)$Component.series[[1]], c(2, 5, 8))
  expect_equal(NNS:::generate.vectors(as.numeric(1:10), 3L, 7L)$Component.series[[1]], c(2, 5))
  expect_error(NNS:::generate.vectors(c(1, 2, 3), 2L, 5L), "len")
  lv <- NNS:::generate.lin.vectors(c(1, 2, 3), 7L, 2L)
  expect_length(lv$Component.series[[1]], 0)
})

test_that("vector-target LPM / UPM reject negative, NA and NaN degrees", {
  x <- c(-2, 1)
  for (t in list(c(0, 1), rep(0, 33))) {
    expect_error(LPM(-1, t, x), "non-negative")
    expect_error(UPM(-1, t, x), "non-negative")
    expect_error(LPM(NaN, t, x), "non-negative")
    expect_error(UPM(NA_real_, t, x), "non-negative")
  }
})

test_that("vectorised UPM / UPM.ratio keep relative accuracy just below the maximum", {
  set.seed(11)
  x <- rnorm(1e4, 100, 1)
  t <- sort(x)[(1e4 - 39):1e4] - 1e-6
  expect_equal(UPM(2, t, x) / direct_upm(2, t, x), rep(1, 40), tolerance = 1e-10)
  direct_ratio <- vapply(t, function(tt) UPM.ratio(2, tt, x), numeric(1))
  expect_equal(UPM.ratio(2, t, x) / direct_ratio, rep(1, 40), tolerance = 1e-10)
})

test_that("discrete NNS.mode keeps tied modes beyond the int range apart", {
  x <- c(3e9, 3e9, -5e9, -5e9, 1)
  expect_identical(NNS.mode(x, discrete = TRUE, multi = TRUE), c(-5e9, 3e9))
  expect_identical(NNS.mode(x, discrete = TRUE, multi = FALSE), -1e9)
})

test_that("generate.vectors gives empty components for empty input, len = 0 and lag > n", {
  for (g in list(NNS:::generate.vectors(numeric(0), 2L),
                 NNS:::generate.vectors(c(1, 2, 3, 4), 2L, 0L),
                 NNS:::generate.vectors(c(1, 2, 3), 7L))) {
    expect_identical(g$Component.series[[1]], numeric(0))
    expect_identical(g$Component.index[[1]], integer(0))
  }
})

test_that("generate.lin.vectors(x, 0) errors instead of killing the R session", {
  # Run in a child process: on the unfixed kernel this call raises SIGFPE.
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script))
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               "tryCatch(NNS:::generate.lin.vectors(1:4, 0L),",
               "         error = function(e) cat(conditionMessage(e)))"), script)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                                  stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"))
  expect_match(paste(out, collapse = "\n"), "'l' must be a positive integer")
})

test_that("equally high mode peaks resolve to the one nearest the gravity", {
  x <- rep(0:19, each = 400)
  m <- NNS.mode(x)
  r <- NNS.mode(x, multi = FALSE)
  expect_true(r %in% m)
  expect_equal(r, m[which.min(abs(m - NNS.gravity(x)))])
  # Tied clusters: always a real peak, never a value in the gap between them
  # (averaging them gave 8.91 for this sample).
  y <- c(rep(1:3, each = 300), rep(10:12, each = 300), rep(20:22, each = 150))
  expect_true(NNS.mode(y, multi = FALSE) %in% NNS.mode(y))
})
