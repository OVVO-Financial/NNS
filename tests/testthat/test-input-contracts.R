test_that("factor thresholds retain their type inside lists", {
  corners <- as.matrix(expand.grid(c(-1, 1), c(-1, 1)))
  expect_error(Co.LPM_nD(corners, list(factor("0"), 0)), "`target` must be numeric, not a factor.", fixed = TRUE)
  expect_error(LPM.VaR(list(factor("0.05")), 0, 1:101), "`percentile` must be numeric, not a factor.", fixed = TRUE)
  expect_identical(LPM.VaR(list(0.05), 0, 1:101), LPM.VaR(.05, 0, 1:101))
})

test_that("dependence compares displayed numeric labels", {
  f <- factor(rep(c(10, 20, 100, 50, 10, 100, 50, 20), 10), levels = c(100, 10, 50, 20))
  x <- as.numeric(as.character(f))
  expect_equal(NNS.dep(f, x), NNS.dep(x, x), tolerance = 1e-12)
  expect_equal(NNS.dep(x, f), NNS.dep(x, x), tolerance = 1e-12)
  expect_equal(NNS.dep(data.frame(x = f), x), NNS.dep(x, x), tolerance = 1e-12)
  expect_equal(NNS.dep(data.frame(a = f, b = x)), NNS.dep(data.frame(a = x, b = x)),
               tolerance = 1e-12)
  nominal <- factor(letters[as.integer(f)])
  expect_equal(NNS.dep(nominal, x), NNS.dep(as.numeric(nominal), x), tolerance = 1e-12)
})

test_that("a one-column bootstrap counts observations and preserves time", {
  x <- stats::ts(c(4, 1, 3, 2, 7, 3), start = c(2000, 1), frequency = 4)
  withr::local_seed(314)
  expected <- NNS.meboot(x, reps = 5, rho = 0)
  set.seed(314)
  df <- data.frame(x = x)
  expect_equal(NNS.meboot(df, reps = 5, rho = 0), expected, tolerance = 1e-12)
  tbl <- structure(df, class = c("tbl_df", "tbl", "data.frame"))
  set.seed(314)
  expect_equal(NNS.meboot(tbl, reps = 5, rho = 0), expected, tolerance = 1e-12)
  dt <- structure(df, class = c("data.table", "data.frame"))
  set.seed(314)
  expect_equal(NNS.meboot(dt, reps = 5, rho = 0), expected, tolerance = 1e-12)
})

test_that("a distance needs an observed coordinate after conversion", {
  rpm <- data.frame(a = c(0, 1), y.hat = c(0, 100))
  expect_error(NNS.distance(rpm, numeric(), k = 1), "Each RPM row and query must share a finite coordinate.", fixed = TRUE)
  expect_error(NNS.distance(rpm, c(a = Inf), k = 1), "Each RPM row and query must share a finite coordinate.", fixed = TRUE)
})

test_that("CDF keeps observation storage and the surviving row labels", {
  d <- data.frame(a = c(1L, NA_integer_, 4L), b = c(4L, 2L, 1L),
                  row.names = c("first", "missing", "last"))
  curve <- NNS.CDF(d, plot = FALSE)$Function
  expect_identical(rownames(curve), c("first", "last"))
  expect_identical(curve$a, c(1L, 4L))
  expect_identical(curve$b, c(4L, 1L))
})

test_that("a distance query needs a measurement shared by its input coordinates", {
  rpm <- data.frame(a = 0:1, y.hat = 0:1)
  expect_error(NNS.distance(rpm, factor("1")), "`dist.estimate` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.distance(rpm, c(a = NA_real_)), "Each RPM row and query must share a finite coordinate.", fixed = TRUE)
  expect_error(NNS.distance(rpm[2:1, ], c(a = NA_real_)), "Each RPM row and query must share a finite coordinate.", fixed = TRUE)
})

test_that("measurement factors are rejected before their codes reach a calculation", {
  f <- factor(1:12)
  x <- seq_along(f)
  expect_error(NNS.ARMA(f, h = 2, plot = FALSE),
               "`variable` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.ARMA.optim(f, training.set = 8),
               "`variable` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(f, x),
               "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(x, f),
               "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SS(ordered(f), x),
               "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.copula(cbind(x, x), target = factor(c(3, 3)), plot = FALSE),
               "`target` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.gravity(f),
               "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.mode(f),
               "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.part(x, f),
               "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.rescale(f, 5, 10),
               "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.rescale(x, 100, .05, "riskneutral", factor("2")),
               "`T` must be numeric, not a factor.", fixed = TRUE)
  expect_error(dy.d_(cbind(x, x), f, 1),
               "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(dy.dx(x, x, factor("3")),
               "eval.point must be numeric point(s) or 'overall'", fixed = TRUE)
  expect_error(NNS.ANOVA(list(a = x, b = f), plot = FALSE),
               "`control` must be numeric, not a factor.", fixed = TRUE)
})

test_that("only overall is a text derivative query", {
  x <- 1:20
  msg <- "eval.point must be numeric point(s) or 'overall'"
  expect_error(dy.dx(x, x^2, "OVERALL"), msg, fixed = TRUE)
  expect_error(dy.dx(x, x^2, as.character(x)), msg, fixed = TRUE)
  expect_error(dy.dx(x, x^2, "1.75"), msg, fixed = TRUE)
  expect_error(dy.dx(x, x^2, "typo"), msg, fixed = TRUE)
})

test_that("series have one measurement per observation", {
  x <- 1:24
  dates <- as.Date("2024-01-01") + x
  tbl <- structure(data.frame(x, dates), class = c("tbl_df", "tbl", "data.frame"))
  msg <- "`variable` must have exactly one column."
  expect_error(NNS.ARMA(tbl, h = 2), msg, fixed = TRUE)
  expect_error(NNS.ARMA(tbl[2:1], h = 2), msg, fixed = TRUE)
  expect_error(NNS.ARMA(cbind(x, x), h = 2), msg, fixed = TRUE)
  expect_error(NNS.seas(cbind(x, x), plot = FALSE), msg, fixed = TRUE)
  expect_error(NNS.meboot(cbind(x, x)), "`x` must have exactly one column.", fixed = TRUE)
  expect_error(NNS.meboot(numeric(), rho = 0), "`x` must contain at least two observations.", fixed = TRUE)
  expect_error(NNS.MC(data.frame(x = 1)), "`x` must contain at least two observations.", fixed = TRUE)
})

test_that("dominance respects column shape and numeric frame labels", {
  x <- c(10, 20, 30)
  df <- data.frame(a = factor(x))
  expect_identical(NNS.FSD.uni(df, 5:6), NNS.FSD.uni(x, 5:6))
  expect_error(NNS.FSD.uni(df$a, 5:6), "x must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.FSD.uni(cbind(df, b = x), 5:6),
               "x must be a single-column object or numeric vector.", fixed = TRUE)
  expect_identical(NNS.FSD.uni(c(-Inf, 0, 1), c(-Inf, -1, 0)),
                   NNS.FSD.uni(c(-100, 0, 1), c(-100, -1, 0)))
})

test_that("minmax propagates infinite measurements", {
  finite <- NNS.rescale(c(2, 4), 5, 10)
  expect_equal(NNS.rescale(c(2, 4, Inf, -Inf), 5, 10), c(finite, Inf, -Inf))
})

test_that("partition predictors keep level order and infinite coordinates", {
  x <- factor(c("b", "a", "b", "c"), levels = c("c", "b", "a"))
  y <- c(-Inf, 1, Inf, NA)
  expect_equal(NNS.part(x, y), NNS.part(as.numeric(x), y), tolerance = 1e-12)
  expect_equal(NNS.part(c(-Inf, 1, Inf, NA), 1:4)$dt$x, c(-Inf, 1, Inf))
})

test_that("mean target aliases use the matrix column means", {
  m <- cbind(a = 1:4, b = 4:1)
  expected <- PM.matrix(1, 1, colMeans(m), m, FALSE)
  expect_equal(PM.matrix(1, 1, NULL, m, FALSE), expected, tolerance = 1e-12)
  expect_equal(PM.matrix(1, 1, "mean", m, FALSE), expected, tolerance = 1e-12)
  expect_equal(PM.matrix(1, 1, "means", m, FALSE), expected, tolerance = 1e-12)
  expect_equal(PM.matrix(1, 1, NULL, m > 2, FALSE), PM.matrix(1, 1, NULL, (m > 2) * 1, FALSE))
})

test_that("only risk-neutral scaling reads maturity", {
  expect_identical(NNS.rescale(1:3, 0, 1, T = factor("2")), NNS.rescale(1:3, 0, 1))
  expect_error(NNS.rescale(95:105, 100, .05, c(m = "RISKNEUTRAL"), T = factor("2")),
               "`T` must be numeric, not a factor.", fixed = TRUE)
})

test_that("unused factor levels do not change dependence", {
  x <- rep(c(10, 20, 100, 50), 12)
  y <- x + rep(c(0, 2, -1), 16)
  f <- factor(x, levels = c(100, 10, 50, 20, "other"))
  expect_identical(NNS.dep(f, y, asym = TRUE), NNS.dep(droplevels(f), y, asym = TRUE))
  expect_identical(NNS.dep(data.frame(x = f, y), asym = TRUE),
                   NNS.dep(data.frame(x = droplevels(f), y), asym = TRUE))
})

test_that("bootstrap evaluates supplied measurements and rho in call order", {
  withr::local_seed(31)
  actual <- NNS.meboot(rnorm(24), reps = 3, rho = runif(1, 0, .2))
  set.seed(31)
  x <- rnorm(24)
  rho <- runif(1, 0, .2)
  expect_identical(actual, NNS.meboot(x, reps = 3, rho = rho))
})

test_that("a bootstrap with no draws accepts a singleton", {
  expect_identical(NNS.meboot(3, reps = 3), list())
  expect_identical(NNS.meboot(3, reps = 3, rho = NULL), list())
  expect_error(NNS.meboot(3, reps = 3, target_drift = 1),
               "`x` must contain at least two observations.", fixed = TRUE)
})

test_that("measurement columns and derivative queries use their numeric conversion", {
  assign("as.double.nns_test_scale", function(x, ...) -unclass(x), envir = .GlobalEnv)
  withr::defer(remove("as.double.nns_test_scale", envir = .GlobalEnv))
  x <- rep(c(10, 20, 100, 50), 12)
  y <- x + rep(c(0, 2, -1), 16)
  z <- structure(-x, class = "nns_test_scale")
  expected <- as.numeric(z)
  tbl <- structure(data.frame(x = I(z)), class = c("tbl_df", "tbl", "data.frame"))
  expect_equal(NNS.dep(data.frame(x = I(z), y), asym = TRUE),
               NNS.dep(data.frame(x = expected, y), asym = TRUE))
  expect_equal(dy.dx(tbl, y, "overall"), dy.dx(expected, y, "overall"))
  expect_equal(NNS.CDF(tbl, plot = FALSE), NNS.CDF(expected, plot = FALSE))
  rpm <- data.frame(a = c(-20, 20), y.hat = c(100, 200))
  query <- structure(-15, class = "nns_test_scale")
  expect_identical(NNS.distance(rpm, data.frame(a = I(query)), k = 1),
                   NNS.distance(rpm, as.numeric(query), k = 1))
  withr::local_seed(31)
  boot <- NNS.meboot(z, reps = 3, rho = 0)
  set.seed(31)
  expect_equal(boot, NNS.meboot(expected, reps = 3, rho = 0))
  expect_equal(NNS.ANOVA(list(z, y), plot = FALSE), NNS.ANOVA(list(expected, y), plot = FALSE))
  t <- seq(0, 2 * pi, by = pi / 100)
  point <- structure(-1.75, class = "nns_test_scale")
  expect_equal(dy.dx(t, sin(t), point), dy.dx(t, sin(t), as.numeric(point)))
})
