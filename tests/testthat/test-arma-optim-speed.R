# NNS.ARMA.optim run time: the same models with less repeated work.

test_that("NNS.ARMA.optim computes each validation forecast once", {
  orig <- NNS::NNS.ARMA
  keys <- character()
  testthat::local_mocked_bindings(
    NNS.ARMA = function(variable, ...) {
      a <- list(...)
      keys <<- c(keys, paste(length(variable), deparse(a$training.set), a$h, a$method,
                             paste(a$seasonal.factor, collapse = ","),
                             paste(a$weights, collapse = ","), isTRUE(a$shrink),
                             sep = "|"))
      orig(variable, ...)
    },
    .package = "NNS"
  )
  set.seed(9)
  t <- 1:160
  x <- 50 + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
  fit <- NNS.ARMA.optim(x, h = 10, seasonal.factor = c(7, 12), ncores = 1,
                        print.trace = FALSE,
                        obj.fn = expression(mean((predicted - actual)^2)))
  # The selected multi-period model and its equal-weight candidate were both
  # requested twice before.
  expect_equal(as.numeric(fit$periods), c(7, 12))
  expect_equal(fit$weights, c(0.5, 0.5))
  expect_false(anyDuplicated(keys) > 0)
})

test_that("the nonlinear step uses the squared-index weighted rise and run", {
  x <- as.numeric(AirPassengers)[1:120]
  for (lag in c(1, 12)) {
    gv <- NNS:::generate.vectors(x, lag)
    ci <- gv$Component.index[[1]]
    cs <- gv$Component.series[[1]]
    rp <- NNS.reg(ci, cs, type = "REG", return.values = FALSE, plot = FALSE,
                  multivariate.call = TRUE)
    xs <- rp$x[nrow(rp)] - rp$x
    ys <- rp$y[nrow(rp)] - rp$y
    xs <- xs[-length(xs)]
    ys <- ys[-length(ys)]
    expected <- tail(cs, 1) +
      mean(rep(ys, seq_along(ys)^2)) / mean(rep(xs, seq_along(xs)^2))
    expect_equal(NNS.ARMA(x, h = 1, seasonal.factor = lag, method = "nonlin",
                          plot = FALSE), expected, tolerance = 1e-12)
  }
})

test_that("partition data frames are the ones as.data.frame() builds", {
  set.seed(3)
  x <- rnorm(60)
  y <- x + rnorm(60)
  for (type in list(NULL, "XONLY")) {
    out <- NNS:::NNS_part_cpp(x, y, type = type, order_in = 3L, obs_req = 8L,
                              min_obs_stop = TRUE, noise_reduction = "off")
    for (nm in c("dt", "regression.points", "segments_h", "segments_v")) {
      expect_identical(out[[nm]],
                       as.data.frame(as.list(out[[nm]]), stringsAsFactors = FALSE))
    }
  }
})

test_that("small partial-moment inputs agree with the threaded path", {
  # 30 x 3 runs on the calling thread; 200 stacked copies are dispatched.
  set.seed(4)
  m <- matrix(rnorm(90), 30, 3)
  big <- m[rep(seq_len(30), 200), ]
  tg <- colMeans(m)
  for (deg in c(0, 1)) {
    expect_equal(DPM_nD(m, tg, deg, TRUE), DPM_nD(big, tg, deg, TRUE), tolerance = 1e-12)
    expect_equal(Co.LPM_nD(m, tg, deg, FALSE), Co.LPM_nD(big, tg, deg, FALSE),
                 tolerance = 1e-12)
    expect_equal(PM.matrix(deg, deg, tg, m, FALSE), PM.matrix(deg, deg, tg, big, FALSE),
                 tolerance = 1e-12)
  }
})

test_that("a plain numeric predictor is encoded as its one-column data frame", {
  enc <- NNS:::.nns_reg_encode_predictors
  x <- c(a = 3, b = 1, c = 2)
  expect_identical(enc(x), enc(data.frame(x = x)))
  expect_identical(enc(1:5), enc(data.frame(x = 1:5)))
  expect_error(enc(c(1, NA)), "finite numeric")
})

test_that("regression points match the data frame construction", {
  # The data frame version the vector code replaced.
  reference <- function(x, y, order, noise.reduction = "off") {
    reducer <- function(z) NNS:::.nns_reg_reduce_value(z, noise.reduction, FALSE)
    if (identical(order, "max")) {
      s <- split(y, x)
      rp <- data.frame(x = as.numeric(names(s)), y = vapply(s, reducer, numeric(1L)))
    } else {
      part <- NNS.part(x, y, type = "XONLY", order = order, obs.req = 0,
                       min.obs.stop = TRUE, noise.reduction = noise.reduction)
      rp <- as.data.frame(part$regression.points[, c("x", "y"), drop = FALSE])
    }
    rp <- rp[is.finite(rp$x) & is.finite(rp$y), , drop = FALSE]
    rp <- rp[order(rp$x, method = "radix"), , drop = FALSE]
    merge <- function(rp) {
      if (!anyDuplicated(rp$x)) return(rp)
      g <- split(rp$y, rp$x)
      rp <- data.frame(x = as.numeric(names(g)), y = vapply(g, reducer, numeric(1L)))
      rp[order(rp$x, method = "radix"), , drop = FALSE]
    }
    rp <- merge(rp)
    rp <- rbind(rp, data.frame(x = c(min(x), max(x)),
                               y = c(reducer(y[x == min(x)]), reducer(y[x == max(x)]))))
    rp <- merge(rp[order(rp$x, method = "radix"), , drop = FALSE])
    rownames(rp) <- NULL
    rp
  }
  set.seed(5)
  cases <- list(list(rnorm(40), rnorm(40)),
                list(round(runif(50) * 4), rnorm(50)),
                list(1:30, cumsum(rnorm(30))))
  for (cs in cases) for (ord in list(2L, 4L, "max")) for (nr in c("off", "median")) {
    got <- NNS:::.nns_reg_build_points(cs[[1]], cs[[2]], ord, nr, FALSE)
    expect_identical(got[c("x", "y")], reference(cs[[1]], cs[[2]], ord, nr))
    expect_false(is.null(attr(got, "nns.order")))
  }
})

test_that("linear.approximation defaults to TRUE", {
  # FALSE scores every candidate set with the nonlinear methods; the default
  # passes the linear search's selection to the nonlinear comparison.
  expect_identical(eval(formals(NNS.ARMA.optim)$linear.approximation), TRUE)
})
