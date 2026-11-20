# Values
x <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
y <- c(10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0)
z <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)
A <- cbind(x, y, z)

test_that(
  "NNS.norm - linear matrix path equalizes means", {
    N <- NNS.norm(A, linear = TRUE)
    expect_equal(dim(N), dim(A))
    expect_equal(as.numeric(colMeans(N)), rep(mean(colMeans(A)), 3), tolerance = 1e-12)
  }
)

test_that(
  "NNS.norm - nonlinear matrix path returns finite matrix", {
    N <- NNS.norm(A)
    expect_equal(dim(N), dim(A))
    expect_true(all(is.finite(N)))
  }
)

test_that(
  "NNS.norm - equal-length list matches matrix path with linear = FALSE", {
    N_matrix <- NNS.norm(A)
    N_list <- NNS.norm(list(x, y, z))
    expect_equal(unname(N_list), unname(N_matrix), tolerance = 1e-12)
  }
)

test_that(
  "NNS.norm - unequal-length list forces linear scaling", {
    vec1 <- c(1, 2, 3, 4, 5, 6, 7)
    vec2 <- c(10, 20, 30, 40, 50, 60)
    vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
    N <- NNS.norm(list(vec1, vec2, vec3))
    expect_true(is.list(N))
    expect_equal(lengths(N, use.names = FALSE), c(7L, 6L, 9L))
    m <- c(mean(vec1), mean(vec2), mean(vec3))
    expect_equal(unname(sapply(N, mean)), rep(mean(m), 3), tolerance = 1e-12)
  }
)

test_that(
  "NNS.norm - zero-sum length diffs still detected as unequal", {
    vec1 <- c(1, 2, 3, 4, 5)
    vec2 <- c(10, 20, 30, 40, 50, 60)
    vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    N <- NNS.norm(list(vec1, vec2, vec3))
    expect_true(is.list(N))
    expect_equal(lengths(N, use.names = FALSE), c(5L, 6L, 5L))
    m <- c(mean(vec1), mean(vec2), mean(vec3))
    expect_equal(unname(sapply(N, mean)), rep(mean(m), 3), tolerance = 1e-12)
  }
)

test_that("NNS matrix footer gives the full print call", {
  x <- NNS.norm(cbind(1:20, 21:40), linear = TRUE)
  footer <- tail(capture.output(print(x)), 1)
  expect_match(footer, "Use print(x, n = Inf) for all.", fixed = TRUE)
  expect_identical(capture.output(print(x, n = Inf)), capture.output(print(unclass(x))))
})

test_that("linear normalization retains a single row", {
  X <- rbind(first = c(a = 1, b = 2, c = 3))
  N <- NNS.norm(X, linear = TRUE)
  expect_equal(dim(N), c(1L, 3L))
  expect_equal(as.numeric(N), c(2, 2, 2))
  expect_equal(dimnames(N), list("first", c("a Normalized", "b Normalized", "c Normalized")))
})

test_that("normalization accepts a data frame with a numeric list column", {
  X <- data.frame(a = I(as.list(1:4)), b = 4:1)
  M <- cbind(a = 1:4, b = 4:1)
  expect_identical(NNS.norm(X), NNS.norm(M))
  expect_identical(NNS.norm(X, linear = TRUE), NNS.norm(M, linear = TRUE))
})

test_that("normalization checks missing observations inside list elements", {
  expect_error(NNS.norm(list(1:4, c(5, NA, 4))), "missing values")
  expect_error(NNS.norm(list(1:3, c(5, NaN, 4))), "missing values")
})
