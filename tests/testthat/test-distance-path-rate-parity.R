test_that("bulk distance keeps the session thread settings", {
  rpm <- data.frame(x = 1:3, y.hat = 3:1)
  query <- data.frame(x = 1.5)
  withr::local_envvar(c(RCPP_PARALLEL_NUM_THREADS = "3", RCPP_PARALLEL_STACK_SIZE = "2097152"))
  NNS:::NNS.distance.path.bulk(rpm, query, 2L, ncores = 1L)
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS"), "3")
  expect_identical(Sys.getenv("RCPP_PARALLEL_STACK_SIZE"), "2097152")

  Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
  NNS:::NNS.distance.path.single.bulk(rpm, query, 2L, ncores = 1L)
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_), NA_character_)
})

test_that("corrected exponential ranks agree across distance paths and threads", {
  set.seed(133)
  rpm <- data.frame(a = runif(18), b = runif(18), y.hat = rnorm(18))
  query <- data.frame(a = c(.23, .41, .71), b = c(.37, .63, .19))
  for (cores in c(1L, 2L)) {
    path <- NNS:::NNS.distance.path.bulk(rpm, query, kmax = 8L, ncores = cores)
    for (k in c(1L, 3L, 5L, 8L)) {
      scalar <- vapply(seq_len(nrow(query)), function(i) {
        NNS.distance(rpm, query[i, ], k = k)
      }, numeric(1))
      one <- NNS:::NNS.distance.path.single.bulk(rpm, query, k, ncores = cores)
      bulk <- NNS:::NNS.distance.bulk(rpm, query, k)
      expect_equal(as.numeric(path[, k]), scalar, tolerance = 1e-12)
      expect_equal(as.numeric(one), scalar, tolerance = 1e-12)
      expect_equal(as.numeric(bulk), scalar, tolerance = 1e-12)
    }
  }
})

test_that("a partly matching query cannot change the distance coordinates", {
  rpm <- data.frame(a = c(0, 1, 2), b = c(0, 9, 1), y.hat = c(10, 20, 30))
  pt <- c(b = 9, a = 0)
  typo <- c(b = 9, aa = 0)
  expect_equal(NNS.distance(rpm, pt, 1), 20)
  expect_equal(NNS:::NNS.distance.bulk(rpm, t(pt), 1), 20)
  expect_error(NNS.distance(rpm, c(b = 9, b = 0), 1), "Duplicate")
  expect_error(NNS:::NNS.distance.bulk(rpm, t(c(b = 9, b = 0)), 1), "Duplicate")
  expect_error(NNS.distance(rpm, typo, 1), "aa")
  expect_error(NNS:::NNS.distance.bulk(rpm, t(typo), 1), "aa")
  expect_error(NNS:::NNS.distance.path.bulk(rpm, t(typo), 2), "aa")
})

test_that("distance comparisons require a shared finite coordinate", {
  rpm <- data.frame(a = c(0, 1), b = c(2, NA), y.hat = c(10, 20))
  expect_equal(NNS.distance(rpm, c(1, NA), 1), 20)
  expect_error(NNS.distance(rpm, c(NA_real_, NA_real_), 1), "finite coordinate")
  expect_error(NNS.distance(rpm, c(NA, 2), 1), "finite coordinate")
  bad <- rbind(rpm, c(NA, NA, 999))
  expect_error(NNS.distance(bad, c(.1, .1), 1), "finite coordinate")
  query <- rbind(c(1, NA), c(NA, 2))
  expect_error(NNS:::NNS.distance.bulk(rpm, query, 1), "finite coordinate")
  expect_error(NNS:::NNS.distance.path.bulk(rpm, query, 2), "finite coordinate")
})

test_that("distance ties use the RPM row order before selecting k", {
  rpm <- expand.grid(a = 1:3, b = 1:3)
  rpm$y.hat <- 1:9
  pt <- c(a = 1.5, b = 1.5)
  # The four nearest rows tie. At k = 3 the weights apply to rows 1, 2 and 4.
  expected <- 2.04922830129691
  expect_equal(NNS.distance(rpm, pt, 3), expected)
  expect_equal(NNS:::NNS.distance.bulk(rpm, t(pt), 3), expected)
  expect_equal(NNS:::NNS.distance.path.bulk(rpm, t(pt), 8)[1, 3], expected)
})

test_that("an RPM row tied at the cutoff keeps its original rank", {
  rpm <- data.frame(a = c(1, 1, 1, 0, 0), y.hat = c(1, 0, 0, 1, 1))
  expect_equal(NNS.distance(rpm, 0, 3), 1)
  expect_equal(NNS:::NNS.distance.bulk(rpm, matrix(0), 3), 1)
  expect_equal(NNS:::NNS.distance.path.bulk(rpm, matrix(0), 3)[1, 3], 1)
})
