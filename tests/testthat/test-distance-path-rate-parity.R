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
