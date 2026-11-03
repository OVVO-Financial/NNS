mc_rhos <- function(...) {
  set.seed(3)
  m <- NNS.MC(rnorm(20), reps = 2, drift = FALSE, ...)
  as.numeric(sub("rho = ", "", names(m$replicates)))
}

test_that("NNS.MC exp keeps the rho grid inside the bounds", {
  for (b in list(c(0.2, 0.8), c(-0.8, -0.2), c(-0.5, 0.5))) {
    r <- mc_rhos(lower_rho = b[1], upper_rho = b[2], by = 0.2, exp = 2)
    expect_equal(range(r), b, info = paste(b, collapse = ", "))
    expect_true(all(r >= b[1] & r <= b[2]))
  }
})

test_that("NNS.MC exp moves the grid towards upper_rho", {
  r <- mc_rhos(lower_rho = 0.2, upper_rho = 0.8, by = 0.2, exp = 2)
  expect_true(all(r[2:3] > c(0.6, 0.4)))
  r <- mc_rhos(lower_rho = -0.8, upper_rho = -0.2, by = 0.2, exp = 2)
  expect_true(all(r[2:3] > c(-0.4, -0.6)))
})

test_that("NNS.MC keeps its default full-range grid", {
  r <- seq(-1, 1, 0.01)
  for (e in c(1, 2)) {
    old <- rev(c(-(abs(r[r <= 0])^e), r[r > 0]^(1 / e)))
    set.seed(3)
    m <- NNS.MC(rnorm(20), reps = 1, drift = FALSE, exp = e)
    expect_identical(names(m$replicates), paste0("rho = ", old))
  }
})
