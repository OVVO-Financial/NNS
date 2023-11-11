library(testthat)
library(NNS)
Sys.setenv("OMP_THREAD_LIMIT" = 2)
test_check("NNS")
