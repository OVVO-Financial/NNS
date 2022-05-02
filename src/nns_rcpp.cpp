// example from: https://github.com/r-pkg-examples/rcpp-headers-src
//[[Rcpp::depends(RcppParallel)]]
//[[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>

// Load directory header files
#include "partial_moments_rcpp.h"
#include "partition_map_rcpp.h"

