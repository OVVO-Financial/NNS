#ifndef NNS_partition_map_RCPP_H
#define NNS_partition_map_RCPP_H

#include <Rcpp.h>
using namespace Rcpp;
List NNS_part_RCPP(
  const NumericVector &x,
  const NumericVector &y,
  bool Voronoi,
  Rcpp::Nullable<String> type,
  Rcpp::Nullable<RObject> order,
  Rcpp::Nullable<int> obs_req,
  bool min_obs_stop,
  String noise_reduction
);
#endif // NNS_partition_map_RCPP_H
