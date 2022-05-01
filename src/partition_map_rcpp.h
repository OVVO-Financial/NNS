#ifndef NNS_partition_map_RCPP_H
#define NNS_partition_map_RCPP_H

#include <Rcpp.h>
using namespace Rcpp;

List NNS_part_RCPP(
  const NumericVector &x,
  const NumericVector &y,
  Rcpp::Nullable<String> type=R_NilValue,
  Rcpp::Nullable<RObject> order=R_NilValue,
  int obs_req=8,
  bool min_obs_stop=true,
  String noise_reduction="off"
);
#endif // NNS_partition_map_RCPP_H
