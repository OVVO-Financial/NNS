#ifndef NNS_partition_map_H
#define NNS_partition_map_H

#include <Rcpp.h>
using namespace Rcpp;
enum class ENUM_NSS_PART_NOISE_REDUCTION {
  NOISE_REDUCTION_OFF, 
  NOISE_REDUCTION_MEAN, 
  NOISE_REDUCTION_MEDIAN, 
  NOISE_REDUCTION_MODE,
  NOISE_REDUCTION_MODE_CLASS
};

size_t NNS_part(
  const NumericVector &x,
  const NumericVector &y,
  CharacterVector &quadrant,
  CharacterVector &prior_quadrant,
  bool &type_xonly,
  int order,
  bool order_null,
  bool order_max,
  size_t obs_req,
  bool min_obs_stop,
  ENUM_NSS_PART_NOISE_REDUCTION noise_reduction,
  CharacterVector &RP_quadrant, 
  CharacterVector &RP_prior_quadrant,
  NumericVector &RP_x, 
  NumericVector &RP_y
);
#endif // NNS_partition_map_H
