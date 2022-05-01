#include <Rcpp.h>
using namespace Rcpp;
#include "partition_map_rcpp.h"
#include "partition_map.h"

//' NNS Partition Map - Internal Use
//'
//' Creates partitions based on partial moment quadrant centroids, iteratively assigning identifications to observations based on those quadrants (unsupervised partitional and hierarchial clustering method).  Basis for correlation, dependence \link{NNS.dep}, regression \link{NNS.reg} routines.
//'
//' @param x a numeric vector.
//' @param y a numeric vector with compatible dimensions to \code{x}.
//' @param type \code{NULL} (default) Controls the partitioning basis.  Set to \code{(type = "XONLY")} for X-axis based partitioning.  Defaults to \code{NULL} for both X and Y-axis partitioning.
//' @param order integer; Number of partial moment quadrants to be generated.  \code{(order = "max")} will institute a perfect fit.
//' @param obs_req integer; (8 default) Required observations per cluster where quadrants will not be further partitioned if observations are not greater than the entered value.  Reduces minimum number of necessary observations in a quadrant to 1 when \code{(obs.req = 1)}.
//' @param min_obs_stop logical; \code{TRUE} (default) Stopping condition where quadrants will not be further partitioned if a single cluster contains less than the entered value of \code{obs.req}.
//' @param noise_reduction the method of determining regression points options for the dependent variable \code{y}: ("mean", "median", "mode", "off"); \code{(noise.reduction = "mean")} uses means for partitions.  \code{(noise.reduction = "median")} uses medians instead of means for partitions, while \code{(noise.reduction = "mode")} uses modes instead of means for partitions.  Defaults to \code{(noise.reduction = "off")} where an overall central tendency measure is used, which is the default for the independent variable \code{x}.
//' @return Returns:
//'  \itemize{
//'   \item{\code{"dt"}} a \link{data.table} of \code{x} and \code{y} observations with their partition assignment \code{"quadrant"} in the 3rd column and their prior partition assignment \code{"prior.quadrant"} in the 4th column.
//'   \item{\code{"regression.points"}} the \link{data.table} of regression points for that given \code{(order = ...)}.
//'   \item{\code{"order"}}  the \code{order} of the final partition given \code{"min.obs.stop"} stopping condition.
//'   }
//'
//' @note \code{min.obs.stop = FALSE} will not generate regression points due to unequal partitioning of quadrants from individual cluster observations.
//'
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' NNS_part_RCPP(x, y)
//'
//' ## Data.table of observations and partitions
//' NNS_part_RCPP(x, y, order = 1)$dt
//'
//' ## Regression points
//' NNS_part_RCPP(x, y, order = 1)$regression.points
//'
//' ## Examine final counts by quadrant
//' DT <- NNS_part_RCPP(x, y)$dt
//' DT[ , counts := .N, by = quadrant]
//' DT
//' @export
// [[Rcpp::export("part_RCPP", rng = false)]]
Rcpp::List NNS_part_RCPP(
  const NumericVector &x,
  const NumericVector &y,
  Rcpp::Nullable<Rcpp::String> type,
  Rcpp::Nullable<RObject> order,
  int obs_req,
  bool min_obs_stop,
  Rcpp::String noise_reduction
) {
  StringVector RP_quadrant(0); 
  StringVector RP_prior_quadrant(0);
  NumericVector RP_x(0);
  NumericVector RP_y(0);
  
  StringVector quadrant(x.size());
  StringVector prior_quadrant(x.size());
  
  bool type_xonly = false;
  if(type.isNotNull()){
	Rcpp::String _t(type);
	type_xonly = (_t=="XONLY");
  }
  int new_order=0;
  bool order_null=false, order_max=false;
  ENUM_NSS_PART_NOISE_REDUCTION nr{};
  if(order.isNull()){
	  order_null=true;
  } else if (is<NumericVector>(order) || is<IntegerVector>(order)) {
	  new_order = as<int>(order);
  } else {
	  order_max = true;
  }
  //noise_reduction = tolower(noise_reduction);   // TODO tolower
  if (noise_reduction=="off")
	  nr=ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_OFF;
  else if (noise_reduction=="mean")
	  nr=ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MEAN;
  else if (noise_reduction=="median")
	  nr=ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MEDIAN;
  else if (noise_reduction=="mode")
	  nr=ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MODE;
  else if (noise_reduction=="mode_class")
	  nr=ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MODE_CLASS;
  else
	  Rcpp::stop("Please ensure noise.reduction is from 'mean', 'median', 'mode' or 'off'");
  
  int ret_order=NNS_part(
	  x,
	  y,
	  quadrant,
	  prior_quadrant,
	  type_xonly,
	  new_order,
	  order_null,
	  order_max,
	  obs_req,
	  min_obs_stop,
	  nr,
	  RP_quadrant, 
      RP_prior_quadrant,
      RP_x, 
      RP_y
  );
  DataFrame dt = DataFrame::create( 
	Named("x") = x,
	Named("y") = y,
	Named("quadrant") = quadrant,
	Named("prior_quadrant") = prior_quadrant
  );
  
  DataFrame rp = DataFrame::create( 
	Named("quadrant") = RP_quadrant,
	Named("x") = RP_x,
	Named("y") = RP_y
//	Named("prior_quadrant") = RP_prior_quadrant
  );
  return Rcpp::List::create(
	Named("order") = ret_order,
	Named("dt") = dt,
	Named("regression.points") = rp
  );
}
