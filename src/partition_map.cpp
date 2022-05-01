#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include "partition_map.h"
using namespace Rcpp;

size_t NNS_part(
  const NumericVector &x,
  const NumericVector &y,
  StringVector &quadrant,
  StringVector &prior_quadrant,
  bool &type_xonly,
  int order,
  bool order_null,
  bool order_max,
  size_t obs_req,
  bool min_obs_stop,
  ENUM_NSS_PART_NOISE_REDUCTION noise_reduction,
  StringVector &RP_quadrant, 
  StringVector &RP_prior_quadrant,
  NumericVector &RP_x, 
  NumericVector &RP_y
) {
	// noise reduction functions
	Function RFunc_gravity("gravity"),
			 RFunc_mean("mean"),
			 RFunc_median("median"),
			 RFunc_mode("mode"),
			 RFunc_gravity_class("gravity_class"),
			 RFunc_mode_class("mode_class");
	Rcpp::String s1("1"), s2("2"), s3("3"), s4("4");

	// count unsorted maps and max/min counters, and greater than required observation vector map
	std::unordered_map<Rcpp::String, size_t> q_counts{}, pq_counts{};
	size_t max_q_counts=0, max_pq_counts=0, min_q_counts=0, min_pq_counts=0;
	size_t x_size=x.size(), y_size=y.size();
	size_t obs_req_2=obs_req/2;
	StringVector quadrants_gt_obs_req(0), 
				 prior_quadrants_gt_obs_req(0);

	
	size_t hard_stop=x_size;
	if (x_size != y_size){
		Rcpp::stop("varible x length != y length != quadrant length != prior_quadrant length != counts legnth != old_counts length");
		return 0;
	}
	// initialize quadrant and prior quadrant
	if (x_size != (size_t)quadrant.size())
		quadrant = StringVector(x_size);
	if (x_size != (size_t)prior_quadrant.size())
		prior_quadrant = StringVector(x_size);
	quadrant.fill("q");
	prior_quadrant.fill("pq");
	//
	if (x_size <= 8) {
		if(order_null){
			order = 1;
			order_null = false;
			order_max = false;
			hard_stop = ceil(log(x_size));
			if(hard_stop<1)
				hard_stop=1;
		} else {
			obs_req = 0;
			hard_stop = x_size;
		}
	}
	if(order_null){
		order = ceil(log(x_size));
		if(order<1)
			order=1;
		order_max = false;
	}
	if(!order_max){
		obs_req = 0;
		hard_stop = ceil(log(x_size)) + 2;
		if(hard_stop<3)
			hard_stop=3;
	} else {
		hard_stop = 2*ceil(log(x_size)) + 2;
		if(hard_stop<4)
			hard_stop=4;
	}
	obs_req_2 = obs_req/2;
	//
	size_t i=0;
	std::unordered_map<std::string, size_t> reverse_RP_prior_quadrant{};
	Rcout << "order=" << order << ", hard_stop="<< hard_stop << ", obs_req=" << obs_req << "\n";
	while (true) {
		Rcout << i << "\n";
		if (i == (size_t)order || i == (size_t)hard_stop)
			break;
		// START COUNT QUADRANTS
		// create unsorted maps of quadrant=>count
		q_counts.clear();
		pq_counts.clear();
		max_q_counts=0, max_pq_counts=0, min_q_counts=0, min_pq_counts=0;
		for(size_t m=0;m<x_size;m++){
			Rcpp::String _q =quadrant[m], _pq=prior_quadrant[m];
			if (q_counts.count(_q)==0)	q_counts[_q] = 1;
			else						q_counts[_q] += 1;
			
			if (pq_counts.count(_pq)==0)pq_counts[_pq] = 1;
			else						pq_counts[_pq] += 1;
			
			if(min_q_counts==0 || q_counts[_q]<min_q_counts)	min_q_counts=q_counts[_q];
			if(q_counts[_q] > max_q_counts)						max_q_counts=q_counts[_q];
			
			if(min_pq_counts==0 || pq_counts[_pq]<min_pq_counts)min_pq_counts=q_counts[_pq];
			if(pq_counts[_pq] > max_pq_counts)					max_pq_counts=pq_counts[_pq];
		}
		
		// COUNT greater than required observations
		quadrants_gt_obs_req = StringVector(0);
		prior_quadrants_gt_obs_req = StringVector(0);
		size_t _obs_r = (type_xonly ? obs_req_2 : obs_req);
		for (auto& it: q_counts) {
			if (it.second >= _obs_r)
				quadrants_gt_obs_req.push_back(it.first);
		}
		for (auto& it: pq_counts) {
			if (it.second >= _obs_r)
				prior_quadrants_gt_obs_req.push_back(it.first);
		}

		// END COUNT QUADRANTS
		
		
		// max_q_counts                      == R"l.PART <- max(PART$counts)"
		// quadrants_gt_obs_req.size()       == "obs.req.rows <- PART[counts >= obs.req, which = TRUE]"
		// prior_quadrants_gt_obs_req.size() == "old.obs.req.rows <- PART[old.counts >= obs.req, which = TRUE]"
		if(!type_xonly){
			if (
				(quadrants_gt_obs_req.size() <= 0) ||
				(min_obs_stop && obs_req > 0 && quadrants_gt_obs_req.size() < prior_quadrants_gt_obs_req.size())
			)
				break;
		}else{
			// X ONLY
			if (
				(quadrants_gt_obs_req.size() <= 0) ||
				(obs_req > 0 && quadrants_gt_obs_req.size() < prior_quadrants_gt_obs_req.size())
			)
				break;
		}

		
		LogicalVector obs_req_rows = in(quadrant, quadrants_gt_obs_req);
		//LogicalVector old_obs_req_rows = in(prior_quadrant, prior_quadrants_gt_obs_req);


		// NOISE REDUCTION
		StringVector _tmp_slice = quadrant[obs_req_rows];
		RP_quadrant = Rcpp::unique(_tmp_slice);
		size_t RP_quadrant_size = RP_quadrant.size();
		RP_x = NumericVector(RP_quadrant_size);
		RP_y = NumericVector(RP_quadrant_size);
		
		reverse_RP_prior_quadrant.clear();
		for(size_t ri=0;ri<RP_quadrant_size;ri++){
		  std::string _s = as<std::string>(RP_quadrant[ri]);
		  reverse_RP_prior_quadrant[_s]=ri;
		  LogicalVector _rows = in(quadrant, as<Rcpp::StringVector>(RP_quadrant[ri]));
		  switch(noise_reduction){
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_OFF:
			  RP_x[ri] = as<double>( RFunc_gravity(x[_rows]) );
			  RP_y[ri] = as<double>( RFunc_gravity(y[_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MEAN:
			  RP_x[ri] = as<double>( RFunc_gravity(x[_rows]) );
			  RP_y[ri] = as<double>( RFunc_mean(y[_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MEDIAN:
			  RP_x[ri] = as<double>( RFunc_gravity(x[_rows]) );
			  RP_y[ri] = as<double>( RFunc_median(y[_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MODE:
			  RP_x[ri] = as<double>( RFunc_gravity(x[_rows]) );
			  RP_y[ri] = as<double>( RFunc_mode(y[_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MODE_CLASS:
			  RP_x[ri] = as<double>( RFunc_gravity_class(x[_rows]) );
			  RP_y[ri] = as<double>( RFunc_mode_class(y[_rows]) );
			  break;
		  }
		}
		// new quadrant id
		for(size_t ri=0; ri<x_size; ri++){
			if(!obs_req_rows[ri])
				continue;
			Rcpp::String old_qs = quadrant[ri];
			std::string _s = as<std::string>(quadrant[ri]);
			size_t rev_RP=reverse_RP_prior_quadrant.at(_s);
			double cur_x = x[ri], 
				   cur_y = y[ri],
			       cur_rp_x = RP_x[rev_RP],
			       cur_rp_y = RP_y[rev_RP];
			int new_q = (
				type_xonly
				?(cur_x <= cur_rp_x ? 0 : 1)
				:(cur_x <= cur_rp_x ? 0 : 1) + (cur_y <= cur_rp_y ? 2 : 3)
			);
			prior_quadrant[ri] = old_qs;
			old_qs += (new_q==0?s1:(new_q==1?s2:(new_q==2?s3:s4)));
			quadrant[ri] = old_qs;
		}
		if((min_q_counts <= obs_req) && i >= 1)
			break;
		i++;
	}
	return i;
}
