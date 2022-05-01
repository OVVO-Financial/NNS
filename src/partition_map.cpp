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
  bool obs_req_null,
  bool min_obs_stop,
  ENUM_NSS_PART_NOISE_REDUCTION noise_reduction,
  bool Voronoi,
  StringVector &RP_quadrant, 
  StringVector &RP_prior_quadrant,
  NumericVector &RP_x, 
  NumericVector &RP_y
) {
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	// noise reduction functions
	Function RFunc_gravity("gravity"),
			 RFunc_mean("mean"),
			 RFunc_median("median"),
			 RFunc_mode("mode"),
			 RFunc_gravity_class("gravity_class"),
			 RFunc_mode_class("mode_class"),
			 RFunc_plot("plot"),
			 RFunc_segments("segments"),
			 RFunc_title("title"),
			 RFunc_points("points"),
			 RFunc_abline("abline"),
			 RFunc_min("min"),
			 RFunc_max("max");
	Rcpp::String s1("1"), s2("2"), s3("3"), s4("4");

	// count unsorted maps and max/min counters, and greater than required observation vector map
	std::unordered_map<Rcpp::String, size_t> q_counts{}, pq_counts{};
	size_t max_q_counts=0, min_q_counts=0;
	size_t x_size=x.size(), y_size=y.size();
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
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	if (obs_req_null)
		obs_req=8;
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	if (!order_null && order == 0 && !order_max){
//Rcout << "(!order_null && order == 0 && order_max)" << "\n";
		order = 1;
		order_max = false;
		order_null = false;
	}
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	if (x_size <= 8) {
		if(order_null){
//Rcout << "x_size<=8 order_null" << "\n";
			order = 1;
			order_null = false;
			order_max = false;
			hard_stop = ceil(log2(x_size));
			if(hard_stop<1)
				hard_stop=1;
		} else {
			obs_req = 0;
			hard_stop = x_size;
		}
	}
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	if(order_null){
//Rcout << "order_null" << "\n";
		order = ceil(log2(x_size));
		if(order<1)
			order=1;
		order_max = false;
	}
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	if (order_max) {
//Rcout << "order_max" << "\n";
		obs_req = 0;
		hard_stop = ceil(log2(x_size)) + 2;
		if(hard_stop<3)
			hard_stop=3;
	} else {
//Rcout << "order numbered " << order << " x_size=" << x_size << ", log=" << log2(x_size) << ", ceil=" << ceil(log2(x_size)) << "\n";
		hard_stop = 2*ceil(log2(x_size)) + 2;
		if(hard_stop<4)
			hard_stop=4;
	}
//Rcout << "type_xonly=" << type_xonly << ", order=" << order << ", order_null=" << order_null << ", order_max=" << order_max << ", obs_req=" << obs_req << ", obs_req_null=" << obs_req_null << ", min_obs_stop="<< min_obs_stop << ", noise_reduction=" << static_cast<int>(noise_reduction) << "\n";
	//
	size_t obs_req_2 = obs_req/2;
	size_t i=0;
	std::unordered_map<std::string, size_t> reverse_RP_quadrant{};
//Rcout << "order=" << order << ", hard_stop="<< hard_stop << ", obs_req=" << obs_req << "\n";

	if (Voronoi) {
	  RFunc_plot(
		x, y, 
		Named("col") = "steelblue",
		Named("cex.lab") = 1.5,
		Named("xlab") = "x", //Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "substitute"), "deparse"), 
		Named("ylab") = "y"  //Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "substitute"), "deparse")
	  );
    }

	while (true) {
		//Rcout << i << "\n";
		if (i == (size_t)order || i == (size_t)hard_stop)
			break;
		// START COUNT QUADRANTS
		// create unsorted maps of quadrant=>count
		q_counts.clear();
		pq_counts.clear();
		max_q_counts=0, min_q_counts=0;
		for(size_t m=0;m<x_size;m++){
			Rcpp::String _q =quadrant[m], _pq=prior_quadrant[m];
			if (q_counts.count(_q)==0)	q_counts[_q] = 1;
			else						q_counts[_q] += 1;
			
			if (pq_counts.count(_pq)==0)pq_counts[_pq] = 1;
			else						pq_counts[_pq] += 1;
			
			if(min_q_counts==0 || q_counts[_q]<min_q_counts)	min_q_counts=q_counts[_q];
			if(q_counts[_q] > max_q_counts)						max_q_counts=q_counts[_q];
			
			//if(min_pq_counts==0 || pq_counts[_pq]<min_pq_counts)min_pq_counts=q_counts[_pq];
			//if(pq_counts[_pq] > max_pq_counts)					max_pq_counts=pq_counts[_pq];
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
		//RP_quadrant = Rcpp::unique(_tmp_slice);
		RP_quadrant = Rcpp::sort_unique(_tmp_slice);
		size_t RP_quadrant_size = RP_quadrant.size();
		RP_x = NumericVector(RP_quadrant_size);
		RP_y = NumericVector(RP_quadrant_size);
		
		reverse_RP_quadrant.clear();
		for(size_t ri=0;ri<RP_quadrant_size;ri++){
		  std::string _s = as<std::string>(RP_quadrant[ri]);
		  reverse_RP_quadrant[_s]=ri;
		  LogicalVector selected_rows = in(quadrant, as<Rcpp::StringVector>(RP_quadrant[ri]));
		  switch (noise_reduction) {
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_OFF:
			  RP_x[ri] = as<double>( RFunc_gravity(x[selected_rows]) );
			  RP_y[ri] = as<double>( RFunc_gravity(y[selected_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MEAN:
			  RP_x[ri] = as<double>( RFunc_gravity(x[selected_rows]) );
			  RP_y[ri] = as<double>( RFunc_mean(y[selected_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MEDIAN:
			  RP_x[ri] = as<double>( RFunc_gravity(x[selected_rows]) );
			  RP_y[ri] = as<double>( RFunc_median(y[selected_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MODE:
			  RP_x[ri] = as<double>( RFunc_gravity(x[selected_rows]) );
			  RP_y[ri] = as<double>( RFunc_mode(y[selected_rows]) );
			  break;
			case ENUM_NSS_PART_NOISE_REDUCTION::NOISE_REDUCTION_MODE_CLASS:
			  RP_x[ri] = as<double>( RFunc_gravity_class(x[selected_rows]) );
			  RP_y[ri] = as<double>( RFunc_mode_class(y[selected_rows]) );
			  break;
		  }

		  if (Voronoi){
		    if(!type_xonly){
			  // l.part = max_q_counts
		      if(max_q_counts > obs_req){
			    //PART[obs.req.rows, {
				//    segments(min(x), gravity(y), max(x), gravity(y), lty = 3)
				//    segments(gravity(x), min(y), gravity(x), max(y), lty = 3)
  				//  }, by = quadrant]
			    double _min_x = as<double>(RFunc_min(x[selected_rows])), 
			           _min_y = as<double>(RFunc_min(y[selected_rows])),
					   _max_x = as<double>(RFunc_max(x[selected_rows])), 
				  	   _max_y = as<double>(RFunc_max(y[selected_rows]));
			    RFunc_segments(_min_x, RP_y[ri], _max_x, RP_y[ri], Named("lty") = 3);
  			    RFunc_segments(RP_x[ri], _min_y, RP_x[ri], _max_y, Named("lty") = 3);
			  }
		    }else{
			  double _min_x = as<double>(RFunc_min(x[selected_rows])), 
				     _max_x = as<double>(RFunc_max(x[selected_rows]));
			  RFunc_abline(Named("v") = NumericVector::create(_min_x, _max_x), Named("lty") = 3);
  		    }
		  }

		}
		// new quadrant id
		for(size_t ri=0; ri<x_size; ri++){
			if(!obs_req_rows[ri])
				continue;
			Rcpp::String old_qs = quadrant[ri];
			std::string _s = as<std::string>(quadrant[ri]);
			size_t rev_RP=reverse_RP_quadrant.at(_s);
			double cur_x = x[ri], 
				   cur_y = y[ri],
			       cur_rp_x = RP_x[rev_RP],
			       cur_rp_y = RP_y[rev_RP];
/*
	xa) x<=rx -> 1
	xb) x>rx  -> 0
	ya) y<=ry -> 1
	yb) y>ry  -> 0

	xb yb = 1 + 0 + 0 = 1
	xa yb = 1 + 1 + 0 = 2
	xb ya = 1 + 0 + 2 = 3
	xa ya = 1 + 1 + 2 = 4
*/
			char new_q = 1 + (cur_x <= cur_rp_x ? 1 : 0) + (type_xonly ? 0 : (cur_y <= cur_rp_y ? 2 : 0));
//Rcout << "ri=" << ri << ", q=" << old_qs.get_cstring() << ", cur_x=" << cur_x << ", cur_y=" << cur_y << ", cur_rp_x=" << cur_rp_x << ", cur_rp_y=" << cur_rp_y << ", new_q=" << new_q;
			prior_quadrant[ri] = old_qs;
			old_qs += (new_q==1?s1:(new_q==2?s2:(new_q==3?s3:s4)));
//Rcout << ", q=" << old_qs.get_cstring() << "\n";
			quadrant[ri] = old_qs;
		}
		if((min_q_counts <= obs_req) && i >= 1)
			break;
		i++;
	}

	if (type_xonly){
		//if(mean(c(length(unique(diff(x))), length(unique(x)))) < .33*length(x)){
		//  RP$x <- ifelse(RP$x%%1 < .5, floor(RP$x), ceiling(RP$x))
		//}
		double _mean = mean(
			NumericVector::create(
				unique(diff(x)).size(), 
				unique(x).size()
			)
		);
	    if (_mean < .33 * x_size){
			// https://stackoverflow.com/questions/27686319/c-armadillo-modulus-function
			// mod = a - floor(a/n)*n
		  RP_x = ifelse(RP_x - floor(RP_x) < .5, floor(RP_x), ceiling(RP_x));
		}
	}


	if (Voronoi) {
      RFunc_title(Named("main") = ("NNS Order = " + std::to_string(i)), Named("cex.main") = 2);
      if(min_obs_stop)
		RFunc_points(RP_x, RP_y, Named("pch")=15, Named("lwd")=2, Named("col") = "red");
    }

	if (!min_obs_stop) {
	  RP_x = NumericVector();
	  RP_y = NumericVector();
	  RP_quadrant = StringVector();
	  RP_prior_quadrant = StringVector();
	}
	return i;
}
