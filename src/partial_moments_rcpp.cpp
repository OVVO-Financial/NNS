// partial_moments_rcpp.cpp
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include "partial_moments.h"
#include "partial_moments_rcpp.h"

using namespace Rcpp;

static inline double repeatMultiplication(double value, int n) {
  double result = 1.0;
  for (int i = 0; i < n; ++i) {
    result *= value;
  }
  return result;
}

//static inline double fastPow(double a, double b) {
//union { double d; int x[2]; } u = { a };
//  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
//  u.x[0] = 0;
//  return u.d;
//}

static inline bool isInteger(double value) {
  return value == static_cast<int>(value);
}

static inline NumericVector moment_target(const RObject &target, const NumericVector &sample) {
  if (target.isNULL()) return NumericVector::create(mean(sample));
  return NumericVector(target);
}

// Same contract as LPM()/UPM(): non-finite values are rejected, not counted in n.
static inline void stop_if_not_finite(const NumericVector &v, const char *what) {
  for (R_xlen_t i = 0; i < v.size(); ++i) {
    if (!R_FINITE(v[i])) Rcpp::stop("`%s` must be finite.", what);
  }
}

// Runs before parallelFor: the workers index targets modulo their length and cannot raise R errors.
static inline void check_co_inputs(const NumericVector &x, const NumericVector &y,
                                    const NumericVector &target_x, const NumericVector &target_y) {
  if (x.size() != y.size()) Rcpp::stop("`x` and `y` must have the same length.");
  if (x.size() == 0) Rcpp::stop("`x` and `y` must not be empty.");
  const R_xlen_t ntx = target_x.size(), nty = target_y.size();
  if (ntx == 0 || nty == 0) Rcpp::stop("`target_x` and `target_y` must not be empty.");
  if (ntx != nty && ntx != 1 && nty != 1)
    Rcpp::stop("`target_x` and `target_y` must have the same length, or one of them length 1.");
  stop_if_not_finite(x, "x");
  stop_if_not_finite(y, "y");
  stop_if_not_finite(target_x, "target_x");
  stop_if_not_finite(target_y, "target_y");
}

// For degree > 0 the ratio is 0/0 when every observation equals the target;
// use the degree -> 0 limit (LPM 1, UPM 0, the ECDF convention x <= target).
static inline NumericVector fill_degenerate_ratio(NumericVector out,
                                                  const double &degree,
                                                  const NumericVector &target,
                                                  const NumericVector &variable,
                                                  const bool lower) {
  if (!(degree > 0) || variable.size() == 0) return out;
  const double c = variable[0];
  for (R_xlen_t j = 1; j < variable.size(); ++j) {
    if (variable[j] != c) return out;
  }
  for (R_xlen_t i = 0; i < out.size(); ++i) {
    if (ISNAN(out[i]) && target[i] == c) out[i] = lower ? 1.0 : 0.0;
  }
  return out;
}


// [[Rcpp::export(rng = false)]]
double CoLPM_nD_RCPP(const NumericMatrix &data,
                     const NumericVector &target,
                     const double &degree,
                     const bool &norm ) {
  nns_check_degree(degree, "degree");
  return clpm_nD_cpp(data, target, degree, norm);
}

// [[Rcpp::export(rng = false)]]
double CoUPM_nD_RCPP(const NumericMatrix &data,
                     const NumericVector &target,
                     const double &degree,
                     const bool &norm ) {
  nns_check_degree(degree, "degree");
  return cupm_nD_cpp(data, target, degree, norm);
}

// [[Rcpp::export(rng = false)]]
double DPM_nD_RCPP(const NumericMatrix &data,
                     const NumericVector &target,
                     const double &degree,
                     const bool &norm ) {
  nns_check_degree(degree, "degree");
  return dpm_nD_cpp(data, target, degree, norm);
} 



// [[Rcpp::export(rng = false)]]
NumericVector LPM_RCPP(const double &degree,
                       const RObject &target,
                       const RObject &variable,
                       const bool &excess_ret) {
  nns_check_degree(degree, "degree");
  NumericVector variable_vec(variable);
  NumericVector target_vec = moment_target(target, variable_vec);

  if (excess_ret) {
    int n = variable_vec.size();
    int tlen = target_vec.size();
    if (!(tlen == 1 || tlen == n))
      Rcpp::stop("When excess_ret=TRUE, target must be length 1 or same length as variable");
    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
      double t    = (tlen == 1 ? target_vec[0] : target_vec[i]);
      double diff = t - variable_vec[i];
      // Degree 0 counts ties (x <= target), as the non-excess path does.
      if (diff > 0 || (degree == 0 && diff == 0)) {
        if      (degree == 0)       out[i] = 1;
        else if (degree == 1)       out[i] = diff;
        else if (isInteger(degree)) out[i] = repeatMultiplication(diff, (int)degree);
        else                         out[i] = std::pow(diff, degree);
      }
    }
    return NumericVector::create(mean(out));
  }
  
  return LPM_CPv(degree, target_vec, variable_vec);
}

// [[Rcpp::export(rng = false)]]
NumericVector UPM_RCPP(const double &degree,
                       const RObject &target,
                       const RObject &variable,
                       const bool &excess_ret) {
  nns_check_degree(degree, "degree");
  NumericVector variable_vec(variable);
  NumericVector target_vec = moment_target(target, variable_vec);

  if (excess_ret) {
    int n = variable_vec.size();
    int tlen = target_vec.size();
    if (!(tlen == 1 || tlen == n))
      Rcpp::stop("When excess_ret=TRUE, target must be length 1 or same length as variable");
    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
      double t    = (tlen == 1 ? target_vec[0] : target_vec[i]);
      double diff = variable_vec[i] - t;
      if (diff > 0) {
        if      (degree == 0)       out[i] = 1;
        else if (degree == 1)       out[i] = diff;
        else if (isInteger(degree)) out[i] = repeatMultiplication(diff, (int)degree);
        else                         out[i] = std::pow(diff, degree);
      }
    }
    return NumericVector::create(mean(out));
  }
  
  return UPM_CPv(degree, target_vec, variable_vec);
}

// [[Rcpp::export(rng = false)]]
 NumericVector LPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable) {
   nns_check_degree(degree, "degree");
   NumericVector variable_vec(variable);
   NumericVector target_vec = moment_target(target, variable_vec);
   stop_if_not_finite(variable_vec, "variable");
   stop_if_not_finite(target_vec, "target");
   return fill_degenerate_ratio(LPM_ratio_CPv(degree, target_vec, variable_vec),
                                degree, target_vec, variable_vec, true);
 }


// [[Rcpp::export(rng = false)]]
 NumericVector UPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable) {
   nns_check_degree(degree, "degree");
   NumericVector variable_vec(variable);
   NumericVector target_vec = moment_target(target, variable_vec);
   stop_if_not_finite(variable_vec, "variable");
   stop_if_not_finite(target_vec, "target");
   return fill_degenerate_ratio(UPM_ratio_CPv(degree, target_vec, variable_vec),
                                degree, target_vec, variable_vec, false);
 }



// [[Rcpp::export(rng = false)]]
NumericVector CoLPM_RCPP(
     const double &degree_lpm,
     const RObject &x, const RObject &y,
     const RObject &target_x, const RObject &target_y,
     const double &degree_y
 ) {
   nns_check_degree(degree_lpm, "degree_lpm");
   nns_check_degree(degree_y, "degree_y");
   NumericVector x_vec(x), y_vec(y);
   NumericVector target_x_vec = moment_target(target_x, x_vec);
   NumericVector target_y_vec = moment_target(target_y, y_vec);
   check_co_inputs(x_vec, y_vec, target_x_vec, target_y_vec);
   return CoLPM_CPv(degree_lpm, degree_y, x_vec, y_vec, target_x_vec, target_y_vec);
 }



// [[Rcpp::export(rng = false)]]
NumericVector CoUPM_RCPP(
     const double &degree_upm,
     const RObject &x, const RObject &y,
     const RObject &target_x, const RObject &target_y,
     const double &degree_y
 ) {
   nns_check_degree(degree_upm, "degree_upm");
   nns_check_degree(degree_y, "degree_y");
   NumericVector x_vec(x), y_vec(y);
   NumericVector target_x_vec = moment_target(target_x, x_vec);
   NumericVector target_y_vec = moment_target(target_y, y_vec);
   check_co_inputs(x_vec, y_vec, target_x_vec, target_y_vec);
   return CoUPM_CPv(degree_upm, degree_y, x_vec, y_vec, target_x_vec, target_y_vec);
 }


// [[Rcpp::export(rng = false)]]
 NumericVector DLPM_RCPP(
     const double &degree_lpm, const double &degree_upm, 
     const RObject &x, const RObject &y, 
     const RObject &target_x, const RObject &target_y
 ) {
   nns_check_degree(degree_lpm, "degree_lpm");
   nns_check_degree(degree_upm, "degree_upm");
   NumericVector x_vec(x), y_vec(y);
   NumericVector target_x_vec = moment_target(target_x, x_vec);
   NumericVector target_y_vec = moment_target(target_y, y_vec);
   check_co_inputs(x_vec, y_vec, target_x_vec, target_y_vec);
   return DLPM_CPv(degree_lpm, degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
 }


// [[Rcpp::export(rng = false)]]
 NumericVector DUPM_RCPP(
     const double &degree_lpm, const double &degree_upm, 
     const RObject &x, const RObject &y, 
     const RObject &target_x, const RObject &target_y
 ) {
   nns_check_degree(degree_lpm, "degree_lpm");
   nns_check_degree(degree_upm, "degree_upm");
   NumericVector x_vec(x), y_vec(y);
   NumericVector target_x_vec = moment_target(target_x, x_vec);
   NumericVector target_y_vec = moment_target(target_y, y_vec);
   check_co_inputs(x_vec, y_vec, target_x_vec, target_y_vec);
   return DUPM_CPv(degree_lpm, degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
 }




// [[Rcpp::export(rng = false)]]
 List PMMatrix_RCPP(
     const double &LPM_degree,
     const double &UPM_degree,
     const RObject &target,
     const RObject &variable,
     const bool pop_adj,
     const bool norm
 ) {
   nns_check_degree(LPM_degree, "LPM_degree");
   nns_check_degree(UPM_degree, "UPM_degree");
   if(variable.isNULL()){
     Rcpp::stop("variable can't be null");
     return List::create();
   }
   NumericMatrix variable_matrix;
   if (is<NumericMatrix>(variable))
     variable_matrix = as<NumericMatrix>(variable);
   else if (is<IntegerMatrix>(variable))
     variable_matrix = as<NumericMatrix>(variable);
   else
     variable_matrix = Rcpp::internal::convert_using_rfunction(variable, "as.matrix");
   
   size_t variable_cols=variable_matrix.cols();
   NumericVector tgt;
   if(!target.isNULL()){
     tgt=as<NumericVector>(target);
   }else{
     tgt = NumericVector(variable_cols);
     for (size_t j = 0; j < variable_cols; ++j) {
       tgt[j] = mean(variable_matrix(_, j));
     }
   }
   
   size_t target_length=tgt.size();
   if(variable_cols != target_length){
     Rcpp::stop("variable matrix cols != target vector length");
     return List::create();
   }
   
   stop_if_not_finite(variable_matrix, "variable");
   stop_if_not_finite(tgt, "target");
   return PMMatrix_CPv(LPM_degree, UPM_degree, tgt, variable_matrix, pop_adj, norm);
 }



// [[Rcpp::export]]
List NNS_bin(NumericVector x, double width, double origin = 0, bool missinglast = false) {
  int bin, nmissing = 0;
  std::vector<int> out;
  
  if (width <= 0)
    stop("width must be positive");
  
  NumericVector::iterator x_it = x.begin();
  for (; x_it != x.end(); ++x_it) {
    double val = *x_it;
    if (ISNAN(val)) {
      ++nmissing;
    } else {
      if (val < origin)
        continue;
      
      bin = (val - origin) / width;
      
      if ((long long unsigned) bin >= out.size()) {
        out.resize(bin + 1);
      }
      ++out[bin];
    }
  }
  
  if (missinglast)
    out.push_back(nmissing);
  
  Rcpp::List RVAL = Rcpp::List::create(Rcpp::Named("counts") = out,
                                       Rcpp::Named("origin") = origin,
                                       Rcpp::Named("width") = width,
                                       Rcpp::Named("missing") = nmissing,
                                       Rcpp::Named("last_bin_is_missing") = missinglast);
  
  return RVAL;
}
