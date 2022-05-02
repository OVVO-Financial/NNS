#ifndef NNS_partial_moments_RCPP_H
#define NNS_partial_moments_RCPP_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector LPM_RCPP(const double &degree, const RObject &target, const RObject &variable);
NumericVector UPM_RCPP(const double &degree, const RObject &target, const RObject &variable);
NumericVector LPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable);
NumericVector UPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable);
NumericVector CoLPM_RCPP(const double &degree_lpm, const RObject &x, const RObject &y, const RObject &target_x, const RObject &target_y);
NumericVector CoUPM_RCPP(const double &degree_upm, const RObject &x, const RObject &y, const RObject &target_x, const RObject &target_y);
NumericVector DLPM_RCPP(const double &degree_lpm, const double &degree_upm, const RObject &x, const RObject &y, const RObject &target_x, const RObject &target_y);
NumericVector DUPM_RCPP(const double &degree_lpm, const double &degree_upm, const RObject &x, const RObject &y, const RObject &target_x, const RObject &target_y);
List PMMatrix_RCPP(const double &LPM_degree, const double &UPM_degree,const RObject &target, const RObject &variable,const bool pop_adj=false);

#endif // NNS_partial_moments_RCPP_H
