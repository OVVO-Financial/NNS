#ifndef NNS_partial_moments_RCPP_H
#define NNS_partial_moments_RCPP_H

#include <Rcpp.h>

// The declarations below describe the R-facing wrappers defined in
// partial_moments_rcpp.cpp. Any user-facing defaults are supplied in the R
// layer (see R/partial_moments.R) while the compiled entry points expose the
// fully expanded signatures that RcppExports.cpp expects when registering the
// native routines.

Rcpp::NumericVector LPM_RCPP(const double &degree,
                             const Rcpp::RObject &target,
                             const Rcpp::RObject &variable,
                             const bool &excess_ret);

Rcpp::NumericVector UPM_RCPP(const double &degree,
                             const Rcpp::RObject &target,
                             const Rcpp::RObject &variable,
                             const bool &excess_ret);

Rcpp::NumericVector LPM_ratio_RCPP(const double &degree,
                                   const Rcpp::RObject &target,
                                   const Rcpp::RObject &variable);
Rcpp::NumericVector UPM_ratio_RCPP(const double &degree,
                                   const Rcpp::RObject &target,
                                   const Rcpp::RObject &variable);
Rcpp::NumericVector CoLPM_RCPP(const double &degree_lpm,
                               const Rcpp::RObject &x,
                               const Rcpp::RObject &y,
                               const Rcpp::RObject &target_x,
                               const Rcpp::RObject &target_y,
                               const double &degree_y);
Rcpp::NumericVector CoUPM_RCPP(const double &degree_upm,
                               const Rcpp::RObject &x,
                               const Rcpp::RObject &y,
                               const Rcpp::RObject &target_x,
                               const Rcpp::RObject &target_y,
                               const double &degree_y);
Rcpp::NumericVector DLPM_RCPP(const double &degree_lpm,
                              const double &degree_upm,
                              const Rcpp::RObject &x,
                              const Rcpp::RObject &y,
                              const Rcpp::RObject &target_x,
                              const Rcpp::RObject &target_y);
Rcpp::NumericVector DUPM_RCPP(const double &degree_lpm,
                              const double &degree_upm,
                              const Rcpp::RObject &x,
                              const Rcpp::RObject &y,
                              const Rcpp::RObject &target_x,
                              const Rcpp::RObject &target_y);

Rcpp::List PMMatrix_RCPP(const double &LPM_degree,
                         const double &UPM_degree,
                         const Rcpp::RObject &target,
                         const Rcpp::RObject &variable,
                         const bool pop_adj,
                         const bool norm);

double DPM_nD_RCPP(const Rcpp::NumericMatrix &data,
                   const Rcpp::NumericVector &target,
                   const double &degree,
                   const bool &norm);

// n-D exported wrappers (declare explicitly for clarity)
double CoLPM_nD_RCPP(const Rcpp::NumericMatrix &data,
                     const Rcpp::NumericVector &target,
                     const double &degree,
                     const bool &norm);

double CoUPM_nD_RCPP(const Rcpp::NumericMatrix &data,
                     const Rcpp::NumericVector &target,
                     const double &degree,
                     const bool &norm);

#endif // NNS_partial_moments_RCPP_H
