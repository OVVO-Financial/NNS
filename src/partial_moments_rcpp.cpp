// partial_moments_rcpp.cpp
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "partial_moments.h"
using namespace Rcpp;

static inline double repeatMultiplication(double value, int n) {
  double result = 1.0;
  for (int i = 0; i < n; ++i) {
    result *= value;
  }
  return result;
}

static inline double fastPow(double a, double b) {
  union { double d; int x[2]; } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

static inline bool isInteger(double value) {
  return value == static_cast<int>(value);
}


// [[Rcpp::export(rng = false)]]
double CoLPM_nD_RCPP(const NumericMatrix &data,
                     const NumericVector &target,
                     const double &degree,
                     const bool &norm = true) {
  return clpm_nD_cpp(data, target, degree, norm);
}

// [[Rcpp::export(rng = false)]]
double CoUPM_nD_RCPP(const NumericMatrix &data,
                     const NumericVector &target,
                     const double &degree,
                     const bool &norm = true) {
  return cupm_nD_cpp(data, target, degree, norm);
}

// [[Rcpp::export(rng = false)]]
double DPM_nD_RCPP(const NumericMatrix &data,
                     const NumericVector &target,
                     const double &degree,
                     const bool &norm = true) {
  return dpm_nD_cpp(data, target, degree, norm);
} 



// [[Rcpp::export(rng = false)]]
NumericVector LPM_RCPP(const double &degree,
                       const RObject &target,
                       const RObject &variable,
                       const bool &excess_ret) {
  NumericVector variable_vec = as<NumericVector>(clone(variable));
  NumericVector target_vec;
  if (is<NumericVector>(target) && !target.isNULL()) {
    target_vec = as<NumericVector>(target);
  } else {
    target_vec = NumericVector::create(mean(variable_vec));
  }
  
  if (excess_ret) {
    int n = variable_vec.size();
    int tlen = target_vec.size();
    if (!(tlen == 1 || tlen == n))
      Rcpp::stop("When excess_ret=TRUE, target must be length 1 or same length as variable");
    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
      double t    = (tlen == 1 ? target_vec[0] : target_vec[i]);
      double diff = t - variable_vec[i];
      if (diff > 0) {
        if      (degree == 0)       out[i] = 1;
        else if (degree == 1)       out[i] = diff;
        else if (isInteger(degree)) out[i] = repeatMultiplication(diff, (int)degree);
        else                         out[i] = fastPow(diff, degree);
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
  NumericVector variable_vec = as<NumericVector>(clone(variable));
  NumericVector target_vec;
  if (is<NumericVector>(target) && !target.isNULL()) {
    target_vec = as<NumericVector>(target);
  } else {
    target_vec = NumericVector::create(mean(variable_vec));
  }
  
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
        else                         out[i] = fastPow(diff, degree);
      }
    }
    return NumericVector::create(mean(out));
  }
  
  return UPM_CPv(degree, target_vec, variable_vec);
}

//' @name LPM.ratio
//' @title Lower Partial Moment Ratio
//' @description
//'   This function generates a standardized univariate lower partial moment
//'   of any non‑negative degree for a given target.
//' @param degree numeric; degree = 0 gives frequency (CDF), degree = 1 gives area.
//' @param target numeric vector; threshold(s). Defaults to mean(variable).
//' @param variable numeric vector or data‑frame column to evaluate.
//' @return Numeric vector of standardized lower partial moments.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @references
//'   Viole, F. (2017) Continuous CDFs and ANOVA with NNS. \doi{10.2139/ssrn.3007373}
//' @examples
//'   set.seed(123)
//'   x <- rnorm(100)
//'   LPM.ratio(0, mean(x), x)
//' \dontrun{
//'   plot(sort(x), LPM.ratio(0, sort(x), x))
//'   plot(sort(x), LPM.ratio(1, sort(x), x))
//' }
//' @export
// [[Rcpp::export("LPM.ratio", rng = false)]]
 NumericVector LPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable) {
   NumericVector target_vec, variable_vec;
   if (is<NumericVector>(variable))
     variable_vec=as<NumericVector>(variable);
   else if (is<IntegerVector>(variable))
     variable_vec=as<NumericVector>(variable);
   else if (is<DataFrame>(variable))
     variable_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(variable, "unlist"), "as.vector");
   else
     Rcpp::stop("variable should be numeric vector, or data table");
   if (is<NumericVector>(target) && !target.isNULL()){
     target_vec = as<NumericVector>(target);
   }else{
     target_vec = NumericVector(1);
     target_vec[0] = mean(variable_vec);
   }
   return LPM_ratio_CPv(degree, target_vec, variable_vec);
 }


//' @name UPM.ratio
//' @title Upper Partial Moment Ratio
//' @description
//'   This function generates a standardized univariate upper partial moment
//'   of any non‑negative degree for a given target.
//' @param degree numeric; degree = 0 gives frequency, degree = 1 gives area.
//' @param target numeric vector; threshold(s). Defaults to mean(variable).
//' @param variable numeric vector or data‑frame column to evaluate.
//' @return Numeric vector of standardized upper partial moments.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @examples
//'   set.seed(123)
//'   x <- rnorm(100)
//'   UPM.ratio(0, mean(x), x)
//' \dontrun{
//'   plot3d(x, y, Co.UPM(0, sort(x), sort(y), x, y), …)
//' }
//' @export
// [[Rcpp::export("UPM.ratio", rng = false)]]
 NumericVector UPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable) {
   NumericVector target_vec, variable_vec;
   if (is<NumericVector>(variable))
     variable_vec=as<NumericVector>(variable);
   else if (is<IntegerVector>(variable))
     variable_vec=as<NumericVector>(variable);
   else if (is<DataFrame>(variable))
     variable_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(variable, "unlist"), "as.vector");
   else
     Rcpp::stop("variable should be numeric vector, or data table");
   if (is<NumericVector>(target) && !target.isNULL()){
     target_vec = as<NumericVector>(target);
   }else{
     target_vec = NumericVector(1);
     target_vec[0] = mean(variable_vec);
   }
   return UPM_ratio_CPv(degree, target_vec, variable_vec);
 }


//' @name Co.LPM
//' @title Co‑Lower Partial Moment
//' @description
//'   Computes the co‑lower partial moment (lower‑left quadrant 4) between two
//'   equal‑length numeric vectors at any degree and target.
//' @param degree_lpm numeric; degree = 0 gives frequency, degree = 1 gives area.
//' @param x numeric vector of observations.
//' @param y numeric vector of the same length as x.
//' @param target_x numeric vector; thresholds for x (defaults to mean(x)).
//' @param target_y numeric vector; thresholds for y (defaults to mean(y)).
//' @return Numeric vector of co‑LPM values.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @examples
//'   set.seed(123)
//'   x <- rnorm(100); y <- rnorm(100)
//'   Co.LPM(0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("Co.LPM", rng = false)]]
 NumericVector CoLPM_RCPP(
     const double &degree_lpm, 
     const RObject &x, const RObject &y, 
     const RObject &target_x, const RObject &target_y
 ) {
   NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
   if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
   else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
   else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
   else                         Rcpp::stop("x should be numeric vector, or data table");
   
   if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
   else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
   else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
   else                         Rcpp::stop("y should be numeric vector, or data table");
   
   if (is<NumericVector>(target_x) && !target_x.isNULL()){
     target_x_vec = as<NumericVector>(target_x);
   }else{
     target_x_vec = NumericVector(1);
     target_x_vec[0] = mean(x_vec);
   }
   if (is<NumericVector>(target_y) && !target_y.isNULL()){
     target_y_vec = as<NumericVector>(target_y);
   }else{
     target_y_vec = NumericVector(1);
     target_y_vec[0] = mean(y_vec);
   }
   return CoLPM_CPv(degree_lpm, x_vec, y_vec, target_x_vec, target_y_vec);
 }


//' @name Co.UPM
//' @title Co‑Upper Partial Moment
//' @description
//'   Computes the co‑upper partial moment (upper‑right quadrant 1) between two
//'   equal‑length numeric vectors at any degree and target.
//' @param degree_upm numeric; degree = 0 gives frequency, degree = 1 gives area.
//' @param x numeric vector of observations.
//' @param y numeric vector of the same length as x.
//' @param target_x numeric vector; thresholds for x (defaults to mean(x)).
//' @param target_y numeric vector; thresholds for y (defaults to mean(y)).
//' @return Numeric vector of co‑UPM values.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @examples
//'   set.seed(123)
//'   x <- rnorm(100); y <- rnorm(100)
//'   Co.UPM(0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("Co.UPM", rng = false)]]
 NumericVector CoUPM_RCPP(
     const double &degree_upm, 
     const RObject &x, const RObject &y, 
     const RObject &target_x, const RObject &target_y
 ) {
   NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
   if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
   else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
   else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
   else                         Rcpp::stop("x should be numeric vector, or data table");
   
   if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
   else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
   else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
   else                         Rcpp::stop("y should be numeric vector, or data table");
   
   if (is<NumericVector>(target_x) && !target_x.isNULL()){
     target_x_vec = as<NumericVector>(target_x);
   }else{
     target_x_vec = NumericVector(1);
     target_x_vec[0] = mean(x_vec);
   }
   if (is<NumericVector>(target_y) && !target_y.isNULL()){
     target_y_vec = as<NumericVector>(target_y);
   }else{
     target_y_vec = NumericVector(1);
     target_y_vec[0] = mean(y_vec);
   }
   return CoUPM_CPv(degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
 }


//' @name D.LPM
//' @title Divergent‑Lower Partial Moment
//' @description
//'   Computes the divergent lower partial moment (lower‑right quadrant 3)
//'   between two equal‑length numeric vectors.
//' @param degree_lpm numeric; LPM degree = 0 gives frequency, = 1 gives area.
//' @param degree_upm numeric; UPM degree = 0 gives frequency, = 1 gives area.
//' @param x numeric vector of observations.
//' @param y numeric vector of the same length as x.
//' @param target_x numeric vector; thresholds for x (defaults to mean(x)).
//' @param target_y numeric vector; thresholds for y (defaults to mean(y)).
//' @return Numeric vector of divergent LPM values.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @examples
//'   set.seed(123)
//'   x <- rnorm(100); y <- rnorm(100)
//'   D.LPM(0, 0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("D.LPM", rng = false)]]
 NumericVector DLPM_RCPP(
     const double &degree_lpm, const double &degree_upm, 
     const RObject &x, const RObject &y, 
     const RObject &target_x, const RObject &target_y
 ) {
   NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
   if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
   else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
   else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
   else                         Rcpp::stop("x should be numeric vector, or data table");
   
   if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
   else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
   else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
   else                         Rcpp::stop("y should be numeric vector, or data table");
   
   if (is<NumericVector>(target_x) && !target_x.isNULL()){
     target_x_vec = as<NumericVector>(target_x);
   }else{
     target_x_vec = NumericVector(1);
     target_x_vec[0] = mean(x_vec);
   }
   if (is<NumericVector>(target_y) && !target_y.isNULL()){
     target_y_vec = as<NumericVector>(target_y);
   }else{
     target_y_vec = NumericVector(1);
     target_y_vec[0] = mean(y_vec);
   }
   return DLPM_CPv(degree_lpm, degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
 }


//' @name D.UPM
//' @title Divergent‑Upper Partial Moment
//' @description
//'   Computes the divergent upper partial moment (upper‑left quadrant 2)
//'   between two equal‑length numeric vectors.
//' @param degree_lpm numeric; LPM degree = 0 gives frequency, = 1 gives area.
//' @param degree_upm numeric; UPM degree = 0 gives frequency, = 1 gives area.
//' @param x numeric vector of observations.
//' @param y numeric vector of the same length as x.
//' @param target_x numeric vector; thresholds for x (defaults to mean(x)).
//' @param target_y numeric vector; thresholds for y (defaults to mean(y)).
//' @return Numeric vector of divergent UPM values.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @examples
//'   set.seed(123)
//'   x <- rnorm(100); y <- rnorm(100)
//'   D.UPM(0, 0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("D.UPM", rng = false)]]
 NumericVector DUPM_RCPP(
     const double &degree_lpm, const double &degree_upm, 
     const RObject &x, const RObject &y, 
     const RObject &target_x, const RObject &target_y
 ) {
   NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
   if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
   else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
   else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
   else                         Rcpp::stop("x should be numeric vector, or data table");
   
   if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
   else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
   else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
   else                         Rcpp::stop("y should be numeric vector, or data table");
   
   if (is<NumericVector>(target_x) && !target_x.isNULL()){
     target_x_vec = as<NumericVector>(target_x);
   }else{
     target_x_vec = NumericVector(1);
     target_x_vec[0] = mean(x_vec);
   }
   if (is<NumericVector>(target_y) && !target_y.isNULL()){
     target_y_vec = as<NumericVector>(target_y);
   }else{
     target_y_vec = NumericVector(1);
     target_y_vec[0] = mean(y_vec);
   }
   return DUPM_CPv(degree_lpm, degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
 }



//' @name PM.matrix
//' @title Partial Moment Matrix
//' @description
//'   Builds a list containing all four quadrant partial‑moment matrices
//'   (CUPM, DUPM, DLPM, CLPM) plus the overall covariance matrix.
//' @param LPM_degree numeric; lower partial moment degree (0 = freq, 1 = area).
//' @param UPM_degree numeric; upper partial moment degree (0 = freq, 1 = area).
//' @param target numeric vector; thresholds for each column (defaults to colMeans).
//' @param variable numeric matrix or data.frame.
//' @param pop_adj logical; TRUE adjusts population vs. sample moments.
//' @param norm logical; \code{FALSE} (default) if TRUE, each of the four quadrant partial-moment matrices (cupm, dupm, dlpm, clpm) is normalized cell-wise so that their sum at each position is 1. The covariance matrix is then recomputed from those normalized quadrants.
//' @return A list with elements $cupm, $dupм, $dlpm, $clpm and $cov.matrix.
//' @author Fred Viole, OVVO Financial Systems
//' @references
//'   Viole, F. & Nawrocki, D. (2013) *Nonlinear Nonparametric Statistics: Using Partial Moments* (ISBN:1490523995)
//' @examples
//'   set.seed(123)
//'   A <- cbind(rnorm(100), rnorm(100), rnorm(100))
//'   PM.matrix(1, 1, NULL, A, TRUE)
//' @export
// [[Rcpp::export("PM.matrix", rng = false)]]
 List PMMatrix_RCPP(
     const double &LPM_degree,
     const double &UPM_degree,
     const RObject &target,
     const RObject &variable,
     const bool pop_adj,
     const bool norm = false
 ) {
   if(variable.isNULL()){
     Rcpp::stop("varible can't be null");
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
   if((is<NumericVector>(target) || is<DataFrame>(target)) && !target.isNULL()){
     tgt=as<NumericVector>(target);
   }else{
     tgt=colMeans(variable_matrix);
   }
   
   size_t target_length=tgt.size();
   if(variable_cols != target_length){
     Rcpp::stop("varible matrix cols != target vector length");
     return List::create();
   }
   
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