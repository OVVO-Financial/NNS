// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "partial_moments.h"



double repeatMultiplication(double value, int n) {
  double result = 1.0;
  for (int i = 0; i < n; ++i) {
    result *= value;
  }
  return result;
}

double fastPow(double a, double b) {
  union {
  double d;
  int x[2];
} u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

// Function to check if a value is an integer
inline bool isInteger(double value) {
  return value == static_cast<int>(value);
}



/////////////////
// UPM / LPM
// single thread
double LPM_C(const double &degree, const double &target, const RVector<double> &variable) {
  size_t n = variable.size();
  double out = 0;
  double value;
  
  for (size_t i = 0; i < n; i++) {
    value = target - variable[i];
    if (value >= 0) {
      if (isInteger(degree)) {
        if (degree == 0) {
          out += 1;
        } else if (degree == 1) {
          out += value;
        } else {
          // Use repeatMultiplication function for integer degrees
          out += repeatMultiplication(value, static_cast<int>(degree));
        }
      } else {
        // Use fastPow for non-integer degrees
        out += fastPow(value, degree);
      }
    } else out+= 0;
  }
  out /= n;
  return out;
}

double UPM_C(const double &degree, const double &target, const RVector<double> &variable) {
  size_t n = variable.size();
  double out = 0;
  double value;
  
  for (size_t i = 0; i < n; i++) {
    value = variable[i] - target;
    if (value > 0) {
      if (isInteger(degree)) {
        if (degree == 0) {
          out += 1;
        } else if (degree == 1) {
          out += value;
        } else {
          // Use repeatMultiplication function for integer degrees
          out += repeatMultiplication(value, static_cast<int>(degree));
        }
      } else {
        // Use fastPow for non-integer degrees
        out += fastPow(value, degree);
      }
    } else out+= 0;
  }
  out /= n;
  return out;
}

// parallelFor
#define NNS_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS)      \
size_t target_size=target.size();                        \
NumericVector output = NumericVector(target_size);       \
WORKER_CLASS tmp_func(degree, target, variable, output); \
parallelFor(0, target_size, tmp_func);                   \
return(output);
NumericVector LPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(LPM_Worker);
}
NumericVector UPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(UPM_Worker);
}

NumericVector LPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  if (degree>0) {
    NumericVector lpm_output = LPM_CPv(degree, target, variable);
    NumericVector upm_output = UPM_CPv(degree, target, variable);
    NumericVector area = lpm_output+upm_output;
    return(lpm_output / area);
  } else {
    return LPM_CPv(degree, target, variable);
  }
}
NumericVector UPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  if (degree>0) {
    NumericVector lpm_output = LPM_CPv(degree, target, variable);
    NumericVector upm_output = UPM_CPv(degree, target, variable);
    NumericVector area = lpm_output+upm_output;
    return(upm_output / area);
  } else {
    return UPM_CPv(degree, target, variable);
  }
}

/////////////////
// CoUPM / CoLPM / DUPM / DLPM
// single thread
double CoUPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
){
  size_t n_x = x.size(), n_y = y.size();
  size_t max_size = (n_x>n_y ? n_x : n_y);
  size_t min_size = (n_x<n_y ? n_x : n_y);
  if (n_x != n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  
  double out=0;
  bool d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    
    double y1=(y[i]-target_y);
    
    if(d_upm_0){
      x1 = (x1 > 0 ? 1 : x1);
      y1 = (y1 > 0 ? 1 : y1);
    }
    
    x1 = (x1 < 0 ? 0 : x1);
    y1 = (y1 < 0 ? 0 : y1);
    
    if(isInteger(degree_upm)){
      if(d_upm_0) out += x1 * y1; 
      else
        // Use repeatMultiplication function for integer degrees
        out += repeatMultiplication(x1, static_cast<int>(degree_upm)) * repeatMultiplication(y1, static_cast<int>(degree_upm));
    } else out += fastPow(x1, degree_upm) * fastPow(y1, degree_upm);
  }
  return out/max_size;
}

double CoLPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  bool d_lpm_0=(degree_lpm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    
    double y1=(target_y-y[i]);
    
    if(d_lpm_0){
      x1 = (x1 >= 0 ? 1 : x1);
      y1 = (y1 >= 0 ? 1 : y1);
    }
    
    x1 = (x1 < 0 ? 0 : x1);
    y1 = (y1 < 0 ? 0 : y1);
    
    if(isInteger(degree_lpm)){
      if(d_lpm_0) out += x1 * y1;
      else
        // Use repeatMultiplication function for integer degrees
        out += repeatMultiplication(x1, static_cast<int>(degree_lpm)) * repeatMultiplication(y1, static_cast<int>(degree_lpm));
    } else out += fastPow(x1, degree_lpm) * fastPow(y1, degree_lpm);
  }
  return out/max_size;
}



double DLPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  bool dont_use_pow_lpm=isInteger(degree_lpm), 
    dont_use_pow_upm=isInteger(degree_upm),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    
    double y1=(target_y-y[i]);
    
    if(d_upm_0) x1 = (x1 > 0 ? 1 : x1);
    if(d_lpm_0) y1 = (y1 >= 0 ? 1 : y1);
    
    x1 = (x1 < 0 ? 0 : x1);
    y1 = (y1 < 0 ? 0 : y1);
    
    
    if(dont_use_pow_lpm && dont_use_pow_upm){
      if(!d_upm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_upm));
      if(!d_lpm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_lpm));
      out += x1 * y1;
    } else if(dont_use_pow_lpm && !dont_use_pow_upm){
      if(!d_lpm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_lpm));
      out += fastPow(x1, degree_upm) * y1;
    } else if(dont_use_pow_upm && !dont_use_pow_lpm){
      if(!d_upm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_upm));
      out += x1 * fastPow(y1, degree_lpm);
    } else out += fastPow(x1, degree_upm) * fastPow(y1, degree_lpm);
  }
  return out/max_size;
}



double DUPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  bool dont_use_pow_lpm=(isInteger(degree_lpm)), 
    dont_use_pow_upm=(isInteger(degree_upm)),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    
    double y1=(y[i]-target_y);
    
    if(d_lpm_0) x1 = (x1 >= 0 ? 1 : x1);
    if(d_upm_0) y1 = (y1 > 0 ? 1 : y1);
    
    x1 = (x1 < 0 ? 0 : x1);
    y1 = (y1 < 0 ? 0 : y1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm){
      if(!d_upm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_lpm));
      if(!d_lpm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_upm));
      out += x1 * y1;
    } else if(dont_use_pow_lpm && !dont_use_pow_upm){
      if(!d_lpm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_upm));
      out += fastPow(x1, degree_lpm) * y1;
    } else if(dont_use_pow_upm && !dont_use_pow_lpm){
      if(!d_upm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_lpm));
      out += x1 * fastPow(y1, degree_upm);
    } else out += fastPow(x1, degree_lpm) * fastPow(y1, degree_upm);
  }
  return out/max_size;
}

// parallelFor
#define NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS, LPM_DEGREE_VARIABLE, UPM_DEGREE_VARIABLE) \
size_t target_x_size=target_x.size();                                                               \
size_t target_y_size=target_y.size();                                                               \
size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);                   \
NumericVector output = NumericVector(max_target_size);                                              \
WORKER_CLASS tmp_func(LPM_DEGREE_VARIABLE, UPM_DEGREE_VARIABLE, x, y, target_x, target_y, output);  \
parallelFor(0, output.size(), tmp_func);                                                            \
return(output);
NumericVector CoLPM_CPv(
    const double &degree_lpm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoLPM_Worker, degree_lpm, degree_lpm);
}
NumericVector CoUPM_CPv(
    const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y 
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoUPM_Worker, degree_upm, degree_upm);
}
NumericVector DLPM_CPv(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(DLPM_Worker, degree_lpm, degree_upm);
}
NumericVector DUPM_CPv(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(DUPM_Worker, degree_lpm, degree_upm);
}



/////////////////
// PM MATRIX
// single thread
void PMMatrix_Cv(
    const double &degree_lpm, 
    const double &degree_upm, 
    const RMatrix<double>::Column &x, 
    const RMatrix<double>::Column &y, 
    const double &target_x,
    const double &target_y, 
    const bool &pop_adj, 
    const double &adjust,
    const size_t &rows, 
    double &coLpm,
    double &coUpm,   
    double &dLpm, 
    double &dUpm,
    double &covMat
){
  
  // Convert RMatrix<double>::Column to RVector<double>
  RVector<double> x_rvec(x);
  RVector<double> y_rvec(y);
  
  coLpm=CoLPM_C(degree_lpm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  coUpm=CoUPM_C(degree_lpm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  dLpm=DLPM_C(degree_lpm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  dUpm=DUPM_C(degree_lpm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  covMat=0;
  if(rows == 0)
    return;
  
  if(pop_adj && rows > 1){
    coLpm *= adjust;
    coUpm *= adjust;
    dLpm *= adjust;
    dUpm *= adjust;
  }
  covMat = coUpm + coLpm - dUpm - dLpm;
}


// parallelFor
List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const NumericVector &target,
    const NumericMatrix &variable,
    const bool &pop_adj
) {
  size_t variable_cols=variable.cols();
  size_t target_length=target.size();
  if(variable_cols != target_length){
    Rcpp::stop("varible matrix cols != target vector length");
    return List::create();
  }
  NumericMatrix coLpm(variable_cols, variable_cols);
  NumericMatrix coUpm(variable_cols, variable_cols);
  NumericMatrix dLpm(variable_cols, variable_cols);
  NumericMatrix dUpm(variable_cols, variable_cols);
  NumericMatrix covMat(variable_cols, variable_cols);
  
  PMMatrix_Worker tmp_func(LPM_degree, UPM_degree, variable, target, pop_adj, coLpm, coUpm, dLpm, dUpm, covMat);
  parallelFor(0, variable_cols, tmp_func);
  
  rownames(coLpm) = colnames(variable);
  colnames(coLpm) = colnames(variable);
  
  rownames(coUpm) = colnames(variable);
  colnames(coUpm) = colnames(variable);
  
  rownames(dLpm) = colnames(variable);
  colnames(dLpm) = colnames(variable);
  
  rownames(dUpm) = colnames(variable);
  colnames(dUpm) = colnames(variable);
  
  rownames(covMat) = colnames(variable);
  colnames(covMat) = colnames(variable);
  
  return(
    List::create(
      Named("cupm") = coUpm,
      Named("dupm") = dUpm,
      Named("dlpm") = dLpm,
      Named("clpm") = coLpm,
      Named("cov.matrix") = covMat
    )
  );
}

