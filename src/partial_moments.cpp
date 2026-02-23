// partial_moments.cpp
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include "partial_moments.h"

using namespace Rcpp;
using namespace RcppParallel;

static double repeatMultiplication(double value, int n) {
  double result = 1.0;
  for (int i = 0; i < n; ++i) result *= value;
  return result;
}

//static double fastPow(double a, double b) {
//  union { double d; int x[2]; } u = { a };
//  u.x[1] = static_cast<int>(b * (u.x[1] - 1072632447) + 1072632447);
//  u.x[0] = 0;
//  return u.d;
//}

inline bool isInteger(double v) {
  return v == static_cast<int>(v);
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
        out += std::pow(value, degree);
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
        out += std::pow(value, degree);
      }
    } else out+= 0;
  }
  out /= n;
  return out;
}


// Lower Partial Moment (LPM) count: degree == 0
struct CoLPM_CountWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  RVector<double> output;
  CoLPM_CountWorker(const NumericMatrix& data_, const NumericVector& target_, NumericVector& output_)
    : data(data_), target(target_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      bool below_all = true;
      for (std::size_t j = 0; j < d; ++j) {
        if (data(i, j) > target[j]) { below_all = false; break; }
      }
      output[i] = below_all ? 1.0 : 0.0;
    }
  }
};

// Lower Partial Moment (LPM) sum: degree > 0
struct CoLPM_SumWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  const double degree;
  RVector<double> output;
  CoLPM_SumWorker(const NumericMatrix& data_, const NumericVector& target_, double degree_, NumericVector& output_)
    : data(data_), target(target_), degree(degree_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      double prod = 1.0;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = target[j] - data(i, j);
        if (diff < 0.0) { prod = 0.0; break; }
        prod *= isInteger(degree)
          ? repeatMultiplication(diff, static_cast<int>(degree))
            : std::pow(diff, degree);
      }
      output[i] = prod;
    }
  }
};

// Upper Partial Moment (UPM) count: degree == 0
struct CoUPM_CountWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  RVector<double> output;
  CoUPM_CountWorker(const NumericMatrix& data_, const NumericVector& target_, NumericVector& output_)
    : data(data_), target(target_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      bool above_all = true;
      for (std::size_t j = 0; j < d; ++j) {
        if (data(i, j) < target[j]) { above_all = false; break; }
      }
      output[i] = above_all ? 1.0 : 0.0;
    }
  }
};

// Upper Partial Moment (UPM) sum: degree > 0
struct CoUPM_SumWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  const double degree;
  RVector<double> output;
  CoUPM_SumWorker(const NumericMatrix& data_, const NumericVector& target_, double degree_, NumericVector& output_)
    : data(data_), target(target_), degree(degree_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      double prod = 1.0;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = data(i, j) - target[j];
        if (diff < 0.0) { prod = 0.0; break; }
        prod *= isInteger(degree)
          ? repeatMultiplication(diff, static_cast<int>(degree))
            : std::pow(diff, degree);
      }
      output[i] = prod;
    }
  }
};

// Discordant Partial Moment (DPM) count: degree == 0
struct DpmCountWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  RVector<double> output;
  DpmCountWorker(const NumericMatrix& data_, const NumericVector& target_, NumericVector& output_)
    : data(data_), target(target_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      bool allBelow = true, allAbove = true;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = data(i, j) - target[j];
        if (diff >= 0.0) allBelow = false;
        if (diff <= 0.0) allAbove = false;
        if (!allBelow && !allAbove) break;
      }
      output[i] = (!allBelow && !allAbove) ? 1.0 : 0.0;
    }
  }
};

// Discordant Partial Moment (DPM) sum: degree > 0
struct DpmSumWorker : public Worker {
  const RMatrix<double> data;
  const RVector<double> target;
  const double degree;
  RVector<double> output;
  DpmSumWorker(const NumericMatrix& data_, const NumericVector& target_, double degree_, NumericVector& output_)
    : data(data_), target(target_), degree(degree_), output(output_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    std::size_t d = target.length();
    for (std::size_t i = begin; i < end; ++i) {
      bool allBelow = true, allAbove = true;
      for (std::size_t j = 0; j < d; ++j) {
        double diff = data(i, j) - target[j];
        if (diff >= 0.0) allBelow = false;
        if (diff <= 0.0) allAbove = false;
        if (!allBelow && !allAbove) break;
      }
      if (allBelow || allAbove) { output[i] = 0.0; continue; }
      double prod = 1.0;
      for (std::size_t j = 0; j < d; ++j) {
        double abs_dev = std::abs(data(i, j) - target[j]);
        prod *= isInteger(degree)
          ? repeatMultiplication(abs_dev, static_cast<int>(degree))
            : std::pow(abs_dev, degree);
      }
      output[i] = prod;
    }
  }
};

double clpm_nD_cpp(const NumericMatrix& data,
                   const NumericVector& target,
                   double degree,
                   bool norm) {
  size_t n = data.nrow();
  size_t d = data.ncol();
  if (static_cast<size_t>(target.size()) != d)
    stop("`target` length must match number of columns in `data`");
  
  if (degree == 0.0) {
    NumericVector counts(n);
    CoLPM_CountWorker countWorker(data, target, counts);
    parallelFor(0, n, countWorker);
    return sum(counts) / double(n);
  }
  
  NumericVector vals(n);
  CoLPM_SumWorker sumWorker(data, target, degree, vals);
  parallelFor(0, n, sumWorker);
  double clpm_un = sum(vals) / double(n);
  double result = clpm_un;
  
  if (norm) {
    double cupm_un = cupm_nD_cpp(data, target, degree, false);
    double dpm_un  = dpm_nD_cpp(data, target, degree, false);
    double norm_const = clpm_un + cupm_un + dpm_un;
    result = norm_const > 0.0 ? (clpm_un / norm_const) : 0.0;
  }
  return result;
}

double cupm_nD_cpp(const NumericMatrix& data,
                   const NumericVector& target,
                   double degree,
                   bool norm) {
  size_t n = data.nrow();
  size_t d = data.ncol();
  if (static_cast<size_t>(target.size()) != d)
    stop("`target` length must match number of columns in `data`");
  
  if (degree == 0.0) {
    NumericVector counts(n);
    CoUPM_CountWorker countWorker(data, target, counts);
    parallelFor(0, n, countWorker);
    return sum(counts) / double(n);
  }
  
  NumericVector vals(n);
  CoUPM_SumWorker sumWorker(data, target, degree, vals);
  parallelFor(0, n, sumWorker);
  double cupm_un = sum(vals) / double(n);
  double result = cupm_un;
  
  if (norm) {
    double clpm_un = clpm_nD_cpp(data, target, degree, false);
    double dpm_un  = dpm_nD_cpp(data, target, degree, false);
    double norm_const = clpm_un + cupm_un + dpm_un;
    result = norm_const > 0.0 ? (cupm_un / norm_const) : 0.0;
  }
  return result;
}

double dpm_nD_cpp(const NumericMatrix& data,
                  const NumericVector& target,
                  double degree,
                  bool norm) {
  size_t n = data.nrow();
  size_t d = data.ncol();
  if (static_cast<size_t>(target.size()) != d)
    stop("`target` length must match number of columns in `data`");
  
  if (degree == 0.0) {
    NumericVector counts(n);
    DpmCountWorker countWorker(data, target, counts);
    parallelFor(0, n, countWorker);
    return sum(counts) / double(n);
  }
  
  NumericVector vals(n);
  DpmSumWorker sumWorker(data, target, degree, vals);
  parallelFor(0, n, sumWorker);
  double dpm_un = sum(vals) / double(n);
  double result = dpm_un;
  
  if (norm) {
    double clpm_un = clpm_nD_cpp(data, target, degree, false);
    double cupm_un = cupm_nD_cpp(data, target, degree, false);
    double norm_const = clpm_un + cupm_un + dpm_un;
    result = norm_const > 0.0 ? (dpm_un / norm_const) : 0.0;
  }
  return result;
}

// parallelFor
#define NNS_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS)      \
size_t target_size=target.size();                        \
NumericVector output = NumericVector(target_size);       \
WORKER_CLASS tmp_func(degree, target, variable, output); \
parallelFor(0, target_size, tmp_func);                   \
return(output);                                         

// [[Rcpp::export]]
NumericVector LPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(LPM_Worker);
}

// [[Rcpp::export]]
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
  if (min_size<=0)
    return 0;
  
  double out=0;
  bool d_x_0=(degree_lpm==0);
  bool d_y_0=(degree_upm==0);
  bool x_is_int=isInteger(degree_lpm);
  bool y_is_int=isInteger(degree_upm);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    double y1=(y[i]-target_y);
    
    if(d_x_0) x1 = (x1 > 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);

    if(d_y_0) y1 = (y1 > 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);

    if(!d_x_0){
      if(x_is_int) x1 = repeatMultiplication(x1, static_cast<int>(degree_lpm));
      else x1 = std::pow(x1, degree_lpm);
    }
    if(!d_y_0){
      if(y_is_int) y1 = repeatMultiplication(y1, static_cast<int>(degree_upm));
      else y1 = std::pow(y1, degree_upm);
    }
    out += x1 * y1;
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
  if (min_size<=0)
    return 0;
  double out=0;
  bool d_x_0=(degree_lpm==0);
  bool d_y_0=(degree_upm==0);
  bool x_is_int=isInteger(degree_lpm);
  bool y_is_int=isInteger(degree_upm);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    double y1=(target_y-y[i]);
    
    if(d_x_0) x1 = (x1 >= 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);

    if(d_y_0) y1 = (y1 >= 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);

    if(!d_x_0){
      if(x_is_int) x1 = repeatMultiplication(x1, static_cast<int>(degree_lpm));
      else x1 = std::pow(x1, degree_lpm);
    }
    if(!d_y_0){
      if(y_is_int) y1 = repeatMultiplication(y1, static_cast<int>(degree_upm));
      else y1 = std::pow(y1, degree_upm);
    }
    out += x1 * y1;
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
  if (min_size<=0)
    return 0;
  double out=0;
  bool dont_use_pow_lpm=isInteger(degree_lpm), 
    dont_use_pow_upm=isInteger(degree_upm),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    double y1=(target_y-y[i]);
    
    if(d_upm_0) x1 = (x1 > 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_lpm_0) y1 = (y1 >= 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm){
      if(!d_upm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_upm));
      if(!d_lpm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_lpm));
      out += x1 * y1;
    } else if(dont_use_pow_lpm && !dont_use_pow_upm){
      if(!d_lpm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_lpm));
      out += std::pow(x1, degree_upm) * y1;
    } else if(dont_use_pow_upm && !dont_use_pow_lpm){
      if(!d_upm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_upm));
      out += x1 * std::pow(y1, degree_lpm);
    } else out += std::pow(x1, degree_upm) * std::pow(y1, degree_lpm);
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
  if (min_size<=0)
    return 0;
  double out=0;
  
  bool dont_use_pow_lpm=(isInteger(degree_lpm)), 
    dont_use_pow_upm=(isInteger(degree_upm)),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    double y1=(y[i]-target_y);
    
    if(d_lpm_0) x1 = (x1 >= 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_upm_0) y1 = (y1 > 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm){
      if(!d_lpm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_lpm));
      if(!d_upm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_upm));
      out += x1 * y1;
    } else if(dont_use_pow_lpm && !dont_use_pow_upm){
      if(!d_upm_0) y1 = repeatMultiplication(y1, static_cast<int>(degree_upm));
      out += std::pow(x1, degree_lpm) * y1;
    } else if(dont_use_pow_upm && !dont_use_pow_lpm){
      if(!d_lpm_0) x1 = repeatMultiplication(x1, static_cast<int>(degree_lpm));
      out += x1 * std::pow(y1, degree_upm);
    } else out += std::pow(x1, degree_lpm) * std::pow(y1, degree_upm);
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
    const double &degree_x, const double &degree_y,
    const NumericVector &x, const NumericVector &y,
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoLPM_Worker, degree_x, degree_y);
}
NumericVector CoUPM_CPv(
    const double &degree_x, const double &degree_y,
    const NumericVector &x, const NumericVector &y,
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoUPM_Worker, degree_x, degree_y);
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
  RVector<double> x_rvec(x);
  RVector<double> y_rvec(y);
  
  coLpm=CoLPM_C(degree_lpm, degree_lpm, x_rvec, y_rvec, target_x, target_y);
  coUpm=CoUPM_C(degree_upm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  dLpm=DLPM_C(degree_lpm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  dUpm=DUPM_C(degree_lpm, degree_upm, x_rvec, y_rvec, target_x, target_y);
  covMat=0;
  if(rows == 0)
    return;
  
  if(pop_adj && rows > 1 && degree_lpm > 0 && degree_upm > 0){
    coLpm *= adjust;
    coUpm *= adjust;
    dLpm *= adjust;
    dUpm *= adjust;
  }
  covMat = coUpm + coLpm - dUpm - dLpm;
}

List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const NumericVector &target,
    const NumericMatrix &variable,
    const bool &pop_adj,
    const bool &norm
) {
  size_t variable_cols=variable.cols();
  size_t target_length=target.size();
  if(variable_cols != target_length){
    Rcpp::stop("variable matrix cols != target vector length");
    return List::create();
  }
  NumericMatrix coLpm(variable_cols, variable_cols);
  NumericMatrix coUpm(variable_cols, variable_cols);
  NumericMatrix dLpm(variable_cols, variable_cols);
  NumericMatrix dUpm(variable_cols, variable_cols);
  NumericMatrix covMat(variable_cols, variable_cols);
  
  PMMatrix_Worker tmp_func(LPM_degree, UPM_degree, variable, target, pop_adj, coLpm, coUpm, dLpm, dUpm, covMat);
  parallelFor(0, variable_cols, tmp_func);
  
  if (norm) {
    // Normalize each quadrant matrix cell-wise so that at each (i,j):
    // cupm + dupm + dlpm + clpm = 1 (if their sum > 0), else leave as zeros.
    for (size_t i = 0; i < variable_cols; ++i) {
      for (size_t j = 0; j < variable_cols; ++j) {
        double cupm_ij = coUpm(i, j);
        double dupm_ij = dUpm(i, j);
        double dlpm_ij = dLpm(i, j);
        double clpm_ij = coLpm(i, j);
        double total = cupm_ij + dupm_ij + dlpm_ij + clpm_ij;
        if (total > 0.0) {
          coUpm(i, j) = cupm_ij / total;
          dUpm(i, j) = dupm_ij / total;
          dLpm(i, j) = dlpm_ij / total;
          coLpm(i, j) = clpm_ij / total;
        } else {
          coUpm(i, j) = dUpm(i, j) = dLpm(i, j) = coLpm(i, j) = 0.0;
        }
      }
    }
  }
  
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
