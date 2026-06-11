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

static inline double lower_component(double diff, double degree, bool degree_is_int) {
  if (degree == 0) return diff >= 0.0 ? 1.0 : 0.0;
  if (diff < 0.0) return 0.0;
  return degree_is_int
  ? repeatMultiplication(diff, static_cast<int>(degree))
    : std::pow(diff, degree);
}

static inline double upper_component(double diff, double degree, bool degree_is_int) {
  if (degree == 0) return diff > 0.0 ? 1.0 : 0.0;
  if (diff < 0.0) return 0.0;
  return degree_is_int
  ? repeatMultiplication(diff, static_cast<int>(degree))
    : std::pow(diff, degree);
}

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
          out += repeatMultiplication(value, static_cast<int>(degree));
        }
      } else {
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
          out += repeatMultiplication(value, static_cast<int>(degree));
        }
      } else {
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


// ============================================================================
// Batched nD CoLPM backend
// ============================================================================
//
// Computes CoLPM_nD(data, target_row, degree, norm) for every row of `targets`
// in one C++ call.  This replaces the R-side pattern:
//   apply(variable, 1, function(row) Co.LPM_nD(variable, row, degree = degree))
//
// Semantics match clpm_nD_cpp():
//   degree == 0 returns the raw lower-count probability, regardless of norm.
//   degree  > 0 returns raw CLPM if norm = false.
//   degree  > 0 returns CLPM / (CLPM + CUPM + DPM) if norm = true.

struct CoLPMnDBatchWorker : public Worker {
  const RMatrix<double> data;
  const RMatrix<double> targets;
  const double degree;
  const bool norm;
  const bool degree_is_int;
  RVector<double> output;
  
  CoLPMnDBatchWorker(const NumericMatrix& data_,
                     const NumericMatrix& targets_,
                     double degree_,
                     bool norm_,
                     NumericVector& output_)
    : data(data_),
      targets(targets_),
      degree(degree_),
      norm(norm_),
      degree_is_int(isInteger(degree_)),
      output(output_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    const std::size_t n_obs = data.nrow();
    const std::size_t d = data.ncol();
    
    for (std::size_t r = begin; r < end; ++r) {
      
      // Match clpm_nD_cpp degree == 0 behavior:
      // it returns raw count probability and does not apply norm.
      if (degree == 0.0) {
        double count = 0.0;
        
        for (std::size_t i = 0; i < n_obs; ++i) {
          bool below_all = true;
          
          for (std::size_t j = 0; j < d; ++j) {
            if (data(i, j) > targets(r, j)) {
              below_all = false;
              break;
            }
          }
          
          if (below_all) count += 1.0;
        }
        
        output[r] = count / static_cast<double>(n_obs);
        continue;
      }
      
      double clpm_sum = 0.0;
      double cupm_sum = 0.0;
      double dpm_sum  = 0.0;
      
      for (std::size_t i = 0; i < n_obs; ++i) {
        double lower_prod = 1.0;
        double upper_prod = 1.0;
        double dpm_prod   = 1.0;
        
        bool all_below_strict = true;
        bool all_above_strict = true;
        
        for (std::size_t j = 0; j < d; ++j) {
          const double diff = data(i, j) - targets(r, j);
          
          // CLPM component: target - data
          lower_prod *= lower_component(-diff, degree, degree_is_int);
          
          // CUPM component: data - target
          upper_prod *= upper_component(diff, degree, degree_is_int);
          
          // Match DpmSumWorker strict all-below/all-above logic.
          if (diff >= 0.0) all_below_strict = false;
          if (diff <= 0.0) all_above_strict = false;
          
          dpm_prod *= degree_is_int
          ? repeatMultiplication(std::abs(diff), static_cast<int>(degree))
            : std::pow(std::abs(diff), degree);
        }
        
        clpm_sum += lower_prod;
        cupm_sum += upper_prod;
        
        if (!(all_below_strict || all_above_strict)) {
          dpm_sum += dpm_prod;
        }
      }
      
      const double inv_n = 1.0 / static_cast<double>(n_obs);
      const double clpm_un = clpm_sum * inv_n;
      
      if (!norm) {
        output[r] = clpm_un;
      } else {
        const double cupm_un = cupm_sum * inv_n;
        const double dpm_un  = dpm_sum  * inv_n;
        const double norm_const = clpm_un + cupm_un + dpm_un;
        
        output[r] = norm_const > 0.0 ? clpm_un / norm_const : 0.0;
      }
    }
  }
};


// [[Rcpp::export]]
NumericVector CoLPM_nD_batch_RCPP(const NumericMatrix& data,
                                  const NumericMatrix& targets,
                                  double degree = 0.0,
                                  bool norm = true) {
  if (data.ncol() != targets.ncol()) {
    stop("`targets` must have the same number of columns as `data`");
  }
  
  if (data.nrow() == 0) {
    stop("`data` must have at least one row");
  }
  
  NumericVector output(targets.nrow());
  
  CoLPMnDBatchWorker worker(data, targets, degree, norm, output);
  parallelFor(0, targets.nrow(), worker);
  
  return output;
}

// parallelFor
#define NNS_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS)      \
size_t target_size=target.size();                        \
NumericVector output = NumericVector(target_size);       \
WORKER_CLASS tmp_func(degree, target, variable, output); \
parallelFor(0, target_size, tmp_func);                   \
return(output);

// Scalar guard: the prefix backend costs O(n log n + n*degree) to build,
// which only amortizes over many targets (crossover ~ log2(n) + degree).
// For few targets, a direct O(n) scan per target via the legacy kernels is
// faster and bit-identical to pre-13.0 semantics.
static const R_xlen_t NNS_DIRECT_PATH_MAX_TARGETS = 32;

// [[Rcpp::export]]
NumericVector LPM_CPv(const double &degree,
                      const NumericVector &target,
                      const NumericVector &variable) {
  if (target.size() <= NNS_DIRECT_PATH_MAX_TARGETS) {
    NumericVector output(target.size());
    RcppParallel::RVector<double> v(variable);
    
    for (R_xlen_t i = 0; i < target.size(); ++i) {
      output[i] = LPM_C(degree, target[i], v);
    }
    
    return output;
  }
  
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(LPM_Worker);
}

// [[Rcpp::export]]
NumericVector UPM_CPv(const double &degree,
                      const NumericVector &target,
                      const NumericVector &variable) {
  if (target.size() <= NNS_DIRECT_PATH_MAX_TARGETS) {
    NumericVector output(target.size());
    RcppParallel::RVector<double> v(variable);
    
    for (R_xlen_t i = 0; i < target.size(); ++i) {
      output[i] = UPM_C(degree, target[i], v);
    }
    
    return output;
  }
  
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
    const double &degree_x, const double &degree_y,
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
  bool d_x_0=(degree_x==0);
  bool d_y_0=(degree_y==0);
  bool x_is_int=isInteger(degree_x);
  bool y_is_int=isInteger(degree_y);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    double y1=(y[i]-target_y);
    
    if(d_x_0) x1 = (x1 > 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_y_0) y1 = (y1 > 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(!d_x_0){
      if(x_is_int) x1 = repeatMultiplication(x1, static_cast<int>(degree_x));
      else x1 = std::pow(x1, degree_x);
    }
    if(!d_y_0){
      if(y_is_int) y1 = repeatMultiplication(y1, static_cast<int>(degree_y));
      else y1 = std::pow(y1, degree_y);
    }
    out += x1 * y1;
  }
  return out/max_size;
}

double CoLPM_C(
    const double &degree_x, const double &degree_y,
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
  bool d_x_0=(degree_x==0);
  bool d_y_0=(degree_y==0);
  bool x_is_int=isInteger(degree_x);
  bool y_is_int=isInteger(degree_y);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    double y1=(target_y-y[i]);
    
    if(d_x_0) x1 = (x1 >= 0 ? 1 : 0);
    else x1 = (x1 < 0 ? 0 : x1);
    
    if(d_y_0) y1 = (y1 >= 0 ? 1 : 0);
    else y1 = (y1 < 0 ? 0 : y1);
    
    if(!d_x_0){
      if(x_is_int) x1 = repeatMultiplication(x1, static_cast<int>(degree_x));
      else x1 = std::pow(x1, degree_x);
    }
    if(!d_y_0){
      if(y_is_int) y1 = repeatMultiplication(y1, static_cast<int>(degree_y));
      else y1 = std::pow(y1, degree_y);
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

// Retained for absolute backward compatibility with internal single-pair calls
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
  
  coLpm = 0.0;
  coUpm = 0.0;
  dLpm = 0.0;
  dUpm = 0.0;
  covMat=0;
  if(rows == 0)
    return;
  
  bool lpm_is_int = isInteger(degree_lpm);
  bool upm_is_int = isInteger(degree_upm);
  for(size_t i=0; i<rows; i++){
    double x_lower = lower_component(target_x - x_rvec[i], degree_lpm, lpm_is_int);
    double x_upper = upper_component(x_rvec[i] - target_x, degree_upm, upm_is_int);
    double y_lower = lower_component(target_y - y_rvec[i], degree_lpm, lpm_is_int);
    double y_upper = upper_component(y_rvec[i] - target_y, degree_upm, upm_is_int);
    
    coLpm += x_lower * y_lower;
    coUpm += x_upper * y_upper;
    dLpm += x_upper * y_lower;
    dUpm += x_lower * y_upper;
  }
  
  double inv_rows = 1.0 / static_cast<double>(rows);
  coLpm *= inv_rows;
  coUpm *= inv_rows;
  dLpm *= inv_rows;
  dUpm *= inv_rows;
  
  if(pop_adj && rows > 1 && degree_lpm > 0 && degree_upm > 0){
    coLpm *= adjust;
    coUpm *= adjust;
    dLpm *= adjust;
    dUpm *= adjust;
  }
  covMat = coUpm + coLpm - dUpm - dLpm;
}

// ============================================================================
// ULTRA-OPTIMIZED TENSORIZED MULTIVARIATE INTERNALS
// ============================================================================

// Worker 1: Compute Deviation Matrices exactly ONCE per element.
// Perfectly column-contiguous, cache-friendly SIMD streaming.
struct PrecomputeDeviationsWorker : public Worker {
  const RMatrix<double> variable;
  const RVector<double> target;
  double degree_lpm;
  double degree_upm;
  bool lpm_is_int;
  bool upm_is_int;
  
  RMatrix<double> D_lower;
  RMatrix<double> D_upper;
  
  PrecomputeDeviationsWorker(const NumericMatrix& variable_, const NumericVector& target_,
                             double degree_lpm_, double degree_upm_,
                             NumericMatrix& D_lower_, NumericMatrix& D_upper_)
    : variable(variable_), target(target_),
      degree_lpm(degree_lpm_), degree_upm(degree_upm_),
      lpm_is_int(isInteger(degree_lpm_)), upm_is_int(isInteger(degree_upm_)),
      D_lower(D_lower_), D_upper(D_upper_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    size_t rows = variable.nrow();
    for (std::size_t j = begin; j < end; ++j) {
      double t_j = target[j];
      for (size_t i = 0; i < rows; ++i) {
        double val = variable(i, j);
        D_lower(i, j) = lower_component(t_j - val, degree_lpm, lpm_is_int);
        D_upper(i, j) = upper_component(val - t_j, degree_upm, upm_is_int);
      }
    }
  }
};

// Worker 2: Blistering Fused Matrix Multiplication (t(D) %*% D)
// Completely stripped of all conditions, branching, and pow() calls.
struct FusedMatrixMultiplicationWorker : public Worker {
  const RMatrix<double> D_lower;
  const RMatrix<double> D_upper;
  bool apply_adj;
  double adjust;
  size_t rows;
  
  RMatrix<double> coLpm;
  RMatrix<double> coUpm;
  RMatrix<double> dLpm;
  RMatrix<double> dUpm;
  RMatrix<double> covMat;
  
  FusedMatrixMultiplicationWorker(const NumericMatrix& D_lower_, const NumericMatrix& D_upper_,
                                  bool apply_adj_, double adjust_, size_t rows_,
                                  NumericMatrix& coLpm_, NumericMatrix& coUpm_,
                                  NumericMatrix& dLpm_, NumericMatrix& dUpm_, NumericMatrix& covMat_)
    : D_lower(D_lower_), D_upper(D_upper_), apply_adj(apply_adj_), adjust(adjust_), rows(rows_),
      coLpm(coLpm_), coUpm(coUpm_), dLpm(dLpm_), dUpm(dUpm_), covMat(covMat_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    size_t cols = D_lower.ncol();
    double inv_rows = 1.0 / static_cast<double>(rows);
    
    for (std::size_t i = begin; i < end; ++i) {
      // PM.matrix quadrant symmetry:
      //   CUPM(i,j) = CUPM(j,i)
      //   CLPM(i,j) = CLPM(j,i)
      //   DUPM(i,j) = DLPM(j,i)
      //   DLPM(i,j) = DUPM(j,i)
      // Therefore compute only the upper triangle and cross-mirror DUPM/DLPM.
      for (std::size_t j = i; j < cols; ++j) {
        double sum_cupm = 0.0;
        double sum_clpm = 0.0;
        double sum_dupm = 0.0;
        double sum_dlpm = 0.0;
        
        // Loop fusion: Compute all 4 co-moment quadrants in a single hot-cache row scan.
        for (size_t k = 0; k < rows; ++k) {
          double u_i = D_upper(k, i);
          double l_i = D_lower(k, i);
          double u_j = D_upper(k, j);
          double l_j = D_lower(k, j);
          
          sum_cupm += u_i * u_j;
          sum_clpm += l_i * l_j;
          sum_dupm += l_i * u_j;
          sum_dlpm += u_i * l_j;
        }
        
        sum_cupm *= inv_rows;
        sum_clpm *= inv_rows;
        sum_dupm *= inv_rows;
        sum_dlpm *= inv_rows;
        
        if (apply_adj) {
          sum_cupm *= adjust;
          sum_clpm *= adjust;
          sum_dupm *= adjust;
          sum_dlpm *= adjust;
        }
        
        double cov_ij = sum_cupm + sum_clpm - sum_dupm - sum_dlpm;
        
        coUpm(i, j) = sum_cupm;
        coLpm(i, j) = sum_clpm;
        dUpm(i, j)  = sum_dupm;
        dLpm(i, j)  = sum_dlpm;
        covMat(i, j) = cov_ij;
        
        if (j != i) {
          coUpm(j, i) = sum_cupm;
          coLpm(j, i) = sum_clpm;
          
          // Crossed mirror, not ordinary symmetry.
          dUpm(j, i)  = sum_dlpm;
          dLpm(j, i)  = sum_dupm;
          covMat(j, i) = cov_ij;
        }
      }
    }
  }
};

// [[Rcpp::export]]
List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const NumericVector &target,
    const NumericMatrix &variable,
    const bool &pop_adj,
    const bool &norm
) {
  size_t variable_cols = variable.cols();
  size_t target_length = target.size();
  if(variable_cols != target_length){
    Rcpp::stop("variable matrix cols != target vector length");
    return List::create();
  }
  
  size_t rows = variable.rows();
  if (rows == 0) return List::create();
  
  // 1. Allocate continuous intermediate deviation matrices
  NumericMatrix D_lower(rows, variable_cols);
  NumericMatrix D_upper(rows, variable_cols);
  
  // 2. Step 1: Precompute all element deviation components in parallel
  PrecomputeDeviationsWorker precalc_engine(variable, target, LPM_degree, UPM_degree, D_lower, D_upper);
  parallelFor(0, variable_cols, precalc_engine);
  
  // 3. Allocate final return matrix structures
  NumericMatrix coLpm(variable_cols, variable_cols);
  NumericMatrix coUpm(variable_cols, variable_cols);
  NumericMatrix dLpm(variable_cols, variable_cols);
  NumericMatrix dUpm(variable_cols, variable_cols);
  NumericMatrix covMat(variable_cols, variable_cols);
  
  // 4. Determine population adjustment configurations
  double adjust = 1.0;
  if (pop_adj && rows > 1) {
    adjust = static_cast<double>(rows) / static_cast<double>(rows - 1);
  }
  bool apply_adj = pop_adj && rows > 1 && LPM_degree > 0 && UPM_degree > 0;
  
  // 5. Step 2: High-speed matrix contraction loops across available cores
  FusedMatrixMultiplicationWorker matrix_engine(D_lower, D_upper, apply_adj, adjust, rows,
                                                coLpm, coUpm, dLpm, dUpm, covMat);
  parallelFor(0, variable_cols, matrix_engine);
  
  // 6. Apply cellular normalization adjustments if requested.
  // Preserve the same cross-transpose relationship for DUPM and DLPM.
  if (norm) {
    for (size_t i = 0; i < variable_cols; ++i) {
      for (size_t j = i; j < variable_cols; ++j) {
        double cupm_ij = coUpm(i, j);
        double dupm_ij = dUpm(i, j);
        double dlpm_ij = dLpm(i, j);
        double clpm_ij = coLpm(i, j);
        double total = cupm_ij + dupm_ij + dlpm_ij + clpm_ij;
        
        if (total > 0.0) {
          cupm_ij /= total;
          dupm_ij /= total;
          dlpm_ij /= total;
          clpm_ij /= total;
        } else {
          cupm_ij = 0.0;
          dupm_ij = 0.0;
          dlpm_ij = 0.0;
          clpm_ij = 0.0;
        }
        
        double cov_ij = cupm_ij + clpm_ij - dupm_ij - dlpm_ij;
        
        coUpm(i, j) = cupm_ij;
        coLpm(i, j) = clpm_ij;
        dUpm(i, j)  = dupm_ij;
        dLpm(i, j)  = dlpm_ij;
        covMat(i, j) = cov_ij;
        
        if (j != i) {
          coUpm(j, i) = cupm_ij;
          coLpm(j, i) = clpm_ij;
          
          // Crossed mirror after normalization too.
          dUpm(j, i)  = dlpm_ij;
          dLpm(j, i)  = dupm_ij;
          covMat(j, i) = cov_ij;
        }
      }
    }
  }
  
  // 7. Shape attribute text allocations
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