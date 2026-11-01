// partial_moments.h
#ifndef NNS_partial_moments_H
#define NNS_partial_moments_H

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

// Backend API for the partial moment computations. These routines operate on
// RcppParallel vector and matrix proxies so they can be reused from serial and
// parallel workers. Higher-level wrappers that accept generic R objects live in
// partial_moments_rcpp.h/cpp.

/////////////////
// UPM / LPM
// single thread
double LPM_C(const double &degree,
             const double &target,
             const RcppParallel::RVector<double> &variable);
double UPM_C(const double &degree,
             const double &target,
             const RcppParallel::RVector<double> &variable);

// Rejects negative and NaN degrees. Called from the R-facing wrappers before
// any parallel work starts.
inline void nns_check_degree(const double degree, const char *name) {
  if (!(degree >= 0.0)) {
    Rcpp::stop("'%s' must be a non-negative number", name);
  }
}

namespace nns_pm_detail {

// Keep very large integer degrees on the legacy full-scan path. This prevents
// accidental allocation of many prefix-power columns while still accelerating
// the hot NNS use cases: degree 0, 1, 2, and other small integer degrees.
static const int PREFIX_MAX_DEGREE = 32;

inline bool prefix_supported_degree(const double degree, int &degree_int) {
  if (!std::isfinite(degree) || degree < 0.0) return false;
  
  const double rounded = std::round(degree);
  if (std::fabs(degree - rounded) > 1e-12) return false;
  if (rounded > static_cast<double>(PREFIX_MAX_DEGREE)) return false;
  
  degree_int = static_cast<int>(rounded);
  return true;
}

// Sorted-anchor prefix backend. Every expansion is anchored at the adjacent
// order statistic, so all terms are non-negative and cannot cancel (anchoring
// at the median instead lets terms cancel from degree 3 up):
//   lower_j(k) = sum_{i < k}  (x_(k-1) - x_(i))^j  (anchor: largest point <= t)
//   upper_j(k) = sum_{i >= k} (x_(i) - x_(k))^j    (anchor: smallest point > t)
// built with the exact recurrence
//   lower_j(k+1) = sum_{l <= j} C(j,l) delta^(j-l) lower_l(k) + [j == 0],
// and a target with k = #{x_i <= t} is evaluated as
//   LPM = sum_j C(d,j) (t - x_(k-1))^(d-j) lower_j(k) / n
//   UPM = sum_j C(d,j) (x_(k) - t)^(d-j)   upper_j(k) / n.
struct PrefixPartialMomentBackend {
  std::vector<double> sorted;
  std::vector<double> lower;   // (n + 1) x (degree + 1), row k = lower_j(k)
  std::vector<double> upper;   // (n + 1) x (degree + 1), row k = upper_j(k)
  std::vector<double> binom;  // (degree + 1) x (degree + 1), C(j, l)
  std::size_t n;
  int degree;
  
  PrefixPartialMomentBackend(const Rcpp::NumericVector &variable,
                             const int degree_)
    : sorted(variable.begin(), variable.end()),
      n(sorted.size()),
      degree(degree_) {
    
    // Match the legacy path for missing/non-finite data by declining the prefix
    // backend. The constructor is only called after this same condition is
    // checked, so this is a defensive guard.
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(sorted[i])) {
        sorted.clear();
        n = 0;
        return;
      }
    }
    
    std::sort(sorted.begin(), sorted.end());
    const std::size_t w = static_cast<std::size_t>(degree) + 1U;
    
    binom.assign(w * w, 0.0);
    for (std::size_t j = 0; j < w; ++j) {
      binom[j * w] = 1.0;
      for (std::size_t l = 1; l <= j; ++l) {
        binom[j * w + l] = binom[(j - 1U) * w + l - 1U] +
          (l < j ? binom[(j - 1U) * w + l] : 0.0);
      }
    }
    
    lower.assign((n + 1U) * w, 0.0);
    upper.assign((n + 1U) * w, 0.0);
    if (n == 0U) return;
    std::vector<double> dpow(w);
    
    lower[1U * w] = 1.0;  // k = 1: the single point equals its anchor
    for (std::size_t k = 1; k < n; ++k) {
      fill_powers(sorted[k] - sorted[k - 1U], dpow);
      step(&lower[k * w], &lower[(k + 1U) * w], dpow);
    }
    upper[(n - 1U) * w] = 1.0;  // k = n - 1: the single point equals its anchor
    for (std::size_t k = n - 1U; k > 0; --k) {
      fill_powers(sorted[k] - sorted[k - 1U], dpow);
      step(&upper[k * w], &upper[(k - 1U) * w], dpow);
    }
  }
  
  void fill_powers(const double base, std::vector<double> &pw) const {
    pw[0] = 1.0;
    for (std::size_t p = 1; p < pw.size(); ++p) pw[p] = pw[p - 1U] * base;
  }
  
  // next_j = sum_{l <= j} C(j,l) delta^(j-l) prev_l + [j == 0]
  void step(const double *prev, double *next, const std::vector<double> &dpow) const {
    const std::size_t w = static_cast<std::size_t>(degree) + 1U;
    for (std::size_t j = 0; j < w; ++j) {
      double acc = (j == 0U) ? 1.0 : 0.0;
      for (std::size_t l = 0; l <= j; ++l) {
        if (prev[l] != 0.0) acc += binom[j * w + l] * dpow[j - l] * prev[l];
      }
      next[j] = acc;
    }
  }
  
  // sum_j C(d,j) gap^(d-j) sums_j, all terms >= 0.
  double combine(const double *sums, const double gap) const {
    const std::size_t d = static_cast<std::size_t>(degree);
    const std::size_t w = d + 1U;
    double out = 0.0;
    double gpow = 1.0;  // gap^(d-j), built from j = d downwards
    for (std::size_t jj = w; jj-- > 0U;) {
      if (sums[jj] != 0.0) out += binom[d * w + jj] * gpow * sums[jj];
      gpow *= gap;
    }
    return out;
  }
  
  bool ok() const {
    return n > 0U;
  }
  
  std::size_t count_leq(const double target) const {
    return static_cast<std::size_t>(
      std::upper_bound(sorted.begin(), sorted.end(), target) - sorted.begin()
    );
  }
  
  double lpm(const double target) const {
    if (!std::isfinite(target)) return R_NaN;
    
    const std::size_t k = count_leq(target);
    const double nd = static_cast<double>(n);
    
    if (degree == 0) return static_cast<double>(k) / nd;
    if (k == 0U) return 0.0;
    
    const std::size_t w = static_cast<std::size_t>(degree) + 1U;
    return combine(&lower[k * w], target - sorted[k - 1U]) / nd;
  }
  
  double upm(const double target) const {
    if (!std::isfinite(target)) return R_NaN;
    
    const std::size_t k = count_leq(target);
    const double nd = static_cast<double>(n);
    
    if (degree == 0) return static_cast<double>(n - k) / nd;
    if (k == n) return 0.0;
    
    const std::size_t w = static_cast<std::size_t>(degree) + 1U;
    return combine(&upper[k * w], sorted[k] - target) / nd;
  }
  
  std::pair<double, double> both(const double target) const {
    return std::make_pair(lpm(target), upm(target));
  }
};

inline std::shared_ptr<const PrefixPartialMomentBackend>
  make_prefix_backend(const double degree, const Rcpp::NumericVector &variable) {
    int degree_int = 0;
    if (variable.size() == 0) {
      return std::shared_ptr<const PrefixPartialMomentBackend>();
    }
    
    if (!prefix_supported_degree(degree, degree_int)) {
      return std::shared_ptr<const PrefixPartialMomentBackend>();
    }
    
    for (R_xlen_t i = 0; i < variable.size(); ++i) {
      const double v = variable[i];
      if (!std::isfinite(v)) {
        return std::shared_ptr<const PrefixPartialMomentBackend>();
      }
    }
    
    return std::shared_ptr<const PrefixPartialMomentBackend>(
        new PrefixPartialMomentBackend(variable, degree_int)
    );
  }

} // namespace nns_pm_detail

// parallelFor
struct LPM_Worker : public RcppParallel::Worker
{
  const double degree;
  const RcppParallel::RVector<double> target;
  const RcppParallel::RVector<double> variable;
  RcppParallel::RVector<double> output;
  std::shared_ptr<const nns_pm_detail::PrefixPartialMomentBackend> prefix;
  
  LPM_Worker(
    const double degree,
    const Rcpp::NumericVector &target,
    const Rcpp::NumericVector &variable,
    Rcpp::NumericVector &output
  ):
    degree(degree), target(target), variable(variable), output(output),
    prefix(nns_pm_detail::make_prefix_backend(degree, variable)) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    if (prefix) {
      for (std::size_t i = begin; i < end; ++i) {
        const double t = target[i];
        output[i] = std::isfinite(t) ? prefix->lpm(t) : LPM_C(degree, t, variable);
      }
    } else {
      for (std::size_t i = begin; i < end; ++i) {
        output[i] = LPM_C(degree, target[i], variable);
      }
    }
  }
};

struct UPM_Worker : public RcppParallel::Worker
{
  const double degree;
  const RcppParallel::RVector<double> target;
  const RcppParallel::RVector<double> variable;
  RcppParallel::RVector<double> output;
  std::shared_ptr<const nns_pm_detail::PrefixPartialMomentBackend> prefix;
  
  UPM_Worker(
    const double degree,
    const Rcpp::NumericVector &target,
    const Rcpp::NumericVector &variable,
    Rcpp::NumericVector &output
  ):
    degree(degree), target(target), variable(variable), output(output),
    prefix(nns_pm_detail::make_prefix_backend(degree, variable)) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    if (prefix) {
      for (std::size_t i = begin; i < end; ++i) {
        const double t = target[i];
        output[i] = std::isfinite(t) ? prefix->upm(t) : UPM_C(degree, t, variable);
      }
    } else {
      for (std::size_t i = begin; i < end; ++i) {
        output[i] = UPM_C(degree, target[i], variable);
      }
    }
  }
};

// Use these workers in LPM_ratio_CPv / UPM_ratio_CPv to avoid computing the
// lower and upper partial moments through two separate vectorized kernels.
struct LPM_Ratio_Worker : public RcppParallel::Worker
{
  const double degree;
  const RcppParallel::RVector<double> target;
  const RcppParallel::RVector<double> variable;
  RcppParallel::RVector<double> output;
  std::shared_ptr<const nns_pm_detail::PrefixPartialMomentBackend> prefix;
  
  LPM_Ratio_Worker(
    const double degree,
    const Rcpp::NumericVector &target,
    const Rcpp::NumericVector &variable,
    Rcpp::NumericVector &output
  ):
    degree(degree), target(target), variable(variable), output(output),
    prefix(nns_pm_detail::make_prefix_backend(degree, variable)) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    if (prefix) {
      for (std::size_t i = begin; i < end; ++i) {
        const double t = target[i];
        if (std::isfinite(t)) {
          const std::pair<double, double> pm = prefix->both(t);
          output[i] = pm.first / (pm.first + pm.second);
        } else {
          const double lpm = LPM_C(degree, t, variable);
          const double upm = UPM_C(degree, t, variable);
          output[i] = lpm / (lpm + upm);
        }
      }
    } else {
      for (std::size_t i = begin; i < end; ++i) {
        const double t = target[i];
        const double lpm = LPM_C(degree, t, variable);
        const double upm = UPM_C(degree, t, variable);
        output[i] = lpm / (lpm + upm);
      }
    }
  }
};

struct UPM_Ratio_Worker : public RcppParallel::Worker
{
  const double degree;
  const RcppParallel::RVector<double> target;
  const RcppParallel::RVector<double> variable;
  RcppParallel::RVector<double> output;
  std::shared_ptr<const nns_pm_detail::PrefixPartialMomentBackend> prefix;
  
  UPM_Ratio_Worker(
    const double degree,
    const Rcpp::NumericVector &target,
    const Rcpp::NumericVector &variable,
    Rcpp::NumericVector &output
  ):
    degree(degree), target(target), variable(variable), output(output),
    prefix(nns_pm_detail::make_prefix_backend(degree, variable)) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    if (prefix) {
      for (std::size_t i = begin; i < end; ++i) {
        const double t = target[i];
        if (std::isfinite(t)) {
          const std::pair<double, double> pm = prefix->both(t);
          output[i] = pm.second / (pm.first + pm.second);
        } else {
          const double lpm = LPM_C(degree, t, variable);
          const double upm = UPM_C(degree, t, variable);
          output[i] = upm / (lpm + upm);
        }
      }
    } else {
      for (std::size_t i = begin; i < end; ++i) {
        const double t = target[i];
        const double lpm = LPM_C(degree, t, variable);
        const double upm = UPM_C(degree, t, variable);
        output[i] = upm / (lpm + upm);
      }
    }
  }
};

Rcpp::NumericVector LPM_CPv(const double &degree,
                            const Rcpp::NumericVector &target,
                            const Rcpp::NumericVector &variable);
Rcpp::NumericVector UPM_CPv(const double &degree,
                            const Rcpp::NumericVector &target,
                            const Rcpp::NumericVector &variable);
Rcpp::NumericVector LPM_ratio_CPv(const double &degree,
                                  const Rcpp::NumericVector &target,
                                  const Rcpp::NumericVector &variable);
Rcpp::NumericVector UPM_ratio_CPv(const double &degree,
                                  const Rcpp::NumericVector &target,
                                  const Rcpp::NumericVector &variable);

/////////////////
// CoUPM / CoLPM / DUPM / DLPM
// single thread
double CoUPM_C(
    const double &degree_x, const double &degree_y,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);
double CoLPM_C(
    const double &degree_x, const double &degree_y,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);
double DLPM_C(
    const double &degree_lpm, const double &degree_upm,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);
double DUPM_C(
    const double &degree_lpm, const double &degree_upm,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);

// parallelFor
#define NNS_PM_TWO_VARIABLES_WORKER(NAME, FUNC)                                             \
struct NAME : public RcppParallel::Worker                                                   \
{                                                                                           \
  const double degree_lpm;                                                                  \
  const double degree_upm;                                                                  \
  const RcppParallel::RVector<double> x;                                                    \
  const RcppParallel::RVector<double> y;                                                    \
  const RcppParallel::RVector<double> target_x;                                             \
  const RcppParallel::RVector<double> target_y;                                             \
  const size_t n_t_x;                                                                       \
  const size_t n_t_y;                                                                       \
  RcppParallel::RVector<double> output;                                                     \
  NAME (                                                                                    \
      const double degree_lpm,                                                              \
      const double degree_upm,                                                              \
      const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,                           \
      const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y,             \
      Rcpp::NumericVector &output                                                           \
  ):                                                                                        \
    degree_lpm(degree_lpm), degree_upm(degree_upm),                                         \
    x(x), y(y), target_x(target_x), target_y(target_y),                                     \
    n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)                          \
  {}                                                                                        \
  void operator()(std::size_t begin, std::size_t end) {                                     \
    for (size_t i = begin; i < end; i++) {                                                  \
      output[i] = FUNC(degree_lpm, degree_upm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); \
    }                                                                                       \
  }                                                                                         \
}

NNS_PM_TWO_VARIABLES_WORKER(CoLPM_Worker, CoLPM_C);
NNS_PM_TWO_VARIABLES_WORKER(CoUPM_Worker, CoUPM_C);
NNS_PM_TWO_VARIABLES_WORKER(DLPM_Worker, DLPM_C);
NNS_PM_TWO_VARIABLES_WORKER(DUPM_Worker, DUPM_C);
Rcpp::NumericVector CoLPM_CPv(
    const double &degree_x, const double &degree_y,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);
Rcpp::NumericVector CoUPM_CPv(
    const double &degree_x, const double &degree_y,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);
Rcpp::NumericVector DLPM_CPv(
    const double &degree_lpm, const double &degree_upm,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);
Rcpp::NumericVector DUPM_CPv(
    const double &degree_lpm, const double &degree_upm,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);

/////////////////
// PM MATRIX
// single thread
void PMMatrix_Cv(
    const double &degree_lpm,
    const double &degree_upm,
    const RcppParallel::RMatrix<double>::Column &x,
    const RcppParallel::RMatrix<double>::Column &y,
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
);
// parallelFor
struct PMMatrix_Worker : public RcppParallel::Worker
{
  const double degree_lpm;
  const double degree_upm;
  const RcppParallel::RMatrix<double> variable;
  const RcppParallel::RVector<double> target;
  const size_t variable_cols;
  const size_t variable_rows;
  const size_t target_length;
  const bool pop_adj;
  double adjust;
  RcppParallel::RMatrix<double> coLpm;
  RcppParallel::RMatrix<double> coUpm;
  RcppParallel::RMatrix<double> dLpm;
  RcppParallel::RMatrix<double> dUpm;
  RcppParallel::RMatrix<double> covMat;
  PMMatrix_Worker(
    const double &degree_lpm, const double &degree_upm,
    const Rcpp::NumericMatrix &variable,
    const Rcpp::NumericVector &target,
    const bool &pop_adj,
    Rcpp::NumericMatrix &coLpm, Rcpp::NumericMatrix &coUpm,
    Rcpp::NumericMatrix &dLpm,  Rcpp::NumericMatrix &dUpm,
    Rcpp::NumericMatrix &covMat
  ):
    degree_lpm(degree_lpm), degree_upm(degree_upm),
    variable(variable), target(target),
    variable_cols(variable.cols()), variable_rows(variable.rows()), target_length(target.size()),
    pop_adj(pop_adj),
    coLpm(coLpm), coUpm(coUpm),
    dLpm(dLpm), dUpm(dUpm),
    covMat(covMat)
  {
    if(variable_cols != target_length)
      Rcpp::stop("variable matrix cols != target vector length");
    adjust = 1;
    if (variable_rows > 1)
      adjust=((double)variable_rows)/((double)variable_rows-1);
  }
  void operator()(std::size_t begin, std::size_t end) {
    for (size_t i = begin; i < end; i++){
      for (size_t l = 0; l < variable_cols; l++){
        PMMatrix_Cv(
          degree_lpm,
          degree_upm,
          variable.column(i),
          variable.column(l),
          target[i],
                target[l],
                      pop_adj,
                      adjust,
                      variable_rows,
                      coLpm(i,l),
                      coUpm(i,l),
                      dLpm(i,l),
                      dUpm(i,l),
                      covMat(i,l)
        );
      }
    }
  }
};
Rcpp::List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const Rcpp::NumericVector &target,
    const Rcpp::NumericMatrix &variable,
    const bool &pop_adj,
    const bool &norm
);

// n-D co-partial-moments prototypes (parallel back-ends)
double clpm_nD_cpp(const Rcpp::NumericMatrix &data,
                   const Rcpp::NumericVector &target,
                   double degree,
                   bool norm);

double cupm_nD_cpp(const Rcpp::NumericMatrix &data,
                   const Rcpp::NumericVector &target,
                   double degree,
                   bool norm);

double dpm_nD_cpp(const Rcpp::NumericMatrix &data,
                  const Rcpp::NumericVector &target,
                  double degree,
                  bool norm);

#endif  //NNS_partial_moments_H
