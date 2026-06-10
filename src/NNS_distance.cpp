// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// simple sample sd/var helpers
static inline double mean_vec(const std::vector<double>& v){
  if (v.empty()) return NA_REAL;
  double s = 0.0;
  for(double x : v) s += x;
  return s / (double)v.size();
}

static inline double sd_vec(const std::vector<double>& v){
  size_t n = v.size(); if (n < 2) return NA_REAL;
  double mu = mean_vec(v), acc = 0.0;
  for(double x : v){ double d = x - mu; acc += d * d; }
  return std::sqrt(acc / (double)(n - 1));
}

static inline double var_vec(const std::vector<double>& v){
  double s = sd_vec(v); return std::isfinite(s) ? s * s : NA_REAL;
}

// OPTIMIZED: Replaced std::unordered_map with a sort-and-count vector strategy.
// Bypasses OS heap-allocation locks during multi-threaded parallel execution.
static double mode_class_weighted(const std::vector<double>& y, const std::vector<double>& w) {
  int n = y.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return y[0];
  
  std::vector<std::pair<double, long long>> items;
  items.reserve(n);
  for (int i = 0; i < n; ++i) {
    long long c = (long long)std::ceil(100.0 * w[i]);
    if (c > 0) items.push_back({y[i], c});
  }
  if (items.empty()) return NA_REAL;
  
  std::sort(items.begin(), items.end(), [](const std::pair<double, long long>& a, const std::pair<double, long long>& b) {
    return a.first < b.first;
  });
  
  double best_val = items[0].first;
  long long best_cnt = items[0].second;
  double cur_val = items[0].first;
  long long cur_cnt = items[0].second;
  
  for (size_t i = 1; i < items.size(); ++i) {
    if (items[i].first == cur_val) {
      cur_cnt += items[i].second;
    } else {
      if (cur_cnt > best_cnt) { best_cnt = cur_cnt; best_val = cur_val; }
      cur_val = items[i].first;
      cur_cnt = items[i].second;
    }
  }
  if (cur_cnt > best_cnt) { best_val = cur_val; }
  return best_val;
}

// [[Rcpp::export]]
SEXP NNS_distance_cpp(NumericMatrix X,
                      NumericVector yhat,
                      NumericVector dest,
                      int k,
                      bool use_class) {
  const int l = X.nrow();
  const int n = X.ncol();
  if (yhat.size() != l) stop("yhat length must equal nrow(X)");
  if (dest.size() != n) stop("dist.estimate length must equal ncol(X)");
  
  // OPTIMIZED: Removed in-place data mutation. Computes scales dynamically to protect original matrix.
  std::vector<double> invR(n, 0.0);
  for (int j = 0; j < n; ++j) {
    double cmin = dest[j], cmax = dest[j];
    for (int i = 0; i < l; ++i) {
      double v = X(i,j);
      if (std::isfinite(v)) { if (v < cmin) cmin = v; if (v > cmax) cmax = v; }
    }
    double range = cmax - cmin;
    if (std::isfinite(range) && range > 0.0) invR[j] = 1.0 / range;
  }
  
  std::vector<double> S(l, 0.0);
  for (int i = 0; i < l; ++i) {
    double acc = 0.0; // Demoted from long double to enable SIMD Vectorization
    for (int j = 0; j < n; ++j) {
      double a = X(i,j), b = dest[j];
      if (std::isfinite(a) && std::isfinite(b) && invR[j] > 0.0) {
        double diff = (a - b) * invR[j];
        acc += diff * diff + std::fabs(diff);
      }
    }
    S[i] = (acc == 0.0 ? 1e-10 : acc);
  }
  
  int ll = std::min(k, l);
  std::vector<int> idx(l);
  std::iota(idx.begin(), idx.end(), 0);
  auto cmp = [&](int a, int b){ return S[a] < S[b]; };
  if (ll < l) std::partial_sort(idx.begin(), idx.begin()+ll, idx.end(), cmp);
  else std::sort(idx.begin(), idx.end(), cmp);
  
  idx.resize(ll);
  std::vector<double> Ssel(ll), ysel(ll);
  for (int t = 0; t < ll; ++t) {
    int i = idx[t];
    Ssel[t] = S[i];
    ysel[t] = yhat[i];
  }
  
  if (ll == 1) return wrap(ysel[0]);
  if (k == 1) {
    double smin = *std::min_element(Ssel.begin(), Ssel.end());
    std::vector<double> yties;
    for (int t = 0; t < ll; ++t) if (Ssel[t] == smin) yties.push_back(ysel[t]);
    if (yties.size() == 1) return wrap(yties[0]);
    std::vector<double> fake_w(yties.size(), 1.0);
    return wrap(mode_class_weighted(yties, fake_w));
  }
  
  std::vector<double> uni(ll, 1.0 / (double)ll);
  
  std::vector<double> tw(ll, 0.0);
  for (int i = 0; i < ll; ++i) {
    double dens = ::Rf_dt(Ssel[i], (double)ll, 0);
    tw[i] = std::isfinite(dens) ? dens : 0.0;
  }
  double twsum = std::accumulate(tw.begin(), tw.end(), 0.0);
  if (twsum > 0) for (double &v: tw) v /= twsum; else std::fill(tw.begin(), tw.end(), 0.0);
  
  std::vector<double> emp(ll, 0.0);
  for (int i = 0; i < ll; ++i){ double v = Ssel[i]; emp[i] = (v>0) ? 1.0/v : 0.0; }
  double empsum = std::accumulate(emp.begin(), emp.end(), 0.0);
  if (empsum > 0) for (double &v: emp) v /= empsum; else std::fill(emp.begin(), emp.end(), 0.0);
  
  std::vector<double> exw(ll, 0.0);
  for (int i = 0; i < ll; ++i){
    double dens = ::Rf_dexp((double)(i+1), 1.0/(double)ll, 0);
    exw[i] = std::isfinite(dens) ? dens : 0.0;
  }
  double exsum = std::accumulate(exw.begin(), exw.end(), 0.0);
  if (exsum > 0) for (double &v: exw) v /= exsum; else std::fill(exw.begin(), exw.end(), 0.0);
  
  std::vector<double> lnorm(ll, 0.0);
  double sd_ranks = NA_REAL;
  if (ll >= 2){
    std::vector<double> ranks(ll); for(int i=0; i<ll; ++i) ranks[i] = (double)(i+1);
    sd_ranks = sd_vec(ranks);
  }
  if (std::isfinite(sd_ranks)){
    for (int i = 0; i < ll; ++i){
      double lp = ::Rf_dlnorm((double)(i+1), 0.0, sd_ranks, 1);
      lnorm[i] = std::fabs(lp);
    }
    std::reverse(lnorm.begin(), lnorm.end());
  } else {
    std::fill(lnorm.begin(), lnorm.end(), 0.0);
  }
  double lnsum = std::accumulate(lnorm.begin(), lnorm.end(), 0.0);
  if (lnsum > 0) for (double &v: lnorm) v /= lnsum; else std::fill(lnorm.begin(), lnorm.end(), 0.0);
  
  std::vector<double> pl(ll, 0.0);
  for (int i = 0; i < ll; ++i){ double r = (double)(i+1); pl[i] = std::pow(r, -2.0); }
  double plsum = std::accumulate(pl.begin(), pl.end(), 0.0);
  if (plsum > 0) for (double &v: pl) v /= plsum; else std::fill(pl.begin(), pl.end(), 0.0);
  
  std::vector<double> normw(ll, 0.0);
  double sdS = sd_vec(Ssel);
  if (std::isfinite(sdS) && sdS > 0){
    for (int i = 0; i < ll; ++i){
      double dens = ::Rf_dnorm4(Ssel[i], 0.0, sdS, 0);
      normw[i] = std::isfinite(dens) ? dens : 0.0;
    }
    double nsum = std::accumulate(normw.begin(), normw.end(), 0.0);
    if (nsum > 0) for (double &v: normw) v /= nsum; else std::fill(normw.begin(), normw.end(), 0.0);
  }
  
  std::vector<double> rbf(ll, 0.0);
  double varS = var_vec(Ssel);
  if (std::isfinite(varS) && varS > 0){
    for (int i = 0; i < ll; ++i) rbf[i] = std::exp(- Ssel[i] / (2.0*varS));
    double rsum = std::accumulate(rbf.begin(), rbf.end(), 0.0);
    if (rsum > 0) for (double &v: rbf) v /= rsum; else std::fill(rbf.begin(), rbf.end(), 0.0);
  }
  
  std::vector<double> w(ll, 0.0);
  double tot = 0.0;
  for (int i = 0; i < ll; ++i){
    double wi = uni[i] + tw[i] + emp[i] + exw[i] + lnorm[i] + pl[i] + normw[i] + rbf[i];
    w[i] = wi; tot += wi;
  }
  if (tot > 0) for (double &v: w) v /= tot; else for (double &v: w) v = 1.0/(double)ll;
  
  if (!use_class){
    double dot = 0.0;
    for (int i = 0; i < ll; ++i) dot += ysel[i] * w[i];
    return wrap(dot);
  } else {
    return wrap( mode_class_weighted(ysel, w) );
  }
}

// ---------- NNS_distance_path_cpp / NNS_distance_bulk_cpp ----------
namespace {
inline double safe_eps() { return 1e-12; }
  
  inline void compute_distances(const double* rpm, int n, int p,
                                const double* test_row,
                                std::vector<double>& dist_out) {
    for (int i = 0; i < n; ++i) {
      const double* xi = rpm + static_cast<std::size_t>(i) * p;
      double acc = 0.0;
      for (int j = 0; j < p; ++j) {
        const double d = xi[j] - test_row[j];
        acc += d * d + std::fabs(d);
      }
      dist_out[i] = (acc == 0.0 ? safe_eps() : acc);
    }
  }
  
  inline void argsort_by_distance(const std::vector<double>& dist,
                                  std::vector<int>& idx) {
    const int n = static_cast<int>(dist.size());
    idx.resize(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&dist](int a, int b){ return dist[a] < dist[b]; });
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix NNS_distance_path_cpp(const Rcpp::NumericMatrix& RPM,
                                          const Rcpp::NumericVector& yhat,
                                          const Rcpp::NumericMatrix& Xtest,
                                          int kmax,
                                          bool is_class) {
  (void)is_class;
  const int n = RPM.nrow(), p = RPM.ncol(), m = Xtest.nrow();
  if (n <= 0 || p <= 0 || m <= 0) Rcpp::stop("RPM/Xtest must be non-empty");
  if (yhat.size() != n)          Rcpp::stop("yhat length must equal nrow(RPM)");
  if (Xtest.ncol() != p)         Rcpp::stop("Xtest and RPM must have same number of columns");
  if (kmax < 1)                  Rcpp::stop("kmax must be >= 1");
  if (kmax > n)                  kmax = n;
  
  Rcpp::NumericMatrix out(m, kmax);
  const double* rpm_ptr  = REAL(RPM);
  const double* y_ptr    = REAL(yhat);
  const double* tst_ptr  = REAL(Xtest);
  
  std::vector<double> dist(n), y_sorted(n), d_sorted(n);
  std::vector<int>    ord(n);
  
  for (int r = 0; r < m; ++r) {
    const double* tr = tst_ptr + static_cast<std::size_t>(r) * p;
    compute_distances(rpm_ptr, n, p, tr, dist);
    argsort_by_distance(dist, ord);
    
    for (int i = 0; i < n; ++i) {
      const int j = ord[i];
      y_sorted[i] = y_ptr[j];
      d_sorted[i] = (dist[j] <= 0.0 ? safe_eps() : dist[j]);
    }
    
    double csum_w  = 0.0, csum_yw = 0.0;
    for (int k = 1; k <= kmax; ++k) {
      const double w = 1.0 / d_sorted[k - 1];
      csum_w  += w;
      csum_yw += w * y_sorted[k - 1];
      out(r, k - 1) = (csum_w > 0.0) ? (csum_yw / csum_w) : 0.0;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector NNS_distance_bulk_cpp(const Rcpp::NumericMatrix& RPM,
                                          const Rcpp::NumericVector& yhat,
                                          const Rcpp::NumericMatrix& Xtest,
                                          int k,
                                          bool is_class) {
  (void)is_class;
  const int n = RPM.nrow(), p = RPM.ncol(), m = Xtest.nrow();
  if (n <= 0 || p <= 0 || m <= 0) Rcpp::stop("RPM/Xtest must be non-empty");
  if (yhat.size() != n)          Rcpp::stop("yhat length must equal nrow(RPM)");
  if (Xtest.ncol() != p)         Rcpp::stop("Xtest and RPM must have same number of columns");
  if (k < 1)                     Rcpp::stop("k must be >= 1");
  if (k > n)                     k = n;
  
  Rcpp::NumericVector out(m);
  const double* rpm_ptr  = REAL(RPM);
  const double* y_ptr    = REAL(yhat);
  const double* tst_ptr  = REAL(Xtest);
  
  std::vector<double> dist(n);
  std::vector<int>    ord(n);
  
  for (int r = 0; r < m; ++r) {
    const double* tr = tst_ptr + static_cast<std::size_t>(r) * p;
    compute_distances(rpm_ptr, n, p, tr, dist);
    argsort_by_distance(dist, ord);
    
    double csum_w = 0.0, csum_yw = 0.0;
    for (int i = 0; i < k; ++i) {
      const int j = ord[i];
      const double dj = (dist[j] <= 0.0 ? safe_eps() : dist[j]);
      const double w  = 1.0 / dj;
      csum_w  += w;
      csum_yw += w * y_ptr[j];
    }
    out[r] = (csum_w > 0.0) ? (csum_yw / csum_w) : 0.0;
  }
  return out;
}

// ---------- worker ----------
struct AllKWorker : public Worker {
  RMatrix<double> RPM;
  RVector<double> yhat;
  RMatrix<double> Xtest;
  std::vector<double> minRPM, maxRPM;
  int l, n, m, kmax;
  bool is_class;
  std::vector< std::vector<double> > uniW, expW, lnormW, plW;
  RMatrix<double> out;
  
  AllKWorker(NumericMatrix RPM_, NumericVector yhat_, NumericMatrix Xtest_,
             const std::vector<double>& minRPM_, const std::vector<double>& maxRPM_,
             int kmax_, bool is_class_,
             const std::vector<std::vector<double>>& uniW_,
             const std::vector<std::vector<double>>& expW_,
             const std::vector<std::vector<double>>& lnormW_,
             const std::vector<std::vector<double>>& plW_,
             NumericMatrix out_)
    : RPM(RPM_), yhat(yhat_), Xtest(Xtest_), minRPM(minRPM_), maxRPM(maxRPM_),
      l(RPM_.nrow()), n(RPM_.ncol()), m(Xtest_.nrow()), kmax(kmax_), is_class(is_class_),
      uniW(uniW_), expW(expW_), lnormW(lnormW_), plW(plW_), out(out_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> invR(n), S(l), topS, topY;
    std::vector<int> idx(l);
    
    for (std::size_t r = begin; r < end; ++r) {
      for (int j = 0; j < n; ++j){
        double t = Xtest(r,j);
        double mn = std::min(minRPM[j], t);
        double mx = std::max(maxRPM[j], t);
        double range = mx - mn;
        invR[j] = (std::isfinite(range) && range > 0.0) ? (1.0 / range) : 0.0;
      }
      
      for (int i = 0; i < l; ++i){
        double acc = 0.0; // Demoted to double for SIMD compatibility
        for (int j = 0; j < n; ++j){
          double a = RPM(i,j), b = Xtest(r,j);
          if (std::isfinite(a) && std::isfinite(b) && invR[j] > 0.0){
            double diff = (a - b) * invR[j];
            acc += diff * diff + std::fabs(diff);
          }
        }
        S[i] = (acc == 0.0 ? 1e-10 : acc);
      }
      
      std::iota(idx.begin(), idx.end(), 0);
      auto cmp = [&](int a, int b){ return S[a] < S[b]; };
      if (kmax < l) std::partial_sort(idx.begin(), idx.begin()+kmax, idx.end(), cmp);
      else          std::sort(idx.begin(), idx.end(), cmp);
      
      auto cmp2 = [&](int a, int b){
        if (S[a] < S[b]) return true;
        if (S[b] < S[a]) return false;
        return a < b;
      };
      std::stable_sort(idx.begin(), idx.begin()+kmax, cmp2);
      
      topS.resize(kmax); topY.resize(kmax);
      for (int t = 0; t < kmax; ++t){ int i = idx[t]; topS[t] = S[i]; topY[t] = yhat[i]; }
      
      for (int k = 1; k <= kmax; ++k){
        const double* Ssel = topS.data();
        const double* Ysel = topY.data();
        if (k == 1){ out(r, k-1) = Ysel[0]; continue; }
        
        std::vector<double> tw(k,0.0), emp(k,0.0), normw(k,0.0), rbf(k,0.0);
        
        // OPTIMIZED: Pure C++ Proportional Densities (Thread-Safe)
        for (int i = 0; i < k; ++i){
          // Proportional Student's T: bypasses R C-API ::Rf_dt
          tw[i] = std::pow(1.0 + (Ssel[i] * Ssel[i]) / (double)k, -(double)(k + 1) / 2.0);
          emp[i] = (Ssel[i] > 0) ? 1.0 / Ssel[i] : 0.0;
        }
        double tws = std::accumulate(tw.begin(), tw.end(), 0.0);
        if (tws > 0) for(double &v: tw) v /= tws; else std::fill(tw.begin(), tw.end(), 0.0);
        
        double emps = std::accumulate(emp.begin(), emp.end(), 0.0);
        if (emps > 0) for(double &v: emp) v /= emps; else std::fill(emp.begin(), emp.end(), 0.0);
        
        double sdS = sd_vec(std::vector<double>(topS.begin(), topS.begin() + k));
        if (std::isfinite(sdS) && sdS > 0){
          for (int i = 0; i < k; ++i){
            // Proportional Normal: bypasses R C-API ::Rf_dnorm4
            double z = Ssel[i] / sdS;
            normw[i] = std::exp(-0.5 * z * z);
          }
          double ns = std::accumulate(normw.begin(), normw.end(), 0.0);
          if (ns > 0) for(double &v: normw) v /= ns; else std::fill(normw.begin(), normw.end(), 0.0);
        }
        
        double vS = var_vec(std::vector<double>(topS.begin(), topS.begin() + k));
        if (std::isfinite(vS) && vS > 0){
          for (int i = 0; i < k; ++i) rbf[i] = std::exp(-Ssel[i] / (2.0 * vS));
          double rs = std::accumulate(rbf.begin(), rbf.end(), 0.0);
          if (rs > 0) for(double &v: rbf) v /= rs; else std::fill(rbf.begin(), rbf.end(), 0.0);
        }
        
        double dot = 0.0, tot = 0.0;
        for (int i = 0; i < k; ++i){
          double wi = uniW[k][i] + expW[k][i] + lnormW[k][i] + plW[k][i]
          + tw[i] + emp[i] + normw[i] + rbf[i];
          tot += wi;
          if (!is_class) dot += Ysel[i] * wi;
        }
        double invTot = (tot > 0.0) ? (1.0 / tot) : (1.0 / (double)k);
        
        if (!is_class){
          out(r, k-1) = (tot > 0.0) ? (dot * invTot) : (
            std::accumulate(topY.begin(), topY.begin()+k, 0.0) / (double)k
          );
        } else {
          std::vector<double> w(k);
          if (tot > 0.0) for (int i = 0; i < k; ++i) w[i] = (uniW[k][i]+expW[k][i]+lnormW[k][i]+plW[k][i]+tw[i]+emp[i]+normw[i]+rbf[i]) * invTot;
          else          std::fill(w.begin(), w.end(), 1.0/(double)k);
          out(r, k-1) = mode_class_weighted(std::vector<double>(topY.begin(), topY.begin()+k), w);
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix NNS_distance_path_parallel_cpp(NumericMatrix RPM,
                                             NumericVector yhat,
                                             NumericMatrix Xtest,
                                             int kmax,
                                             bool is_class,
                                             int nthreads = -1) {
  const int l = RPM.nrow(), n = RPM.ncol(), m = Xtest.nrow();
  if (yhat.size() != l) stop("yhat length must equal nrow(RPM)");
  if (kmax <= 0) kmax = l;
  if (kmax > l) kmax = l;
  
  std::vector<double> minRPM(n, R_PosInf), maxRPM(n, R_NegInf);
  for (int j = 0; j < n; ++j){
    for (int i = 0; i < l; ++i){
      double v = RPM(i,j);
      if (std::isfinite(v)) { if(v < minRPM[j]) minRPM[j] = v; if(v > maxRPM[j]) maxRPM[j] = v; }
    }
    if (!std::isfinite(minRPM[j])) { minRPM[j] = 0.0; maxRPM[j] = 0.0; }
  }
  
  std::vector<std::vector<double>> uniW(kmax+1), expW(kmax+1), lnormW(kmax+1), plW(kmax+1);
  for (int k = 1; k <= kmax; ++k){
    uniW[k].assign(k, 1.0 / (double)k);
    
    std::vector<double> ex(k);
    for (int r = 1; r <= k; ++r) ex[r-1] = ::Rf_dexp((double)r, 1.0 / (double)k, 0);
    double exs = std::accumulate(ex.begin(), ex.end(), 0.0);
    if (exs > 0) for (double &v: ex) v /= exs; else std::fill(ex.begin(), ex.end(), 0.0);
    expW[k] = std::move(ex);
    
    std::vector<double> pl(k);
    for (int r = 1; r <= k; ++r) pl[r-1] = std::pow((double)r, -2.0);
    double pls = std::accumulate(pl.begin(), pl.end(), 0.0);
    if (pls > 0) for (double &v: pl) v /= pls; else std::fill(pl.begin(), pl.end(), 0.0);
    plW[k] = std::move(pl);
    
    std::vector<double> ln(k, 0.0);
    if (k >= 2){
      double sdlog = std::sqrt(((double)k * (double)k - 1.0) / 12.0);
      for (int r = 1; r <= k; ++r){ double lp = ::Rf_dlnorm((double)r, 0.0, sdlog, 1); ln[r-1] = std::fabs(lp); }
      std::reverse(ln.begin(), ln.end());
      double lns = std::accumulate(ln.begin(), ln.end(), 0.0);
      if (lns > 0) for (double &v: ln) v /= lns; else std::fill(ln.begin(), ln.end(), 0.0);
    }
    lnormW[k] = std::move(ln);
  }
  
  NumericMatrix out(m, kmax);
  (void)nthreads;
  
  AllKWorker w(RPM, yhat, Xtest, minRPM, maxRPM, kmax, is_class,
               uniW, expW, lnormW, plW, out);
  
  RcppParallel::parallelFor(0, m, w);
  
  return out;
}


// ---------- single-k path ensemble worker ----------
// Computes exactly the kept column produced by NNS_distance_path_parallel_cpp(..., kmax = k)[, k]
// without evaluating the discarded path for 1:(k - 1).
struct SingleKWorker : public Worker {
  RMatrix<double> RPM;
  RVector<double> yhat;
  RMatrix<double> Xtest;
  std::vector<double> minRPM, maxRPM;
  int l, n, m, k;
  bool is_class;
  std::vector<double> uniW, expW, lnormW, plW;
  RVector<double> out;
  
  SingleKWorker(NumericMatrix RPM_, NumericVector yhat_, NumericMatrix Xtest_,
                const std::vector<double>& minRPM_, const std::vector<double>& maxRPM_,
                int k_, bool is_class_,
                const std::vector<double>& uniW_,
                const std::vector<double>& expW_,
                const std::vector<double>& lnormW_,
                const std::vector<double>& plW_,
                NumericVector out_)
    : RPM(RPM_), yhat(yhat_), Xtest(Xtest_), minRPM(minRPM_), maxRPM(maxRPM_),
      l(RPM_.nrow()), n(RPM_.ncol()), m(Xtest_.nrow()), k(k_), is_class(is_class_),
      uniW(uniW_), expW(expW_), lnormW(lnormW_), plW(plW_), out(out_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> invR(n), S(l), topS(k), topY(k);
    std::vector<int> idx(l);
    
    for (std::size_t r = begin; r < end; ++r) {
      for (int j = 0; j < n; ++j) {
        double t = Xtest(r, j);
        double mn = std::min(minRPM[j], t);
        double mx = std::max(maxRPM[j], t);
        double range = mx - mn;
        invR[j] = (std::isfinite(range) && range > 0.0) ? (1.0 / range) : 0.0;
      }
      
      for (int i = 0; i < l; ++i) {
        double acc = 0.0;
        for (int j = 0; j < n; ++j) {
          double a = RPM(i, j), b = Xtest(r, j);
          if (std::isfinite(a) && std::isfinite(b) && invR[j] > 0.0) {
            double diff = (a - b) * invR[j];
            acc += diff * diff + std::fabs(diff);
          }
        }
        S[i] = (acc == 0.0 ? 1e-10 : acc);
      }
      
      std::iota(idx.begin(), idx.end(), 0);
      auto cmp = [&](int a, int b) { return S[a] < S[b]; };
      if (k < l) std::partial_sort(idx.begin(), idx.begin() + k, idx.end(), cmp);
      else       std::sort(idx.begin(), idx.end(), cmp);
      
      auto cmp2 = [&](int a, int b) {
        if (S[a] < S[b]) return true;
        if (S[b] < S[a]) return false;
        return a < b;
      };
      std::stable_sort(idx.begin(), idx.begin() + k, cmp2);
      
      for (int t = 0; t < k; ++t) {
        int i = idx[t];
        topS[t] = S[i];
        topY[t] = yhat[i];
      }
      
      if (k == 1) {
        out[r] = topY[0];
        continue;
      }
      
      std::vector<double> tw(k, 0.0), emp(k, 0.0), normw(k, 0.0), rbf(k, 0.0);
      
      for (int i = 0; i < k; ++i) {
        tw[i] = std::pow(1.0 + (topS[i] * topS[i]) / (double)k,
                         -(double)(k + 1) / 2.0);
        emp[i] = (topS[i] > 0.0) ? 1.0 / topS[i] : 0.0;
      }
      
      double tws = std::accumulate(tw.begin(), tw.end(), 0.0);
      if (tws > 0.0) for (double &v : tw) v /= tws;
      else std::fill(tw.begin(), tw.end(), 0.0);
      
      double emps = std::accumulate(emp.begin(), emp.end(), 0.0);
      if (emps > 0.0) for (double &v : emp) v /= emps;
      else std::fill(emp.begin(), emp.end(), 0.0);
      
      double sdS = sd_vec(topS);
      if (std::isfinite(sdS) && sdS > 0.0) {
        for (int i = 0; i < k; ++i) {
          double z = topS[i] / sdS;
          normw[i] = std::exp(-0.5 * z * z);
        }
        double ns = std::accumulate(normw.begin(), normw.end(), 0.0);
        if (ns > 0.0) for (double &v : normw) v /= ns;
        else std::fill(normw.begin(), normw.end(), 0.0);
      }
      
      double vS = var_vec(topS);
      if (std::isfinite(vS) && vS > 0.0) {
        for (int i = 0; i < k; ++i) rbf[i] = std::exp(-topS[i] / (2.0 * vS));
        double rs = std::accumulate(rbf.begin(), rbf.end(), 0.0);
        if (rs > 0.0) for (double &v : rbf) v /= rs;
        else std::fill(rbf.begin(), rbf.end(), 0.0);
      }
      
      double dot = 0.0, tot = 0.0;
      for (int i = 0; i < k; ++i) {
        double wi = uniW[i] + expW[i] + lnormW[i] + plW[i]
        + tw[i] + emp[i] + normw[i] + rbf[i];
        tot += wi;
        if (!is_class) dot += topY[i] * wi;
      }
      
      double invTot = (tot > 0.0) ? (1.0 / tot) : (1.0 / (double)k);
      
      if (!is_class) {
        if (tot > 0.0) {
          out[r] = dot * invTot;
        } else {
          out[r] = std::accumulate(topY.begin(), topY.end(), 0.0) / (double)k;
        }
      } else {
        std::vector<double> w(k);
        if (tot > 0.0) {
          for (int i = 0; i < k; ++i) {
            w[i] = (uniW[i] + expW[i] + lnormW[i] + plW[i]
                      + tw[i] + emp[i] + normw[i] + rbf[i]) * invTot;
          }
        } else {
          std::fill(w.begin(), w.end(), 1.0 / (double)k);
        }
        out[r] = mode_class_weighted(topY, w);
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector NNS_distance_path_single_parallel_cpp(NumericMatrix RPM,
                                                    NumericVector yhat,
                                                    NumericMatrix Xtest,
                                                    int k,
                                                    bool is_class,
                                                    int nthreads = -1) {
  const int l = RPM.nrow(), n = RPM.ncol(), m = Xtest.nrow();
  if (yhat.size() != l) stop("yhat length must equal nrow(RPM)");
  if (Xtest.ncol() != n) stop("Xtest and RPM must have same number of columns");
  if (k <= 0) k = l;
  if (k > l) k = l;
  
  std::vector<double> minRPM(n, R_PosInf), maxRPM(n, R_NegInf);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < l; ++i) {
      double v = RPM(i, j);
      if (std::isfinite(v)) {
        if (v < minRPM[j]) minRPM[j] = v;
        if (v > maxRPM[j]) maxRPM[j] = v;
      }
    }
    if (!std::isfinite(minRPM[j])) {
      minRPM[j] = 0.0;
      maxRPM[j] = 0.0;
    }
  }
  
  std::vector<double> uniW(k, 1.0 / (double)k);
  
  std::vector<double> expW(k);
  for (int r = 1; r <= k; ++r) expW[r - 1] = ::Rf_dexp((double)r, 1.0 / (double)k, 0);
  double exs = std::accumulate(expW.begin(), expW.end(), 0.0);
  if (exs > 0.0) for (double &v : expW) v /= exs;
  else std::fill(expW.begin(), expW.end(), 0.0);
  
  std::vector<double> plW(k);
  for (int r = 1; r <= k; ++r) plW[r - 1] = std::pow((double)r, -2.0);
  double pls = std::accumulate(plW.begin(), plW.end(), 0.0);
  if (pls > 0.0) for (double &v : plW) v /= pls;
  else std::fill(plW.begin(), plW.end(), 0.0);
  
  std::vector<double> lnormW(k, 0.0);
  if (k >= 2) {
    double sdlog = std::sqrt(((double)k * (double)k - 1.0) / 12.0);
    for (int r = 1; r <= k; ++r) {
      double lp = ::Rf_dlnorm((double)r, 0.0, sdlog, 1);
      lnormW[r - 1] = std::fabs(lp);
    }
    std::reverse(lnormW.begin(), lnormW.end());
    double lns = std::accumulate(lnormW.begin(), lnormW.end(), 0.0);
    if (lns > 0.0) for (double &v : lnormW) v /= lns;
    else std::fill(lnormW.begin(), lnormW.end(), 0.0);
  }
  
  NumericVector out(m);
  (void)nthreads;
  
  SingleKWorker w(RPM, yhat, Xtest, minRPM, maxRPM, k, is_class,
                  uniW, expW, lnormW, plW, out);
  RcppParallel::parallelFor(0, m, w);
  
  return out;
}

