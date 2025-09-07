#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
using namespace Rcpp;
using namespace RcppParallel;

// simple sample sd/var helpers
static inline double mean_vec(const std::vector<double>& v){
  if (v.empty()) return NA_REAL;
  long double s=0; for(double x: v) s+=x; return (double)(s / (long double)v.size());
}
static inline double sd_vec(const std::vector<double>& v){
  size_t n = v.size(); if (n<2) return NA_REAL;
  long double mu = mean_vec(v), acc=0.0; 
  for(double x: v){ long double d=x-mu; acc += d*d; }
  return std::sqrt((double)(acc / (long double)(n-1)));
}
static inline double var_vec(const std::vector<double>& v){
  double s = sd_vec(v); return R_finite(s)? s*s : NA_REAL;
}

// weighted discrete mode via integer replication counts (ceil(100*w))
static double mode_class_weighted(const std::vector<double>& y,
                                  const std::vector<double>& w){
  std::unordered_map<double, long long> cnt;
  cnt.reserve(y.size()*2);
  for (size_t i=0;i<y.size();++i){
    long long c = (long long) std::ceil(100.0 * w[i]);
    if (c <= 0) continue;
    cnt[y[i]] += c;
  }
  if (cnt.empty()) return NA_REAL;
  double best_val = cnt.begin()->first;
  long long best_cnt = cnt.begin()->second;
  for (auto &kv : cnt){
    if (kv.second > best_cnt){ best_cnt = kv.second; best_val = kv.first; }
  }
  return best_val;
}

// [[Rcpp::export]]
SEXP NNS_distance_cpp(NumericMatrix X,      // features, l x n (NO y.hat column)
                      NumericVector yhat,   // length l
                      NumericVector dest,   // length n (dist.estimate)
                      int k,                // resolved k (k == l if "all")
                      bool use_class) {     // TRUE => mode_class path
  const int l = X.nrow();
  const int n = X.ncol();
  if (yhat.size() != l) stop("yhat length must equal nrow(X)");
  if (dest.size() != n) stop("dist.estimate length must equal ncol(X)");
  
  // --- 1) Min–max rescale each column jointly over X and dest -> [0,1] ----
  std::vector<double> d(n);
  for (int j=0;j<n;++j){
    double cmin = dest[j], cmax = dest[j];
    for (int i=0;i<l;++i){
      double v = X(i,j);
      if (R_finite(v)) { if (v<cmin) cmin=v; if (v>cmax) cmax=v; }
      else { /* keep NA as-is; it will propagate in distance safely */ }
    }
    double range = cmax - cmin;
    double dj = R_finite(dest[j]) ? (range==0.0 ? 0.0 : (dest[j]-cmin)/range) : NA_REAL;
    d[j] = dj;
    if (range != 0.0){
      for (int i=0;i<l;++i){
        double v = X(i,j);
        if (R_finite(v)) X(i,j) = (v - cmin) / range;
        else X(i,j) = NA_REAL;
      }
    } else {
      for (int i=0;i<l;++i) X(i,j) = (double)0.0;
    }
  }
  
  // --- 2) distance Sum_i = sum_j ((x_ij - d_j)^2 + |x_ij - d_j|) -----------
  std::vector<double> S(l, 0.0);
  for (int i=0;i<l;++i){
    long double acc = 0.0L;
    for (int j=0;j<n;++j){
      double a = X(i,j);
      double b = d[j];
      double diff = (R_finite(a) && R_finite(b)) ? (a - b) : 0.0; // treat NA as 0 contribution
      acc += (long double)(diff*diff + std::fabs(diff));
    }
    double si = (double)acc;
    S[i] = (si == 0.0 ? 1e-10 : si); // match your epsilon for zeros
  }
  
  // --- 3) take the best ll = min(k, l) rows via partial_sort ---------------
  int ll = std::min(k, l);
  std::vector<int> idx(l);
  for (int i=0;i<l;++i) idx[i]=i;
  auto cmp = [&](int a, int b){ return S[a] < S[b]; };
  if (ll < l) std::partial_sort(idx.begin(), idx.begin()+ll, idx.end(), cmp);
  else std::sort(idx.begin(), idx.end(), cmp);
  
  idx.resize(ll);
  std::vector<double> Ssel; Ssel.reserve(ll);
  std::vector<double> ysel; ysel.reserve(ll);
  for (int t=0;t<ll;++t){
    int i = idx[t];
    Ssel.push_back(S[i]);
    ysel.push_back(yhat[i]);
  }
  
  // --- 4) k==1 fast path (with tie mode like your R) -----------------------
  if (ll == 1) return wrap(ysel[0]);
  if (k == 1){
    // find all with minimal S and take (discrete) mode; else first
    double smin = *std::min_element(Ssel.begin(), Ssel.end());
    std::vector<double> yties;
    for (int t=0;t<ll;++t) if (Ssel[t] == smin) yties.push_back(ysel[t]);
    if (yties.size() == 1) return wrap(yties[0]);
    // mode among ties
    std::unordered_map<double,int> cnt; cnt.reserve(yties.size()*2);
    for (double v: yties) ++cnt[v];
    double best = yties[0]; int bestc = cnt[best];
    for (auto &kv: cnt) if (kv.second > bestc){ bestc = kv.second; best = kv.first; }
    return wrap(best);
  }
  
  // --- 5) build the 8 weight families --------------------------------------
  std::vector<double> uni(ll, 1.0 / (double)ll);
  
  // t weights: dt(S, df=ll)
  std::vector<double> tw(ll, 0.0);
  for (int i=0;i<ll;++i){
    double dens = ::Rf_dt(Ssel[i], (double)ll, 0);
    tw[i] = R_finite(dens) ? dens : 0.0;
  }
  double twsum = std::accumulate(tw.begin(), tw.end(), 0.0);
  if (twsum>0) for (double &v: tw) v /= twsum; else std::fill(tw.begin(), tw.end(), 0.0);
  
  // empirical 1/S
  std::vector<double> emp(ll, 0.0);
  for (int i=0;i<ll;++i){ double v = Ssel[i]; emp[i] = (v>0) ? 1.0/v : 0.0; }
  double empsum = std::accumulate(emp.begin(), emp.end(), 0.0);
  if (empsum>0) for (double &v: emp) v /= empsum; else std::fill(emp.begin(), emp.end(), 0.0);
  
  // exponential on ranks 1:ll with rate = 1/ll
  std::vector<double> exw(ll, 0.0);
  for (int i=0;i<ll;++i){
    double dens = ::Rf_dexp((double)(i+1), 1.0/(double)ll, 0);
    exw[i] = R_finite(dens) ? dens : 0.0;
  }
  double exsum = std::accumulate(exw.begin(), exw.end(), 0.0);
  if (exsum>0) for (double &v: exw) v /= exsum; else std::fill(exw.begin(), exw.end(), 0.0);
  
  // lognormal on ranks (abs(rev(logpdf))) with sdlog = sd(1:ll)
  std::vector<double> lnorm(ll, 0.0);
  double sd_ranks = NA_REAL;
  if (ll >= 2){
    std::vector<double> ranks(ll); for(int i=0;i<ll;++i) ranks[i]= (double)(i+1);
    sd_ranks = sd_vec(ranks);
  }
  if (R_finite(sd_ranks)){
    for (int i=0;i<ll;++i){
      double lp = ::Rf_dlnorm((double)(i+1), 0.0, sd_ranks, 1 /*log*/);
      lnorm[i] = std::fabs(lp);
    }
    std::reverse(lnorm.begin(), lnorm.end());
  } else {
    std::fill(lnorm.begin(), lnorm.end(), 0.0);
  }
  double lnsum = std::accumulate(lnorm.begin(), lnorm.end(), 0.0);
  if (lnsum>0) for (double &v: lnorm) v /= lnsum; else std::fill(lnorm.begin(), lnorm.end(), 0.0);
  
  // power-law on ranks: (1:i)^(-2)
  std::vector<double> pl(ll, 0.0);
  for (int i=0;i<ll;++i){ double r = (double)(i+1); pl[i] = std::pow(r, -2.0); }
  double plsum = std::accumulate(pl.begin(), pl.end(), 0.0);
  if (plsum>0) for (double &v: pl) v /= plsum; else std::fill(pl.begin(), pl.end(), 0.0);
  
  // normal on S: dnorm(S, 0, sd(S))
  std::vector<double> normw(ll, 0.0);
  double sdS = sd_vec(Ssel);
  if (R_finite(sdS) && sdS>0){
    for (int i=0;i<ll;++i){
      double dens = ::Rf_dnorm4(Ssel[i], 0.0, sdS, 0);
      normw[i] = R_finite(dens) ? dens : 0.0;
    }
    double nsum = std::accumulate(normw.begin(), normw.end(), 0.0);
    if (nsum>0) for (double &v: normw) v /= nsum; else std::fill(normw.begin(), normw.end(), 0.0);
  } // else already zeros
  
  // RBF on S: exp(-S / (2*var(S)))
  std::vector<double> rbf(ll, 0.0);
  double varS = var_vec(Ssel);
  if (R_finite(varS) && varS>0){
    for (int i=0;i<ll;++i) rbf[i] = std::exp(- Ssel[i] / (2.0*varS));
    double rsum = std::accumulate(rbf.begin(), rbf.end(), 0.0);
    if (rsum>0) for (double &v: rbf) v /= rsum; else std::fill(rbf.begin(), rbf.end(), 0.0);
  } // else zeros
  
  // combined weights
  std::vector<double> w(ll, 0.0);
  double tot = 0.0;
  for (int i=0;i<ll;++i){
    double wi = uni[i] + tw[i] + emp[i] + exw[i] + lnorm[i] + pl[i] + normw[i] + rbf[i];
    w[i] = wi; tot += wi;
  }
  if (tot > 0) for (double &v: w) v /= tot; else for (double &v: w) v = 1.0/(double)ll;
  
  // --- 6) prediction ---------------------------------------------------------
  if (!use_class){
    long double dot=0.0L;
    for (int i=0;i<ll;++i) dot += (long double)(ysel[i]*w[i]);
    return wrap((double)dot);
  } else {
    return wrap( mode_class_weighted(ysel, w) );
  }
}


// ---------- NNS_distance_path_cpp / NNS_distance_bulk_cpp (guarded) ----------

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
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(),
              [&dist](int a, int b){ return dist[a] < dist[b]; });
  }
} // anon

// [[Rcpp::export]]
Rcpp::NumericMatrix NNS_distance_path_cpp(const Rcpp::NumericMatrix& RPM,
                                          const Rcpp::NumericVector& yhat,
                                          const Rcpp::NumericMatrix& Xtest,
                                          int kmax,
                                          bool is_class) {
  (void)is_class;
  
  const int n = RPM.nrow();
  const int p = RPM.ncol();
  const int m = Xtest.nrow();
  
  if (n <= 0 || p <= 0 || m <= 0) Rcpp::stop("RPM/Xtest must be non-empty");
  if (yhat.size() != n)          Rcpp::stop("yhat length must equal nrow(RPM)");
  if (Xtest.ncol() != p)         Rcpp::stop("Xtest and RPM must have same number of columns");
  if (kmax < 1)                  Rcpp::stop("kmax must be >= 1");
  if (kmax > n)                  kmax = n; // clamp
  
  Rcpp::NumericMatrix out(m, kmax);
  
  const double* rpm_ptr  = REAL(RPM);
  const double* y_ptr    = REAL(yhat);
  const double* tst_ptr  = REAL(Xtest);
  
  std::vector<double> dist(n);
  std::vector<int>    ord; ord.reserve(n);
  std::vector<double> y_sorted(n), d_sorted(n);
  
  for (int r = 0; r < m; ++r) {
    const double* tr = tst_ptr + static_cast<std::size_t>(r) * p;
    
    compute_distances(rpm_ptr, n, p, tr, dist);
    argsort_by_distance(dist, ord);
    
    for (int i = 0; i < n; ++i) {
      const int j = ord[i];
      y_sorted[i] = y_ptr[j];
      d_sorted[i] = (dist[j] <= 0.0 ? safe_eps() : dist[j]);
    }
    
    double csum_w  = 0.0;
    double csum_yw = 0.0;
    
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
  
  const int n = RPM.nrow();
  const int p = RPM.ncol();
  const int m = Xtest.nrow();
  
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
  std::vector<int>    ord; ord.reserve(n);
  
  for (int r = 0; r < m; ++r) {
    const double* tr = tst_ptr + static_cast<std::size_t>(r) * p;
    
    compute_distances(rpm_ptr, n, p, tr, dist);
    argsort_by_distance(dist, ord);
    
    double csum_w  = 0.0;
    double csum_yw = 0.0;
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
  // read-only inputs
  RMatrix<double> RPM;
  RVector<double> yhat;
  RMatrix<double> Xtest;
  std::vector<double> minRPM, maxRPM;
  int l, n, m, kmax;
  bool is_class;
  
  // rank-only caches (size kmax+1)
  std::vector< std::vector<double> > uniW, expW, lnormW, plW;
  
  // output
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
    std::iota(idx.begin(), idx.end(), 0);
    
    for (std::size_t r = begin; r < end; ++r) {
      // per-feature scales for this row (leak-safe)
      for (int j=0;j<n;++j){
        double t = Xtest(r,j);
        double mn = std::min(minRPM[j], t);
        double mx = std::max(maxRPM[j], t);
        double range = mx - mn;
        invR[j] = (R_finite(range) && range>0.0) ? (1.0/range) : 0.0;
      }
      
      // distances to all RPM rows
      for (int i=0;i<l;++i){
        long double acc=0.0L;
        for (int j=0;j<n;++j){
          double a = RPM(i,j), b = Xtest(r,j);
          if (R_finite(a) && R_finite(b) && invR[j]>0.0){
            double diff = (a-b)*invR[j];
            acc += (long double)(diff*diff + std::fabs(diff));
          }
        }
        double si = (double)acc;
        S[i] = (si==0.0 ? 1e-10 : si);
      }
      
      // select top kmax
      std::iota(idx.begin(), idx.end(), 0);
      auto cmp = [&](int a,int b){ return S[a] < S[b]; };
      if (kmax < l) std::partial_sort(idx.begin(), idx.begin()+kmax, idx.end(), cmp);
      else          std::sort(idx.begin(), idx.end(), cmp);
      
      // stable order within top-kmax for deterministic ties
      auto cmp2 = [&](int a,int b){
        if (S[a] < S[b]) return true;
        if (S[b] < S[a]) return false;
        return a < b;
      };
      std::stable_sort(idx.begin(), idx.begin()+kmax, cmp2);
      
      topS.resize(kmax); topY.resize(kmax);
      for (int t=0;t<kmax;++t){ int i = idx[t]; topS[t]=S[i]; topY[t]=yhat[i]; }
      
      // sweep k = 1..kmax
      for (int k=1;k<=kmax;++k){
        // views into prefix
        const double* Ssel = topS.data();
        const double* Ysel = topY.data();
        
        if (k==1){ out(r, k-1) = Ysel[0]; continue; }
        
        // S-dependent families
        std::vector<double> tw(k,0.0), emp(k,0.0), normw(k,0.0), rbf(k,0.0);
        
        for (int i=0;i<k;++i){
          double d_t = ::Rf_dt(Ssel[i], (double)k, 0);
          tw[i] = R_finite(d_t) ? d_t : 0.0;
          
          double v = Ssel[i];
          emp[i] = (v>0) ? 1.0/v : 0.0;
        }
        double tws = std::accumulate(tw.begin(), tw.end(), 0.0);
        if (tws>0) for(double &v: tw) v/=tws; else std::fill(tw.begin(), tw.end(), 0.0);
        
        double emps = std::accumulate(emp.begin(), emp.end(), 0.0);
        if (emps>0) for(double &v: emp) v/=emps; else std::fill(emp.begin(), emp.end(), 0.0);
        
        double sdS = sd_vec(std::vector<double>(topS.begin(), topS.begin()+k));
        if (R_finite(sdS) && sdS>0){
          for (int i=0;i<k;++i){
            double d_n = ::Rf_dnorm4(Ssel[i], 0.0, sdS, 0);
            normw[i] = R_finite(d_n) ? d_n : 0.0;
          }
          double ns = std::accumulate(normw.begin(), normw.end(), 0.0);
          if (ns>0) for(double &v: normw) v/=ns; else std::fill(normw.begin(), normw.end(), 0.0);
        }
        
        double vS = var_vec(std::vector<double>(topS.begin(), topS.begin()+k));
        if (R_finite(vS) && vS>0){
          for (int i=0;i<k;++i) rbf[i] = std::exp(- Ssel[i] / (2.0*vS));
          double rs = std::accumulate(rbf.begin(), rbf.end(), 0.0);
          if (rs>0) for(double &v: rbf) v/=rs; else std::fill(rbf.begin(), rbf.end(), 0.0);
        }
        
        // combine rank-only + S families + uniform
        long double dot = 0.0L, tot = 0.0L;
        for (int i=0;i<k;++i){
          double wi = uniW[k][i] + expW[k][i] + lnormW[k][i] + plW[k][i]
          + tw[i] + emp[i] + normw[i] + rbf[i];
          tot += wi;
          if (!is_class) dot += (long double)(Ysel[i] * wi);
        }
        double invTot = (tot>0.0) ? (1.0/(double)tot) : (1.0/(double)k);
        
        if (!is_class){
          out(r, k-1) = (tot>0.0) ? (double)(dot * invTot) : (double)(
            std::accumulate(topY.begin(), topY.begin()+k, 0.0) / (double)k
          );
        } else {
          std::vector<double> w(k);
          if (tot>0.0) for (int i=0;i<k;++i) w[i] = (uniW[k][i]+expW[k][i]+lnormW[k][i]+plW[k][i]+tw[i]+emp[i]+normw[i]+rbf[i]) * invTot;
          else          std::fill(w.begin(), w.end(), 1.0/(double)k);
          out(r, k-1) = mode_class_weighted(std::vector<double>(topY.begin(), topY.begin()+k), w);
        }
      } // k
    }   // r
  }     // operator()
};

// [[Rcpp::export]]
NumericMatrix NNS_distance_path_parallel_cpp(NumericMatrix RPM,
                                             NumericVector yhat,
                                             NumericMatrix Xtest,
                                             int kmax,
                                             bool is_class,
                                             int nthreads = -1) {
  const int l = RPM.nrow(), n = RPM.ncol(), m = Xtest.nrow();
  if (yhat.size()!=l) stop("yhat length must equal nrow(RPM)");
  if (kmax<=0) kmax=l;
  if (kmax>l) kmax=l;
  
  // RPM min/max per feature (shared, read-only)
  std::vector<double> minRPM(n, R_PosInf), maxRPM(n, R_NegInf);
  for (int j=0;j<n;++j){
    for (int i=0;i<l;++i){
      double v = RPM(i,j);
      if (R_finite(v)) { if(v<minRPM[j]) minRPM[j]=v; if(v>maxRPM[j]) maxRPM[j]=v; }
    }
    if (!R_finite(minRPM[j])) { minRPM[j]=0.0; maxRPM[j]=0.0; }
  }
  
  // rank-only caches up to kmax
  std::vector<std::vector<double>> uniW(kmax+1), expW(kmax+1), lnormW(kmax+1), plW(kmax+1);
  for (int k=1;k<=kmax;++k){
    uniW[k].assign(k, 1.0/(double)k);
    
    std::vector<double> ex(k);
    for (int r=1;r<=k;++r) ex[r-1] = ::Rf_dexp((double)r, 1.0/(double)k, 0);
    double exs = std::accumulate(ex.begin(), ex.end(), 0.0);
    if (exs>0) for (double &v: ex) v/=exs; else std::fill(ex.begin(), ex.end(), 0.0);
    expW[k] = std::move(ex);
    
    std::vector<double> pl(k);
    for (int r=1;r<=k;++r) pl[r-1] = std::pow((double)r, -2.0);
    double pls = std::accumulate(pl.begin(), pl.end(), 0.0);
    if (pls>0) for (double &v: pl) v/=pls; else std::fill(pl.begin(), pl.end(), 0.0);
    plW[k] = std::move(pl);
    
    std::vector<double> ln(k,0.0);
    if (k>=2){
      double sdlog = std::sqrt(((double)k*(double)k - 1.0)/12.0); // sd(1..k)
      for (int r=1;r<=k;++r){ double lp = ::Rf_dlnorm((double)r, 0.0, sdlog, 1); ln[r-1]=std::fabs(lp); }
      std::reverse(ln.begin(), ln.end());
      double lns = std::accumulate(ln.begin(), ln.end(), 0.0);
      if (lns>0) for (double &v: ln) v/=lns; else std::fill(ln.begin(), ln.end(), 0.0);
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

