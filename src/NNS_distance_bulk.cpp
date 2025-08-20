#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace Rcpp;

// small helpers
static inline double mean_vec(const std::vector<double>& v){
  if (v.empty()) return NA_REAL;
  long double s=0; for(double x: v) s+=x; return (double)(s / (long double)v.size());
}
static inline double sd_vec(const std::vector<double>& v){
  size_t n=v.size(); if (n<2) return NA_REAL;
  long double mu=mean_vec(v), acc=0.0L; 
  for(double x: v){ long double d=x-mu; acc += d*d; }
  return std::sqrt((double)(acc / (long double)(n-1)));
}
static inline double var_vec(const std::vector<double>& v){
  double s=sd_vec(v); return R_finite(s)? s*s : NA_REAL;
}

// weighted mode via integer replication ~ ceil(100*w)
static double mode_class_weighted(const std::vector<double>& y,
                                  const std::vector<double>& w) {
  std::unordered_map<double, long long> cnt;
  cnt.reserve(y.size() * 2);
  for (size_t i = 0; i < y.size(); ++i) {
    long long c = (long long) std::ceil(100.0 * w[i]);
    if (c <= 0) continue;
    cnt[y[i]] += c;
  }
  if (cnt.empty()) return NA_REAL;
  
  auto it = cnt.begin();
  double best = it->first;
  long long bestc = it->second;
  ++it;
  
  for (; it != cnt.end(); ++it) {
    if (it->second > bestc) {
      bestc = it->second;
      best  = it->first;
    }
  }
  return best;
}

// [[Rcpp::export]]
NumericVector NNS_distance_bulk_cpp(NumericMatrix RPM,   // l x n (features only)
                                    NumericVector yhat,  // length l
                                    NumericMatrix Xtest, // m x n_test (same or overlapping cols)
                                    int k,
                                    bool is_class) {
  const int l = RPM.nrow();
  const int n = RPM.ncol();
  const int m = Xtest.nrow();
  if (yhat.size() != l) stop("yhat length must equal nrow(RPM)");
  if (k <= 0) k = l;
  const int ll = std::min(k, l);
  
  // Precompute RPM-only min/max per feature (used, then expanded by each test row value)
  std::vector<double> minRPM(n, R_PosInf), maxRPM(n, R_NegInf);
  for (int j=0;j<n;++j){
    for (int i=0;i<l;++i){
      double v = RPM(i,j);
      if (R_finite(v)) {
        if (v < minRPM[j]) minRPM[j] = v;
        if (v > maxRPM[j]) maxRPM[j] = v;
      }
    }
    if (!R_finite(minRPM[j])) { minRPM[j] = 0.0; maxRPM[j] = 0.0; } // empty/NA guard
  }
  
  NumericVector out(m);
  
  // buffers reused per test row
  std::vector<double> S(l), Ssel, ysel, ranks, uni, tw, emp, exw, lnorm, pl, normw, rbf, w;
  
  for (int r=0; r<m; ++r){
    // 1) per-feature scale using this test row only
    std::vector<double> invR(n, 0.0);
    for (int j=0;j<n;++j){
      double t = Xtest(r,j);
      double mn = std::min(minRPM[j], t);
      double mx = std::max(maxRPM[j], t);
      double range = mx - mn;
      invR[j] = (R_finite(range) && range > 0.0) ? (1.0 / range) : 0.0;
    }
    
    // 2) distances to all RPM rows
    for (int i=0;i<l;++i){
      long double acc = 0.0L;
      for (int j=0;j<n;++j){
        double a = RPM(i,j);
        double b = Xtest(r,j);
        if (R_finite(a) && R_finite(b) && invR[j] > 0.0){
          double diff = (a - b) * invR[j];
          acc += (long double)(diff*diff + std::fabs(diff));
        }
        // else NA or zero-range -> contributes 0, matching "constant column" behavior
      }
      double si = (double)acc;
      S[i] = (si == 0.0 ? 1e-10 : si);
    }
    
    // 3) pick top-ll via partial sort
    std::vector<int> idx(l); for (int i=0;i<l;++i) idx[i]=i;
    auto cmp = [&](int a, int b){ return S[a] < S[b]; };
    if (ll < l) std::partial_sort(idx.begin(), idx.begin()+ll, idx.end(), cmp);
    else std::sort(idx.begin(), idx.end(), cmp);
    idx.resize(ll);
    
    Ssel.clear(); Ssel.reserve(ll);
    ysel.clear(); ysel.reserve(ll);
    for (int t=0;t<ll;++t){ Ssel.push_back(S[idx[t]]); ysel.push_back(yhat[idx[t]]); }
    
    // k==1 shortcut (preserve tie logic: mode among equal minima)
    if (ll == 1) { out[r] = ysel[0]; continue; }
    if (k == 1){
      double smin = *std::min_element(Ssel.begin(), Ssel.end());
      std::vector<double> ties;
      for (int t=0;t<ll;++t) if (Ssel[t]==smin) ties.push_back(ysel[t]);
      if (ties.size()==1) { out[r]=ties[0]; continue; }
      // discrete mode among ties
      std::unordered_map<double,int> cnt; cnt.reserve(ties.size()*2);
      for (double v: ties) ++cnt[v];
      double best = ties[0]; int bestc = cnt[best];
      for (auto &kv: cnt) if (kv.second > bestc){ bestc = kv.second; best = kv.first; }
      out[r] = best; continue;
    }
    
    // 4) build weight families (same as your R)
    uni.assign(ll, 1.0 / (double)ll);
    
    tw.assign(ll, 0.0);
    for (int i=0;i<ll;++i) {
      double dens = ::Rf_dt(Ssel[i], (double)ll, 0);
      tw[i] = R_finite(dens) ? dens : 0.0;
    }
    double twsum = std::accumulate(tw.begin(), tw.end(), 0.0);
    if (twsum>0) for (double &v: tw) v/=twsum; else std::fill(tw.begin(), tw.end(), 0.0);
    
    emp.assign(ll, 0.0);
    for (int i=0;i<ll;++i) { double v=Ssel[i]; emp[i]=(v>0)?1.0/v:0.0; }
    double empsum = std::accumulate(emp.begin(), emp.end(), 0.0);
    if (empsum>0) for (double &v: emp) v/=empsum; else std::fill(emp.begin(), emp.end(), 0.0);
    
    exw.assign(ll, 0.0);
    for (int i=0;i<ll;++i){
      double dens = ::Rf_dexp((double)(i+1), 1.0/(double)ll, 0);
      exw[i] = R_finite(dens) ? dens : 0.0;
    }
    double exsum = std::accumulate(exw.begin(), exw.end(), 0.0);
    if (exsum>0) for (double &v: exw) v/=exsum; else std::fill(exw.begin(), exw.end(), 0.0);
    
    lnorm.assign(ll, 0.0);
    if (ll >= 2){
      std::vector<double> ranksv(ll); for (int i=0;i<ll;++i) ranksv[i]=(double)(i+1);
      double sdr = sd_vec(ranksv);
      if (R_finite(sdr) && sdr>0){
        for (int i=0;i<ll;++i){
          double lp = ::Rf_dlnorm((double)(i+1), 0.0, sdr, 1);
          lnorm[i] = std::fabs(lp);
        }
        std::reverse(lnorm.begin(), lnorm.end());
        double lns = std::accumulate(lnorm.begin(), lnorm.end(), 0.0);
        if (lns>0) for (double &v: lnorm) v/=lns; else std::fill(lnorm.begin(), lnorm.end(), 0.0);
      }
    }
    
    pl.assign(ll, 0.0);
    for (int i=0;i<ll;++i){ double rnk=(double)(i+1); pl[i]=std::pow(rnk,-2.0); }
    double plsum = std::accumulate(pl.begin(), pl.end(), 0.0);
    if (plsum>0) for (double &v: pl) v/=plsum; else std::fill(pl.begin(), pl.end(), 0.0);
    
    normw.assign(ll, 0.0);
    { double sdS = sd_vec(Ssel);
      if (R_finite(sdS) && sdS>0){
        for (int i=0;i<ll;++i){
          double dens = ::Rf_dnorm4(Ssel[i], 0.0, sdS, 0);
          normw[i] = R_finite(dens) ? dens : 0.0;
        }
        double nsum = std::accumulate(normw.begin(), normw.end(), 0.0);
        if (nsum>0) for (double &v: normw) v/=nsum; else std::fill(normw.begin(), normw.end(), 0.0);
      }
    }
    
    rbf.assign(ll, 0.0);
    { double vS = var_vec(Ssel);
      if (R_finite(vS) && vS>0){
        for (int i=0;i<ll;++i) rbf[i] = std::exp(- Ssel[i] / (2.0*vS));
        double rsum = std::accumulate(rbf.begin(), rbf.end(), 0.0);
        if (rsum>0) for (double &v: rbf) v/=rsum; else std::fill(rbf.begin(), rbf.end(), 0.0);
      }
    }
    
    w.assign(ll, 0.0);
    double tot=0.0;
    for (int i=0;i<ll;++i){
      double wi = uni[i]+tw[i]+emp[i]+exw[i]+lnorm[i]+pl[i]+normw[i]+rbf[i];
      w[i]=wi; tot+=wi;
    }
    if (tot>0) for (double &v: w) v/=tot; else for (double &v: w) v = 1.0/(double)ll;
    
    // 5) prediction
    if (!is_class){
      long double dot=0.0L; for (int i=0;i<ll;++i) dot += (long double)(ysel[i]*w[i]);
      out[r] = (double)dot;
    } else {
      out[r] = mode_class_weighted(ysel, w);
    }
  }
  
  return out;
}

