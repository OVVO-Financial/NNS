#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <unordered_map>
using namespace Rcpp;

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
