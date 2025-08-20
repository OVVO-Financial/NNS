#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
using namespace Rcpp;
using namespace RcppParallel;

// ---------- helpers ----------
static inline double mean_vec(const std::vector<double>& v){
  if (v.empty()) return NA_REAL;
  long double s=0; for(double x: v) s+=x; return (double)(s / (long double)v.size());
}
static inline double sd_vec(const std::vector<double>& v){
  size_t n=v.size(); if (n<2) return NA_REAL;
  long double mu=mean_vec(v), acc=0.0L; 
  for(double x: v){ long double d=x-mu; acc+=d*d; }
  return std::sqrt((double)(acc / (long double)(n-1)));
}
static inline double var_vec(const std::vector<double>& v){ double s=sd_vec(v); return R_finite(s)? s*s : NA_REAL; }

// weighted discrete mode via ceil(100*w)
static inline double mode_class_weighted(const std::vector<double>& y,
                                         const std::vector<double>& w){
  std::unordered_map<double, long long> cnt;
  cnt.reserve(y.size()*2);
  for (size_t i=0;i<y.size();++i){
    long long c = (long long) std::ceil(100.0 * w[i]);
    if (c <= 0) continue;
    cnt[y[i]] += c;
  }
  if (cnt.empty()) return NA_REAL;
  auto it = cnt.begin();
  double best = it->first; long long mx = it->second; ++it;
  for (; it != cnt.end(); ++it) if (it->second > mx) { mx = it->second; best = it->first; }
  return best;
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
