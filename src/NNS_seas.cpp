// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <set>
#include <algorithm>
#include <cfloat>

using namespace Rcpp;


// --- small utilities (no plotting here) ---
inline bool any_na_or_inf(const NumericVector& x){
  int n = x.size();
  for(int i=0;i<n;++i){
    double v = x[i];
    if (Rcpp::NumericVector::is_na(v) || !R_finite(v)) return true;
  }
  return false;
}

inline double vec_mean(const NumericVector& x){
  double s = 0.0; int n = x.size();
  for(int i=0;i<n;++i) s += x[i];
  return s / (double)n;
}

// c4(k): E[s] = c4(k) * sigma for k normal observations. Dividing by it removes
// the small-sample downward bias of the sd, which otherwise favours the lags
// with the shortest components (#69).
inline double c4(int k){
  if (k < 2) return 1.0;
  const double kd = static_cast<double>(k);
  return std::sqrt(2.0 / (kd - 1.0)) * std::exp(std::lgamma(kd / 2.0) - std::lgamma((kd - 1.0) / 2.0));
}

inline double vec_sd(const NumericVector& x){
  int n = x.size();
  if(n < 2) return NA_REAL;
  double m = vec_mean(x);
  double ss = 0.0;
  for(int i=0;i<n;++i){
    double d = x[i]-m; ss += d*d;
  }
  return std::sqrt(ss / (double)(n-1));
}

// [[Rcpp::export]]
Rcpp::List NNS_seas_cpp(NumericVector variable,
                        Nullable<IntegerVector> modulo = R_NilValue,
                        bool mod_only = true){
  if (variable.size() == 0)  stop("Variable must be numeric and non-empty");
  if (any_na_or_inf(variable)) stop("You have some missing or infinite values, please address.");
  
  const int n = variable.size();
  if (n < 5){
    DataFrame M = DataFrame::create(
      _["Period"] = IntegerVector::create(0),
      _["Coefficient.of.Variation"] = NumericVector::create(0.0),
      _["Variable.Coefficient.of.Variation"] = NumericVector::create(0.0)
    );
    return List::create(
      _["all.periods"] = M,
      _["best.period"] = 0,
      _["periods"]     = IntegerVector::create(0)
    );
  }
  
  const int half_n = n/2;
  const double mean_var = vec_mean(variable);
  const double sd_var = vec_sd(variable) / c4(n);
  // CV is meaningless for a mean that is zero up to rounding; compare sds instead
  const double tol = std::sqrt(DBL_EPSILON);
  const bool use_cv = std::fabs(mean_var) > tol * sd_var;
  const double var_cov = use_cv ? std::fabs(sd_var / mean_var) : sd_var;
  
  // Minimum component length (#69): every component series of lag m must have
  // at least MIN_COMPONENT_LENGTH points, i.e. m <= n / MIN_COMPONENT_LENGTH.
  // With 2-3 points per component the sample sd is biased low (E[s] = c4(k) sigma,
  // c4(2) = 0.80, c4(3) = 0.89), so the longest lags looked tightest even in
  // pure noise and were picked as the best period.
  const int MIN_COMPONENT_LENGTH = 5;
  const int max_m = std::min(half_n, n / MIN_COMPONENT_LENGTH);

  std::vector<int> periods_vec;
  std::vector<double> cvmean_vec;
  const double* x = variable.begin();
  auto score = [&](int m){
    double wsum = 0.0;
    for(int j=0;j<m;++j){
      int len = (n-1-j)/m + 1;
      double mu = 0.0, ss = 0.0;
      for(int i=j;i<n;i+=m) mu += x[i];
      mu /= (double)len;
      for(int i=j;i<n;i+=m){ double d = x[i]-mu; ss += d*d; }
      double sd = std::sqrt(ss / (double)(len-1)) / c4(len);
      // a constant component has CV 0; a near-zero mean caps CV at 1/tol, where the series switches to sd
      wsum += len * (use_cv ? (sd == 0.0 ? 0.0 : sd / std::max(std::fabs(mu), tol * sd)) : sd);
    }
    return wsum / (double)n;
  };
  // A period the user adds through modulo is scored even beyond the scan limit,
  // so ARMA.seas.weighting doesn't give it zero weight.
  auto added_cv = [&](int m){
    if (m > max_m && m <= half_n){ double c = score(m); if (c < var_cov) return c; }
    return var_cov;
  };
  for(int m=1;m<=max_m;++m){
    double cv_m = score(m);
    if (cv_m < var_cov){
      periods_vec.push_back(m);
      cvmean_vec.push_back(cv_m);
    }
  }
  
  IntegerVector Period;
  NumericVector CoefVar;
  NumericVector VarCoefVar;
  
  if(!periods_vec.empty()){
    int m = (int)periods_vec.size();
    Period = IntegerVector(m);
    CoefVar = NumericVector(m);
    VarCoefVar = NumericVector(m);
    for(int k=0;k<m;++k){ Period[k]=periods_vec[k]; CoefVar[k]=cvmean_vec[k]; VarCoefVar[k]=var_cov; }
    IntegerVector ord = seq(0, m-1);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
    Period = Period[ord]; CoefVar = CoefVar[ord]; VarCoefVar = VarCoefVar[ord];
  }else{
    Period = IntegerVector::create(1);
    CoefVar = NumericVector::create(var_cov);
    VarCoefVar = NumericVector::create(var_cov);
  }
  
  // modulo handling (no plotting)
  if (modulo.isNotNull()){
    IntegerVector mods(modulo.get());
    std::set<int> per_set;
    for(int i=0;i<Period.size();++i){
      for(int j=0;j<mods.size();++j){
        int m = mods[j]; if(m<=0) continue;
        int minus = Period[i] - (Period[i] % m);
        int plus  = Period[i] + (m - (Period[i] % m));
        if (minus > 0) per_set.insert(minus);
        if (plus  > 0) per_set.insert(plus);
      }
    }
    if (mod_only){
      std::set<int> curr;
      for(int i=0;i<Period.size();++i) curr.insert(Period[i]);
      std::vector<int> keptP; std::vector<double> keptCV;
      for(int i=0;i<Period.size();++i){
        if (per_set.count(Period[i])){ keptP.push_back(Period[i]); keptCV.push_back(CoefVar[i]); }
      }
      for(int s: per_set) if(!curr.count(s)){ keptP.push_back(s); keptCV.push_back(added_cv(s)); }
      if(keptP.empty()){
        Period   = IntegerVector::create(1);
        CoefVar  = NumericVector::create(var_cov);
        VarCoefVar = NumericVector::create(var_cov);
      } else {
        int m = (int)keptP.size();
        Period = IntegerVector(m); CoefVar = NumericVector(m); VarCoefVar = NumericVector(m);
        for(int i=0;i<m;++i){ Period[i]=keptP[i]; CoefVar[i]=keptCV[i]; VarCoefVar[i]=var_cov; }
        IntegerVector ord = seq(0, m-1);
        std::sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
        Period = Period[ord]; CoefVar = CoefVar[ord]; VarCoefVar = VarCoefVar[ord];
      }
    } else {
      per_set.insert(1);
      std::set<int> curr;
      for(int i=0;i<Period.size();++i) curr.insert(Period[i]);
      std::vector<int> add;
      for(int s: per_set) if(!curr.count(s)) add.push_back(s);
      if(!add.empty()){
        int oldm = Period.size(), addm = (int)add.size();
        IntegerVector P2(oldm+addm); NumericVector CV2(oldm+addm); NumericVector VCV2(oldm+addm);
        for(int i=0;i<oldm;++i){ P2[i]=Period[i]; CV2[i]=CoefVar[i]; VCV2[i]=VarCoefVar[i]; }
        for(int i=0;i<addm;++i){ P2[oldm+i]=add[i]; CV2[oldm+i]=added_cv(add[i]); VCV2[oldm+i]=var_cov; }
        Period=P2; CoefVar=CV2; VarCoefVar=VCV2;
        IntegerVector ord = seq(0, Period.size()-1);
        std::sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
        Period = Period[ord]; CoefVar = CoefVar[ord]; VarCoefVar = VarCoefVar[ord];
      }
    }
  }
  
  // strict cap: Period < n/2
  {
    std::vector<int> P; std::vector<double> CV; std::vector<double> VCV;
    for(int i=0;i<Period.size();++i){
      if(Period[i] < n/2){ P.push_back(Period[i]); CV.push_back(CoefVar[i]); VCV.push_back(VarCoefVar[i]); }
    }
    if(!P.empty()){
      Period = IntegerVector(P.begin(), P.end());
      CoefVar = NumericVector(CV.begin(), CV.end());
      VarCoefVar = NumericVector(VCV.begin(), VCV.end());
      IntegerVector ord = seq(0, Period.size()-1);
      std::sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
      Period = Period[ord]; CoefVar = CoefVar[ord]; VarCoefVar = VarCoefVar[ord];
    } else {
      Period = IntegerVector::create(1);
      CoefVar = NumericVector::create(var_cov);
      VarCoefVar = NumericVector::create(var_cov);
    }
  }
  
  DataFrame M = DataFrame::create(
    _["Period"] = Period,
    _["Coefficient.of.Variation"] = CoefVar,
    _["Variable.Coefficient.of.Variation"] = VarCoefVar
  );
  
  return List::create(
    _["all.periods"] = M,
    _["best.period"] = Period[0],
    _["periods"]     = Period
  );
}

