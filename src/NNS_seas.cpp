// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <set>
#include <algorithm>

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

// lag-1 Pearson autocorrelation
inline double acf1(const NumericVector& x){
  int n = x.size();
  if(n < 2) return NA_REAL;
  double m = vec_mean(x);
  double num = 0.0;
  double den = 0.0;
  for(int t=1;t<n;++t) num += (x[t]-m)*(x[t-1]-m);
  for(int t=0;t<n;++t){ double d=x[t]-m; den += d*d; }
  if(den == 0.0) return NA_REAL;
  return num / den;
}

inline double cv_or_fallback(const NumericVector& x, bool use_cv, double var_cov){
  int n = x.size();
  if(n < 2) return var_cov;
  double z;
  if(use_cv){
    double m = vec_mean(x);
    double s = vec_sd(x);
    z = std::fabs(s / m);
  } else {
    double a1 = acf1(x);
    if (NumericVector::is_na(a1)) return var_cov;
    z = std::pow(std::fabs(a1), -1.0);
  }
  if(!R_finite(z)) return var_cov;
  return z;
}

inline IntegerVector rev_step_indices(int n, int step){
  int len = (n-1) / step + 1;
  IntegerVector out(len);
  int v = n;
  for(int k=0;k<len;++k, v -= step) out[k] = v; // 1-based
  return out;
}

inline NumericVector take_by_index(const NumericVector& x, const IntegerVector& idx){
  int m = idx.size();
  NumericVector out(m);
  for(int i=0;i<m;++i){
    int j = idx[i]-1;
    out[i] = (j>=0 && j<x.size()) ? x[j] : NA_REAL;
  }
  return out;
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
  
  NumericVector variable_1(n-1);
  for(int i=0;i<n-1;++i) variable_1[i] = variable[i];
  NumericVector variable_2 = (n-1>=2) ? NumericVector(n-2) : NumericVector(0);
  for(int i=0;i<variable_2.size();++i) variable_2[i] = variable_1[i];
  
  const int half_n = n/2;
  const double mean_var = vec_mean(variable);
  const bool use_cv = (mean_var != 0.0);
  double var_cov = use_cv ? std::fabs(vec_sd(variable)/mean_var) : std::pow(std::fabs(acf1(variable)), -1.0);
  if(!R_finite(var_cov)) var_cov = R_PosInf;
  
  NumericVector out(half_n), out1(half_n), out2(half_n);
  IntegerVector inst(half_n), inst1(half_n), inst2(half_n);
  
  const int n1 = n-1;
  const int n2 = variable_2.size();
  
  for(int i=1;i<=half_n;++i){
    IntegerVector idx  = rev_step_indices(n , i);
    IntegerVector idx1 = rev_step_indices(n1, i);
    IntegerVector idx2 = (n2>0) ? rev_step_indices(n2, i) : IntegerVector(0);
    
    double t  = cv_or_fallback(take_by_index(variable  , idx ), use_cv, var_cov);
    double t1 = cv_or_fallback(take_by_index(variable_1, idx1), use_cv, var_cov);
    double t2 = cv_or_fallback(take_by_index(variable_2, idx2), use_cv, var_cov);
    
    if (t  <= var_cov){ inst[i-1]  = i; out[i-1]  = t;  }
    if (t1 <= var_cov){ inst1[i-1] = i; out1[i-1] = t1; }
    if (t2 <= var_cov){ inst2[i-1] = i; out2[i-1] = t2; }
  }
  
  // build passing set and average CV across the three staggered series
  std::vector<int> periods_vec;
  std::vector<double> cvmean_vec;
  for(int i=0;i<half_n;++i){
    if (inst[i] > 0 && inst1[i] > 0 && inst2[i] > 0){
      periods_vec.push_back(inst[i]);
      cvmean_vec.push_back( (out[i] + out1[i] + out2[i]) / 3.0 );
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
      for(int s: per_set) if(!curr.count(s)){ keptP.push_back(s); keptCV.push_back(var_cov); }
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
        for(int i=0;i<addm;++i){ P2[oldm+i]=add[i]; CV2[oldm+i]=var_cov; VCV2[oldm+i]=var_cov; }
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

