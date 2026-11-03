// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <set>
#include <algorithm>
#include <cfloat>
#include <numeric>

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

// [[Rcpp::export]]
Rcpp::List NNS_seas_cpp(NumericVector variable,
                        Nullable<IntegerVector> modulo = R_NilValue,
                        bool mod_only = true,
                        Nullable<int> max_period = R_NilValue){
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
  const double* x = variable.begin();

  // Dispersion is the Gini mean difference, the degree-1 lower plus upper
  // partial moment about each observation in turn: as the U-statistic
  // sum_{i != j} |x_i - x_j| / (k (k - 1)) it is unbiased under iid sampling with a finite first moment
  // at any component length k, so no normal-theory factor such as c4(k) is
  // needed (#97). One sort of the series gives every component's order
  // statistics: walking the sorted values, each lands in its component at
  // that component's next rank.
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return x[a] < x[b]; });

  // The coefficient of variation needs a ratio scale: it is used only when the
  // series stays on one side of zero, and dispersions are compared otherwise.
  const bool use_cv = (x[ord[0]] > 0.0) || (x[ord[n - 1]] < 0.0);

  // Minimum component length (#69): every component series of lag m must have
  // at least MIN_COMPONENT_LENGTH points, i.e. m <= n / MIN_COMPONENT_LENGTH.
  // max_period caps the scan further; the scan is O(n) per lag.
  const int MIN_COMPONENT_LENGTH = 5;
  int max_m = std::min(half_n, n / MIN_COMPONENT_LENGTH);
  if (max_period.isNotNull()) max_m = std::max(1, std::min(max_m, as<int>(max_period)));

  const int buf = half_n + 1;
  std::vector<int> rank(buf), len(buf);
  std::vector<double> sum(buf), acc(buf), first(buf), last(buf);
  auto score = [&](int m){
    const double inv_m = 1.0 / (double)m;
    for(int j=0;j<m;++j){ rank[j] = 0; len[j] = (n-1-j)/m + 1; sum[j] = 0.0; acc[j] = 0.0; }
    for(int t=0;t<n;++t){
      const int i = ord[t];
      int j = i - m * (int)((double)i * inv_m);   // i mod m without a division
      if (j >= m) j -= m; else if (j < 0) j += m;
      const double v = x[i];
      if (rank[j] == 0) first[j] = v;
      last[j] = v;
      // The weights sum to 0, so measuring from the component minimum avoids cancelling a large level.
      acc[j] += (2.0 * rank[j] + 1.0 - len[j]) * (v - first[j]);
      sum[j] += v;
      ++rank[j];
    }
    double wsum = 0.0;
    for(int j=0;j<m;++j){
      const double k = (double)len[j];
      // a constant component (smallest value = largest) has dispersion 0
      const double g = (len[j] < 2 || first[j] == last[j]) ? 0.0 : std::max(0.0, 2.0 * acc[j] / (k * (k - 1.0)));
      wsum += k * (use_cv ? (g == 0.0 ? 0.0 : g / std::fabs(sum[j] / k)) : g);
    }
    return wsum / (double)n;
  };
  // The series itself is the single component of lag 1.
  const double var_cov = score(1);

  std::vector<int> periods_vec;
  std::vector<double> cvmean_vec;
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
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
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
        std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
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
        std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
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
      std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return CoefVar[a] < CoefVar[b]; });
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

