// SD_prefix_refactor.cpp
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <string>
using namespace Rcpp;
using namespace RcppParallel;

// =====================================================================
// Per-column precompute: sorted values, prefix sums, basic stats
// =====================================================================
struct ColPre {
  std::vector<double> vals;      // sorted ascending, length m
  std::vector<double> P1;        // prefix sum of vals; length m+1, P1[0]=0
  std::vector<double> P2;        // prefix sum of vals^2; length m+1
  double S1{0.0}, S2{0.0};
  double mn{R_PosInf}, mean{NA_REAL};
  int m{0};
};

static ColPre precompute_col(const NumericMatrix& X, int j){
  ColPre c; c.m = X.nrow();
  c.vals.resize(c.m);
  for(int i=0;i<c.m;++i) c.vals[i] = X(i,j);
  std::sort(c.vals.begin(), c.vals.end());
  c.P1.assign(c.m+1, 0.0);
  c.P2.assign(c.m+1, 0.0);
  for(int k=1;k<=c.m;++k){
    double v = c.vals[k-1];
    c.P1[k] = c.P1[k-1] + v;
    c.P2[k] = c.P2[k-1] + v*v;
  }
  c.S1   = c.P1[c.m];
  c.S2   = c.P2[c.m];
  c.mn   = c.vals.front();
  c.mean = c.S1 / double(c.m);
  return c;
}

static ColPre precompute_vec(const NumericVector& x){
  ColPre c; c.m = x.size();
  c.vals.assign(x.begin(), x.end());
  std::sort(c.vals.begin(), c.vals.end());
  c.P1.assign(c.m+1, 0.0);
  c.P2.assign(c.m+1, 0.0);
  for(int k=1;k<=c.m;++k){
    double v = c.vals[k-1];
    c.P1[k] = c.P1[k-1] + v;
    c.P2[k] = c.P2[k-1] + v*v;
  }
  c.S1   = c.P1[c.m];
  c.S2   = c.P2[c.m];
  c.mn   = c.vals.front();
  c.mean = c.S1 / double(c.m);
  return c;
}

inline bool identical_samples(const ColPre& a, const ColPre& b){
  if(a.m != b.m) return false;
  for(int i=0;i<a.m;++i) if(a.vals[i] != b.vals[i]) return false;
  return true;
}

// =====================================================================
// O(1) evaluators from prefix sums
// =====================================================================
inline void lpm_upm_deg1(const ColPre& c, int k, double t, double& L1, double& U1){
  // L1 = mean(max(t - x,0)) = (k*t - P1[k]) / m
  // U1 = mean(max(x - t,0)) = (S1 - P1[k] - (m-k)*t) / m
  double m = double(c.m);
  L1 = (k*t - c.P1[k]) / m;
  U1 = ((c.S1 - c.P1[k]) - (c.m - k)*t) / m;
}
inline double lpm_deg2(const ColPre& c, int k, double t){
  // L2 = mean(max(t-x,0)^2) = (k*t^2 - 2t*P1[k] + P2[k]) / m
  double m = double(c.m);
  return (k*t*t - 2.0*t*c.P1[k] + c.P2[k]) / m;
}

// Walk merged grid and apply a functor at each threshold t
template <class F>
inline void for_each_threshold(const ColPre& a, const ColPre& b, F f){
  int ia=0, ib=0, m=a.m; // assume same m
  while(ia<m || ib<m){
    double next_a = (ia<m ? a.vals[ia] : R_PosInf);
    double next_b = (ib<m ? b.vals[ib] : R_PosInf);
    double t = (next_a < next_b ? next_a : next_b);
    while(ia<m && a.vals[ia] <= t) ++ia; // k_a = ia
    while(ib<m && b.vals[ib] <= t) ++ib; // k_b = ib
    f(t, ia, ib);
  }
}

// =====================================================================
// Pairwise dominance via prefix sums (O(m) per pair)
// degree: 1=FSD, 2=SSD, 3=TSD. 'discrete' only matters for FSD.
// Returns 1 iff X dominates Y, else 0.
// =====================================================================
inline int sd_dom_pair(const ColPre& X, const ColPre& Y, int degree, bool discrete){
  if (degree==1){ // FSD
    if (!(X.mn >= Y.mn)) return 0;                 // FSD gate
    if (identical_samples(X, Y)) return 0;         // identical series -> 0 (matches R's identical(LPM_x, LPM_y))
    
    bool x_gt_y = false;
    int deg = (discrete ? 0 : 1); // discrete->0, continuous->1
    for_each_threshold(X, Y, [&](double t, int kx, int ky){
      double Rx, Ry;
      if (deg==0){
        // L0/(L0+U0) == ECDF
        Rx = double(kx)/double(X.m);
        Ry = double(ky)/double(Y.m);
      } else {
        double Lx, Ux, Ly, Uy;
        lpm_upm_deg1(X, kx, t, Lx, Ux);
        lpm_upm_deg1(Y, ky, t, Ly, Uy);
        double Ax = Lx+Ux, Ay = Ly+Uy;
        Rx = (Ax>0.0 ? Lx/Ax : 0.0);
        Ry = (Ay>0.0 ? Ly/Ay : 0.0);
      }
      if (Rx > Ry) x_gt_y = true;
    });
    return x_gt_y ? 0 : 1;                         // 1 iff "X FSD Y"
  }
  
  // SSD/TSD gates
  if (!(X.mn >= Y.mn) || (Y.mean > X.mean)) return 0;
  if (identical_samples(X, Y)) return 0;           // identical series -> 0
  
  if (degree==2){ // SSD: compare LPM degree 1
    bool x_gt_y = false;
    for_each_threshold(X, Y, [&](double t, int kx, int ky){
      double Lx, Ux, Ly, Uy; (void)Ux; (void)Uy;   // not used beyond calc
      lpm_upm_deg1(X, kx, t, Lx, Ux);
      lpm_upm_deg1(Y, ky, t, Ly, Uy);
      if (Lx > Ly) x_gt_y = true;
    });
    return x_gt_y ? 0 : 1;                         // 1 iff "X SSD Y"
  }
  
  // TSD: compare LPM degree 2
  bool x_gt_y = false;
  for_each_threshold(X, Y, [&](double t, int kx, int ky){
    double Lx2 = lpm_deg2(X, kx, t);
    double Ly2 = lpm_deg2(Y, ky, t);
    if (Lx2 > Ly2) x_gt_y = true;
  });
  return x_gt_y ? 0 : 1;                           // 1 iff "X TSD Y"
}

// =====================================================================
// Parallel dominance matrix (rows i)
// =====================================================================
struct DomWorker : public Worker {
  const std::vector<ColPre>& cols;
  const int degree;
  const bool discrete;
  RMatrix<int> D;
  DomWorker(const std::vector<ColPre>& cols_, int degree_, bool discrete_, IntegerMatrix& D_)
    : cols(cols_), degree(degree_), discrete(discrete_), D(D_) {}
  void operator()(std::size_t begin, std::size_t end) override {
    const int n = D.nrow();
    for (std::size_t i = begin; i < end; ++i) {
      const int ii = static_cast<int>(i);  // avoid signed/unsigned compare and index mismatch
      for (int j = 0; j < n; ++j) {
        D(ii, j) = (ii == j ? 0 : sd_dom_pair(cols[ii], cols[j], degree, discrete));
      }
    }
  }
};

// =====================================================================
// Export: dominance matrix in parallel (prefix-sum version)
// =====================================================================
// [[Rcpp::export]]
IntegerMatrix sd_dom_matrix_prefix_parallel(const NumericMatrix& X, int degree, std::string type="discrete"){
  if (!(degree==1 || degree==2 || degree==3)) stop("degree must be 1, 2, or 3");
  for (R_xlen_t k=0; k<X.size(); ++k) if (NumericVector::is_na(X[k]))
    stop("You have some missing values, please address.");
  
  std::transform(type.begin(), type.end(), type.begin(), ::tolower);
  bool discrete = true;
  if (degree==1){
    if (!(type=="discrete" || type=="continuous"))
      warning("type needs to be either discrete or continuous");
    discrete = (type != "continuous");
  }
  
  int n = X.ncol();
  std::vector<ColPre> cols; cols.reserve(n);
  for (int j=0;j<n;++j) cols.push_back(precompute_col(X, j));
  
  IntegerMatrix D(n, n);
  DomWorker w(cols, degree, discrete, D);
  parallelFor(0, n, w);
  return D;
}


// small inline helper: repeated multiplication for integer exponents
static inline double repeatMultiplication(double value, int n) {
  double result = 1.0;
  for (int i = 0; i < n; ++i) result *= value;
  return result;
}

// =====================================================================
// Export: Efficient set (LPM(1, tmax, ·) order + single scan)
// =====================================================================
// [[Rcpp::export]]
CharacterVector NNS_SD_efficient_set_parallel_cpp(NumericMatrix X, int degree, std::string type="discrete", bool status=true){
  int n = X.ncol(); if (!n) return CharacterVector();
  
  CharacterVector cn = colnames(X);
  if (cn.size()!=n){ cn = CharacterVector(n); for(int j=0;j<n;++j) cn[j] = "X_"+std::to_string(j+1); colnames(X)=cn; }
  
  // global max for ordering key
  double tmax = R_NegInf; for (R_xlen_t k=0;k<X.size();++k) if (!NumericVector::is_na(X[k]) && X[k]>tmax) tmax = X[k];
  
  // precompute columns
  std::vector<ColPre> cols; cols.reserve(n);
  for (int j=0;j<n;++j) cols.push_back(precompute_col(X, j));
  
  // ===== order by LPM(degree, tmax, ·) =====
  std::vector<double> lpm_vals(n, 0.0);
  std::vector<int> non_na_counts(n, 0);
  
  // compute LPM_r for each column (simple O(nrows * ncols) pass)
  for (int j = 0; j < n; ++j) {
    double sum = 0.0;
    int cnt = 0;
    for (int i = 0; i < X.nrow(); ++i) {
      double xv = X(i, j);
      if (!NumericVector::is_na(xv)) {
        double diff = tmax - xv;
        if (diff > 0.0) {
          // use integer repeated multiplication (faster/more accurate than std::pow for integer exponents)
          sum += repeatMultiplication(diff, degree);
        }
        cnt++;
      }
    }
    if (cnt > 0) lpm_vals[j] = sum / (double) cnt;
    else lpm_vals[j] = R_PosInf;   // push all-NA columns to the end
    non_na_counts[j] = cnt;
  }
  
  // Build index vector and sort by lpm_vals (ascending: lower LPM earlier)
  std::vector<int> ord(n);
  for (int j = 0; j < n; ++j) ord[j] = j;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) {
              if (lpm_vals[a] == lpm_vals[b]) return a < b; // stable tie-break by index
              return lpm_vals[a] < lpm_vals[b];
            });
  
  // build dominance matrix in the sorted order
  NumericMatrix Xo(X.nrow(), n); CharacterVector names_sorted(n);
  for (int k=0;k<n;++k){ Xo(_,k)=X(_,ord[k]); names_sorted[k]=cn[ord[k]]; }
  
  IntegerMatrix D = sd_dom_matrix_prefix_parallel(Xo, degree, type);
  
  // single pass to keep maximal elements
  LogicalVector keep(n);
  for (int k=0;k<n;++k){
    if (status){ if (k<n-1) Rcpp::Rcout << "Checking " << (k+1) << " of " << (n-1) << "\r";
    if (k==n-1) Rcpp::Rcout << std::string(40,' ') << "\n"; }
    bool dominated=false;
    for (int i=0;i<k;++i) if (keep[i] && D(i,k)==1){ dominated=true; break; }
    keep[k] = !dominated;
  }
  std::vector<std::string> out; out.reserve(n);
  for (int k=0;k<n;++k) if (keep[k]) out.push_back( as<std::string>(names_sorted[k]) );
  return wrap(out);
}

// =====================================================================
// Minimal one-pair exports (for R wrappers NNS.*.uni)
// =====================================================================
// [[Rcpp::export]]
int NNS_FSD_uni_cpp(const NumericVector& x, const NumericVector& y, std::string type = "discrete"){
  if (is_true(any(is_na(x))) || is_true(any(is_na(y))))
    stop("You have some missing values, please address.");
  std::transform(type.begin(), type.end(), type.begin(), ::tolower);
  bool discrete = (type != "continuous");
  ColPre X = precompute_vec(x), Y = precompute_vec(y);
  return sd_dom_pair(X, Y, 1, discrete);
}

// [[Rcpp::export]]
int NNS_SSD_uni_cpp(const NumericVector& x, const NumericVector& y){
  if (is_true(any(is_na(x))) || is_true(any(is_na(y))))
    stop("You have some missing values, please address.");
  ColPre X = precompute_vec(x), Y = precompute_vec(y);
  return sd_dom_pair(X, Y, 2, true);
}

// [[Rcpp::export]]
int NNS_TSD_uni_cpp(const NumericVector& x, const NumericVector& y){
  if (is_true(any(is_na(x))) || is_true(any(is_na(y))))
    stop("You have some missing values, please address.");
  ColPre X = precompute_vec(x), Y = precompute_vec(y);
  return sd_dom_pair(X, Y, 3, true);
}