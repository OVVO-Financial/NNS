// SD_prefix_refactor.cpp
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <string>
#include <cfloat>
#include <cmath>
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

inline void neumaier_add(double& sum, double& err, double v){
  double t = sum + v;
  err += (std::fabs(sum) >= std::fabs(v)) ? (sum - t) + v : (v - t) + sum;
  sum = t;
}

static void finish_col(ColPre& c, bool sorted = false){
  if (c.m == 0) stop("Each sample needs at least one observation.");
  if (!sorted) std::sort(c.vals.begin(), c.vals.end());
  c.P1.assign(c.m+1, 0.0);
  c.P2.assign(c.m+1, 0.0);
  // Compensated sums: a plain running sum drifts by up to ~m*eps*|x|,
  // which for large samples exceeds the tie tolerance in sd_dom_pair.
  double s1 = 0.0, e1 = 0.0, s2 = 0.0, e2 = 0.0;
  for(int k=1;k<=c.m;++k){
    double v = c.vals[k-1];
    neumaier_add(s1, e1, v);
    neumaier_add(s2, e2, v*v);
    c.P1[k] = s1 + e1;
    c.P2[k] = s2 + e2;
  }
  c.S1   = c.P1[c.m];
  c.S2   = c.P2[c.m];
  c.mn   = c.vals.front();
  c.mean = c.S1 / double(c.m);
}

// Same sample shifted by s; subtraction keeps the sort order.
static ColPre shifted_col(const ColPre& c, double s){
  ColPre o; o.m = c.m; o.vals.resize(c.m);
  for(int i=0;i<c.m;++i) o.vals[i] = c.vals[i] - s;
  finish_col(o, true);
  return o;
}

static ColPre precompute_col(const NumericMatrix& X, int j){
  ColPre c; c.m = X.nrow();
  c.vals.resize(c.m);
  for(int i=0;i<c.m;++i) c.vals[i] = X(i,j);
  finish_col(c);
  return c;
}

static ColPre precompute_vec(const NumericVector& x){
  ColPre c; c.m = x.size();
  c.vals.assign(x.begin(), x.end());
  finish_col(c);
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
// Each sample is walked with its own length (the samples may differ in size).
template <class F>
inline void for_each_threshold(const ColPre& a, const ColPre& b, F f){
  int ia=0, ib=0, ma=a.m, mb=b.m;
  while(ia<ma || ib<mb){
    double next_a = (ia<ma ? a.vals[ia] : R_PosInf);
    double next_b = (ib<mb ? b.vals[ib] : R_PosInf);
    double t = (next_a < next_b ? next_a : next_b);
    while(ia<ma && a.vals[ia] <= t) ++ia; // k_a = ia
    while(ib<mb && b.vals[ib] <= t) ++ib; // k_b = ib
    // Next knot: the counts (ia, ib) are constant on (t, t_next).
    double t_next = std::min(ia<ma ? a.vals[ia] : R_PosInf,
                             ib<mb ? b.vals[ib] : R_PosInf);
    f(t, ia, ib, t_next);
  }
}

// Continuous FSD ratio L1/(L1+U1) on an interval with fixed count k:
// L1 = alpha*t + beta, L1+U1 = gamma*t + delta (both linear in t).
inline void ratio_deg1_coefs(const ColPre& c, int k, double& alpha, double& beta,
                             double& gamma, double& delta){
  double m = double(c.m);
  alpha = k / m;
  beta  = -c.P1[k] / m;
  gamma = (2.0*k - c.m) / m;
  delta = (c.S1 - 2.0*c.P1[k]) / m;
}

// =====================================================================
// Pairwise dominance via prefix sums (O(m) per pair)
// degree: 1=FSD, 2=SSD, 3=TSD. 'discrete' only matters for FSD.
// Returns 1 iff X dominates Y, else 0.
// =====================================================================
inline int sd_dom_pair(const ColPre& X, const ColPre& Y, int degree, bool discrete, double level = 0.0);

// Dominance is translation invariant but the tolerances below scale with
// |value|^degree, so compare the pair centred on its joint midrange. The data
// still carry rounding at their own level |c|, which sd_dom_pair adds back.
inline int sd_dom_pair_centred(const ColPre& X, const ColPre& Y, int degree, bool discrete){
  if (degree==1 && discrete) return sd_dom_pair(X, Y, degree, discrete);
  if (!(X.mn >= Y.mn) || identical_samples(X, Y)) return 0;
  const double lo = std::min(X.vals.front(), Y.vals.front());
  const double hi = std::max(X.vals.back(), Y.vals.back());
  // A mean gap this large fails the centred gate too; skip the copies.
  if (degree >= 2 && Y.mean - X.mean > 2.0 * 64.0 * DBL_EPSILON * std::max(std::fabs(lo), std::fabs(hi))) return 0;
  const double c = 0.5*lo + 0.5*hi;
  if (!R_FINITE(c)) return sd_dom_pair(X, Y, degree, discrete);  // an infinite value can't be shifted
  return sd_dom_pair(shifted_col(X, c), shifted_col(Y, c), degree, discrete, std::fabs(c));
}

inline int sd_dom_pair(const ColPre& X, const ColPre& Y, int degree, bool discrete, double level){
  // Prefix-sum rounding makes curves that are equal in exact arithmetic differ
  // by about eps * S^degree (S = largest |value|), so differences up to
  // 64 * DBL_EPSILON * S^degree count as ties. Discrete FSD uses no tolerance.
  const double S = std::max(std::max(std::fabs(X.vals.front()), std::fabs(X.vals.back())),
                            std::max(std::fabs(Y.vals.front()), std::fabs(Y.vals.back())));
  const double eps64 = 64.0 * DBL_EPSILON;
  // Error bounds: summation rounding plus the inputs' own rounding at |level|
  // (a few ulps per value). For LPM2 at t that rounding moves each term by about
  // 2 (t - x) din, so its bound uses the LPM1 at t.
  const double din = 4.0 * DBL_EPSILON * level;
  const double E1 = eps64 * S + din;
  auto E2 = [&](int kx, int ky, double t){
    double Lx, Ux, Ly, Uy;
    lpm_upm_deg1(X, kx, t, Lx, Ux);
    lpm_upm_deg1(Y, ky, t, Ly, Uy);
    return eps64 * S * S + 2.0 * din * (std::fabs(Lx) + std::fabs(Ly) + din);
  };
  if (degree==1){ // FSD
    if (!(X.mn >= Y.mn)) return 0;                 // FSD gate
    if (identical_samples(X, Y)) return 0;         // identical series -> 0 (matches R's identical(LPM_x, LPM_y))
    
    bool x_gt_y = false, curves_differ = false;
    int deg = (discrete ? 0 : 1); // discrete->0, continuous->1
    for_each_threshold(X, Y, [&](double t, int kx, int ky, double t_next){
      double Rx, Ry, tol = 0.0;
      if (deg==0){
        // L0/(L0+U0) == ECDF
        Rx = double(kx)/double(X.m);
        Ry = double(ky)/double(Y.m);
      } else {
        double Lx, Ux, Ly, Uy;
        lpm_upm_deg1(X, kx, t, Lx, Ux);
        lpm_upm_deg1(Y, ky, t, Ly, Uy);
        double Ax = Lx+Ux, Ay = Ly+Uy;
        // 0/0 (every observation equals t): the degree-0 ratio k/m, as in LPM.ratio.
        Rx = (Ax>0.0 ? Lx/Ax : double(kx)/double(X.m));
        Ry = (Ay>0.0 ? Ly/Ay : double(ky)/double(Y.m));
        double Amin = std::min(Ax, Ay);
        tol = (Amin > 0.0) ? std::max(eps64, E1 / Amin) : eps64;
      }
      if (Rx - Ry > tol) x_gt_y = true;
      if (std::fabs(Rx - Ry) > tol) curves_differ = true;
      
      // Continuous ratios are linear-fractional between knots, so their
      // difference can peak strictly inside (t, t_next). Its stationary
      // points solve Ax / Dx(t)^2 = Ay / Dy(t)^2 with Dx, Dy > 0, i.e.
      // Dy(t) = sqrt(Ay / Ax) * Dx(t) (only when Ax, Ay have the same sign).
      if (deg==1 && R_FINITE(t_next) && !x_gt_y){
        double ax, bx, gx, dx, ay, by, gy, dy;
        ratio_deg1_coefs(X, kx, ax, bx, gx, dx);
        ratio_deg1_coefs(Y, ky, ay, by, gy, dy);
        double Ax = ax*dx - bx*gx, Ay = ay*dy - by*gy;
        if (Ax != 0.0 && Ay != 0.0 && ((Ax > 0.0) == (Ay > 0.0))){
          double r = std::sqrt(Ay / Ax);
          double den = gy - r*gx;
          if (den != 0.0){
            double ts = (r*dx - dy) / den;
            if (ts > t && ts < t_next){
              double Dx = gx*ts + dx, Dy = gy*ts + dy;
              if (Dx > 0.0 && Dy > 0.0){
                double Rxs = (ax*ts + bx) / Dx, Rys = (ay*ts + by) / Dy;
                double tol_r = std::max(eps64, E1 / std::min(Dx, Dy));
                if (Rxs - Rys > tol_r) { x_gt_y = true; curves_differ = true; }
              }
            }
          }
        }
      }
      // Ratios equal at both knots can still differ inside; equal at the
      // midpoint too means equal on the whole interval.
      if (deg==1 && R_FINITE(t_next) && !curves_differ){
        double s = 0.5*(t + t_next), Lx, Ux, Ly, Uy;
        lpm_upm_deg1(X, kx, s, Lx, Ux);
        lpm_upm_deg1(Y, ky, s, Ly, Uy);
        double Dx = Lx+Ux, Dy = Ly+Uy;
        double d = Lx/Dx - Ly/Dy;
        double tol_m = std::max(eps64, E1 / std::min(Dx, Dy));
        if (std::fabs(d) > tol_m) curves_differ = true;
        if (d > tol_m) x_gt_y = true;
      }
    });
    // Equal curves: no dominance, as R's !identical(LPM_x_sort, LPM_y_sort).
    if (!curves_differ) return 0;
    return x_gt_y ? 0 : 1;                         // 1 iff "X FSD Y"
  }
  
  // SSD/TSD gates
  if (!(X.mn >= Y.mn) || (Y.mean - X.mean > E1)) return 0;
  if (identical_samples(X, Y)) return 0;           // identical series -> 0
  
  if (degree==2){ // SSD: compare LPM degree 1
    bool x_gt_y = false, curves_differ = false;
    // LPM1 is piecewise linear with knots at the data: checking the knots is exact.
    for_each_threshold(X, Y, [&](double t, int kx, int ky, double){
      double Lx, Ux, Ly, Uy; (void)Ux; (void)Uy;   // not used beyond calc
      lpm_upm_deg1(X, kx, t, Lx, Ux);
      lpm_upm_deg1(Y, ky, t, Ly, Uy);
      if (Lx - Ly > E1) x_gt_y = true;
      if (std::fabs(Lx - Ly) > E1) curves_differ = true;
    });
    if (!curves_differ) return 0;                  // equal curves -> no SSD
    return x_gt_y ? 0 : 1;                         // 1 iff "X SSD Y"
  }
  
  // TSD: compare LPM degree 2
  bool x_gt_y = false, curves_differ = false;
  for_each_threshold(X, Y, [&](double t, int kx, int ky, double t_next){
    double Lx2 = lpm_deg2(X, kx, t);
    double Ly2 = lpm_deg2(Y, ky, t);
    const double tol2 = E2(kx, ky, t);
    if (Lx2 - Ly2 > tol2) x_gt_y = true;
    if (std::fabs(Lx2 - Ly2) > tol2) curves_differ = true;
    // LPM2 is piecewise quadratic: on (t, t_next) the difference is
    // a*s^2 - 2*b*s + c with a = kx/mx - ky/my, b = P1x/mx - P1y/my. When
    // a < 0 it has an interior maximum at s* = b / a, which the knots miss.
    if (R_FINITE(t_next) && !x_gt_y){
      double a = double(kx)/double(X.m) - double(ky)/double(Y.m);
      if (a < 0.0){
        double b = X.P1[kx]/double(X.m) - Y.P1[ky]/double(Y.m);
        double ts = b / a;
        if (ts > t && ts < t_next){
          double Dts = lpm_deg2(X, kx, ts) - lpm_deg2(Y, ky, ts);
          if (Dts > E2(kx, ky, ts)) { x_gt_y = true; curves_differ = true; }
        }
      }
    }
  });
  if (!curves_differ) return 0;                    // equal curves -> no TSD
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
        D(ii, j) = (ii == j ? 0 : sd_dom_pair_centred(cols[ii], cols[j], degree, discrete));
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
  if (cn.size()!=n){ cn = CharacterVector(n); for(int j=0;j<n;++j) cn[j] = "X_"+std::to_string(j+1); }
  
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
  return sd_dom_pair_centred(X, Y, 1, discrete);
}

// [[Rcpp::export]]
int NNS_SSD_uni_cpp(const NumericVector& x, const NumericVector& y){
  if (is_true(any(is_na(x))) || is_true(any(is_na(y))))
    stop("You have some missing values, please address.");
  ColPre X = precompute_vec(x), Y = precompute_vec(y);
  return sd_dom_pair_centred(X, Y, 2, true);
}

// [[Rcpp::export]]
int NNS_TSD_uni_cpp(const NumericVector& x, const NumericVector& y){
  if (is_true(any(is_na(x))) || is_true(any(is_na(y))))
    stop("You have some missing values, please address.");
  ColPre X = precompute_vec(x), Y = precompute_vec(y);
  return sd_dom_pair_centred(X, Y, 3, true);
}