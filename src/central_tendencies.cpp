// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ---------- helpers ----------

static inline double frac_part(double x) {
  return x - std::floor(x);
}

static inline double mean_vec(const std::vector<double>& v) {
  if (v.empty()) return NA_REAL;
  long double s = 0.0L;
  for (double x : v) s += x;
  return static_cast<double>(s / v.size());
}

static inline double nearest_int_half_up(double x) {
  double f = std::floor(x);
  return ( (x - f) < 0.5 ) ? f : std::ceil(x);
}

// Given a sorted vector xs, reproduce the q1, q2, q3 *exactly* as in the R code.
// - If length is even: q1 = xs[l*.25], q2 = xs[l*.5], q3 = xs[l*.75]  (1-based indexing)
// - If length is odd : q1 and q3 via linear interpolation at l*p, q2 = average of floor/ceil
static void quartiles_like_R_code(const std::vector<double>& xs, double& q1, double& q2, double& q3) {
  const int l = static_cast<int>(xs.size());
  const double l25 = l * 0.25;
  const double l50 = l * 0.50;
  const double l75 = l * 0.75;
  
  if (l % 2 == 0) {
    // 1-based positions in R -> convert to 0-based for C++
    int i25 = std::max(1, (int)std::floor(l25)) - 1;
    int i50 = std::max(1, (int)std::floor(l50)) - 1;
    int i75 = std::max(1, (int)std::floor(l75)) - 1;
    q1 = xs[i25];
    q2 = xs[i50];
    q3 = xs[i75];
  } else {
    // q1 interpolated
    int f25 = (int)std::floor(l25);
    int c25 = (int)std::ceil(l25);
    f25 = std::max(1, f25); // guard
    c25 = std::max(1, c25);
    f25 = std::min(f25, l);
    c25 = std::min(c25, l);
    double w25 = frac_part(l25);
    q1 = xs[f25 - 1] + w25 * (xs[c25 - 1] - xs[f25 - 1]);
    
    // q2 average of floor and ceil positions
    int f50 = (int)std::floor(l50);
    int c50 = (int)std::ceil(l50);
    f50 = std::max(1, f50);
    c50 = std::max(1, c50);
    f50 = std::min(f50, l);
    c50 = std::min(c50, l);
    q2 = 0.5 * (xs[f50 - 1] + xs[c50 - 1]);
    
    // q3 interpolated
    int f75 = (int)std::floor(l75);
    int c75 = (int)std::ceil(l75);
    f75 = std::max(1, f75);
    c75 = std::max(1, c75);
    f75 = std::min(f75, l);
    c75 = std::min(c75, l);
    double w75 = frac_part(l75);
    q3 = xs[f75 - 1] + w75 * (xs[c75 - 1] - xs[f75 - 1]);
  }
}

// Minimal replacement for NNS_bin used by mode/gravity:
// Given sorted x, fixed bin width, origin, return counts and the width.
// We align with the R usage where z_names <- seq(x1, xl, width) and length(counts) matches z_names length.
//
// We assign each x to idx = floor((x - origin) / width), clipped to [0, nbins-1].
//
static void simple_bin_counts(const std::vector<double>& xs,
                              double width, double origin,
                              std::vector<double>& bin_names,
                              std::vector<int>& counts) {
  const int l = (int)xs.size();
  if (l == 0) { bin_names.clear(); counts.clear(); return; }
  

  const double xmax = xs.back();
  // Number of bins so that last bin_name <= xmax and bin_names[k] = origin + k*width
  int nbins = (int)std::floor( (xmax - origin) / width + 1e-12 ) + 1;
  if (nbins < 1) nbins = 1;
  
  bin_names.resize(nbins);
  for (int k = 0; k < nbins; ++k) bin_names[k] = origin + k * width;
  
  counts.assign(nbins, 0);
  for (double v : xs) {
    int idx = (int)std::floor((v - origin) / width);
    if (idx < 0) idx = 0;
    if (idx >= nbins) idx = nbins - 1;
    counts[idx] += 1;
  }
}

// ---------- NNS.mode ----------

// [[Rcpp::export]]
SEXP NNS_mode_cpp(SEXP xSEXP, bool discrete = false, bool multi = true) {
  NumericVector xR(xSEXP);
  std::vector<double> x(xR.begin(), xR.end());
  
  // Coerce to numeric & handle trivial cases
  std::vector<double> xnum;
  xnum.reserve(x.size());
  for (double v : x) if (R_finite(v)) xnum.push_back((double)v);
  
  const int l = (int)xnum.size();
  if (l == 0) return Rf_ScalarReal(NA_REAL);
  if (l <= 3) {
    // median(x)
    std::vector<double> tmp = xnum;
    std::sort(tmp.begin(), tmp.end());
    double med;
    if (l % 2 == 1) med = tmp[l/2];
    else med = 0.5*(tmp[l/2 - 1] + tmp[l/2]);
    if (discrete) return Rf_ScalarReal( nearest_int_half_up(med) );
    return Rf_ScalarReal(med);
  }
  
  // All-equal?
  bool all_eq = true;
  for (int i = 1; i < l; ++i) if (xnum[i] != xnum[0]) { all_eq = false; break; }
  if (all_eq) return Rf_ScalarReal(xnum[0]);
  
  // Sort
  std::sort(xnum.begin(), xnum.end());
  double range = std::fabs(xnum.back() - xnum.front());
  if (range == 0.0) return Rf_ScalarReal(xnum.front());
  
  // Quartiles and width
  double q1, q2, q3;
  quartiles_like_R_code(xnum, q1, q2, q3);
  
  double width = (q3 - q1) * std::pow((double)l, -0.5);
  if (!(width > 0.0) || !R_finite(width)) {
    width = range / 128.0;
  }
  
  // Bin; fallback if degenerate
  std::vector<double> z_names;
  std::vector<int> counts;
  if (width <= 0.0 || !R_finite(width)) width = range / 128.0;
  simple_bin_counts(xnum, width, xnum.front(), z_names, counts);
  const int lz = (int)counts.size();
  
  // Find max count
  int maxc = 0;
  for (int c : counts) if (c > maxc) maxc = c;
  int ties = 0;
  for (int c : counts) if (c == maxc) ++ties;
  
  if (ties > 1 && multi) {
    // Return all bin "names" at max count
    NumericVector out;
    for (int i = 0; i < lz; ++i) if (counts[i] == maxc) out.push_back(z_names[i]);
    if (discrete) {
      for (int i = 0; i < out.size(); ++i) out[i] = nearest_int_half_up(out[i]);
    }
    return out;
  } else {
    // Weighted center around winning bin ±1
    int zc = 0;
    for (int i = 0; i < lz; ++i) if (counts[i] == maxc) { zc = i; break; }
    int lo = std::max(0, zc - 1);
    int hi = std::min(lz - 1, zc + 1);
    long double num = 0.0L, den = 0.0L;
    for (int i = lo; i <= hi; ++i) { num += (long double)z_names[i] * (long double)counts[i]; den += (long double)counts[i]; }
    double finalv = (den > 0.0L) ? (double)(num / den) : z_names[zc];
    if (discrete) finalv = nearest_int_half_up(finalv);
    if (multi) return Rf_ScalarReal(finalv);
    // The original code returns mean(final) if multi == FALSE, but final is scalar; keep scalar.
    return Rf_ScalarReal(finalv);
  }
}

// ---------- NNS.gravity ----------

// [[Rcpp::export]]
SEXP NNS_gravity_cpp(SEXP xSEXP, bool discrete = false) {
  NumericVector xR(xSEXP);
  std::vector<double> x;
  x.reserve(xR.size());
  for (double v : xR) if (R_finite(v)) x.push_back(v);
  
  const int l = (int)x.size();
  if (l == 0) return Rf_ScalarReal(NA_REAL);
  if (l <= 3) {
    // median(x)
    std::vector<double> t = x;
    std::sort(t.begin(), t.end());
    double med = (l % 2 ? t[l/2] : 0.5*(t[l/2 - 1] + t[l/2]));
    if (discrete) return Rf_ScalarReal( nearest_int_half_up(med) );
    return Rf_ScalarReal(med);
  }
  
  bool all_eq = true;
  for (int i = 1; i < l; ++i) if (x[i] != x[0]) { all_eq = false; break; }
  if (all_eq) return Rf_ScalarReal(x[0]);
  
  std::sort(x.begin(), x.end());
  double range = std::fabs(x.back() - x.front());
  if (range == 0.0) return Rf_ScalarReal(x.front());
  
  double q1, q2, q3;
  quartiles_like_R_code(x, q1, q2, q3);
  
  double width = (q3 - q1) * std::pow((double)l, -0.5);
  if (!(width > 0.0) || !R_finite(width)) width = range / 128.0;
  
  std::vector<double> z_names;
  std::vector<int> counts;
  simple_bin_counts(x, width, x.front(), z_names, counts);
  const int lz = (int)counts.size();
  
  // If unique max, use neighborhood; else use all bins
  int maxc = 0;
  for (int c : counts) if (c > maxc) maxc = c;
  int ties = 0;
  for (int c : counts) if (c == maxc) ++ties;
  
  int lo = 0, hi = lz - 1;
  if (ties == 1) {
    int zc = 0; for (int i = 0; i < lz; ++i) if (counts[i] == maxc) { zc = i; break; }
    lo = std::max(0, zc - 1);
    hi = std::min(lz - 1, zc + 1);
  }
  
  long double num = 0.0L, den = 0.0L;
  for (int i = lo; i <= hi; ++i) { num += (long double)z_names[i] * (long double)counts[i]; den += (long double)counts[i]; }
  double m = (den > 0.0L) ? (double)(num / den) : z_names[ (lo+hi)/2 ];
  
  double mu = mean_vec(x);
  double mid = 0.25 * ( q2 + m + mu + 0.5*(q1 + q3) );
  
  double out = R_finite(mid) ? mid : q2;
  if (discrete) out = nearest_int_half_up(out);
  return Rf_ScalarReal(out);
}

// ---------- NNS.rescale ----------

// [[Rcpp::export]]
NumericVector NNS_rescale_cpp(SEXP xSEXP, double a, double b,
                              std::string method = "minmax",
                              Rcpp::Nullable<double> T_ = R_NilValue,
                              std::string type = "Terminal") {
  NumericVector xR(xSEXP);
  int n = xR.size();
  NumericVector out(n);
  
  std::transform(method.begin(), method.end(), method.begin(), ::tolower);
  std::transform(type.begin(), type.end(), type.begin(), ::tolower);
  
  if (method == "minmax") {
    double xmin = R_PosInf, xmax = R_NegInf;
    for (int i = 0; i < n; ++i) {
      if (R_finite(xR[i])) {
        if (xR[i] < xmin) xmin = xR[i];
        if (xR[i] > xmax) xmax = xR[i];
      }
    }
    if (!R_finite(xmin) || !R_finite(xmax) || xmax == xmin) {
      Rcpp::warning("All x identical: returning midpoint values");
      for (int i = 0; i < n; ++i) out[i] = (a + b) / 2.0;
      return out;
    }
    for (int i = 0; i < n; ++i) {
      out[i] = a + (b - a) * ( (xR[i] - xmin) / (xmax - xmin) );
    }
    return out;
  }
  
  if (method == "riskneutral") {
    if (T_.isNull()) stop("T (time to maturity) must be provided for riskneutral method");
    double T = Rcpp::as<double>(T_);
    if (!(a > 0.0)) stop("S_0 (a) must be positive for riskneutral method");
    double S0 = a;
    double r = b;
    
    // Compute scaling theta so that mean(out) matches target
    long double s = 0.0L; int cnt = 0;
    for (int i = 0; i < n; ++i) if (R_finite(xR[i])) { s += xR[i]; ++cnt; }
    double mx = (cnt > 0) ? (double)(s / cnt) : NA_REAL;
    
    if (!R_finite(mx) || mx <= 0.0)
      stop("Mean(x) must be positive/finite for riskneutral scaling");
    
    double target = (type == "discounted") ? S0 : (S0 * std::exp(r * T));
    double theta = std::log(target / mx);
    
    for (int i = 0; i < n; ++i) out[i] = xR[i] * std::exp(theta);
    return out;
  }
  
  stop("Invalid method: use 'minmax' or 'riskneutral'");
  return out; // never reached
}
