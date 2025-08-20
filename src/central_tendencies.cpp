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


// ---------- NNS.mode ----------

// --- Triangular smoothing helper: 7-tap [1,2,3,4,3,2,1] with mirrored edges ---
static void smooth_counts_tri7(const std::vector<int>& counts, std::vector<double>& smooth) {
  static const int w[7] = {1,2,3,4,3,2,1};
  static const int Wsum = 16; // 1+2+3+4+3+2+1
  const int n = (int)counts.size();
  smooth.assign(n, 0.0);
  if (n == 0) return;
  
  // Mirror at edges (symmetric extension)
  auto at = [&](int idx)->int{
    if (idx < 0)   return counts[-idx];            // reflect: -1 -> 1, -2 -> 2, ...
    if (idx >= n)  return counts[2*n - 2 - idx];   // reflect: n -> n-2, n+1 -> n-3, ...
    return counts[idx];
  };
  
  for (int i = 0; i < n; ++i) {
    int acc = 0;
    acc += w[0]*at(i-3); acc += w[1]*at(i-2); acc += w[2]*at(i-1);
    acc += w[3]*at(i  );
    acc += w[4]*at(i+1); acc += w[5]*at(i+2); acc += w[6]*at(i+3);
    smooth[i] = (double)acc / (double)Wsum;
  }
}

// [[Rcpp::export]]
SEXP NNS_mode_cpp(SEXP xSEXP, bool discrete = false, bool multi = true) {
  NumericVector xR(xSEXP);
  std::vector<double> x(xR.begin(), xR.end());
  
  // Coerce to numeric & drop non-finite
  std::vector<double> xnum; xnum.reserve(x.size());
  for (double v : x) if (R_finite(v)) xnum.push_back((double)v);
  
  const int l = (int)xnum.size();
  if (l == 0) return Rf_ScalarReal(NA_REAL);
  
  // ====================== DISCRETE PATH ======================
  if (discrete) {
    if (l <= 3) {
      // For tiny samples, integerized median
      std::vector<double> tmp = xnum; std::sort(tmp.begin(), tmp.end());
      double med = (l % 2 == 1) ? tmp[l/2] : 0.5*(tmp[l/2 - 1] + tmp[l/2]);
      return Rf_ScalarReal(nearest_int_half_up(med));
    }
    
    // Integerize and count exact frequencies
    std::unordered_map<int,int> freq; freq.reserve(l * 2u);
    for (double v : xnum) ++freq[ nearest_int_half_up(v) ];
    
    int maxf = 0; for (auto &kv : freq) if (kv.second > maxf) maxf = kv.second;
    
    std::vector<int> modes_int;
    for (auto &kv : freq) if (kv.second == maxf) modes_int.push_back(kv.first);
    std::sort(modes_int.begin(), modes_int.end());
    
    if (multi) {
      NumericVector out((int)modes_int.size());
      for (int i = 0; i < (int)modes_int.size(); ++i) out[i] = (double)modes_int[i];
      return out;            // e.g., 2 3 4 for c(1,2,2,3,3,4,4,5)
    } else {
      // Return the arithmetic mean of all tied modes
      long double sum = 0.0L;
      for (int m : modes_int) sum += (long double)m;
      double mean_modes = (modes_int.empty() ? NA_REAL
                             : (double)(sum / (long double)modes_int.size()));
      return Rf_ScalarReal(mean_modes);
    }
  }
  
  // ====================== CONTINUOUS PATH ======================
  if (l <= 3) {
    std::vector<double> tmp = xnum; std::sort(tmp.begin(), tmp.end());
    double med = (l % 2 == 1) ? tmp[l/2] : 0.5*(tmp[l/2 - 1] + tmp[l/2]);
    return Rf_ScalarReal(med);
  }
  
  // All-equal?
  bool all_eq = true;
  for (int i = 1; i < l; ++i) if (xnum[i] != xnum[0]) { all_eq = false; break; }
  if (all_eq) return Rf_ScalarReal(xnum[0]);
  
  // Sort & basic stats
  std::sort(xnum.begin(), xnum.end());
  double range = std::fabs(xnum.back() - xnum.front());
  if (range == 0.0) return Rf_ScalarReal(xnum.front());
  
  // Quartiles & default bin width
  double q1, q2, q3;
  quartiles_like_R_code(xnum, q1, q2, q3);
  double width = (q3 - q1) * std::pow((double)l, -0.5);
  if (!(width > 0.0) || !R_finite(width)) width = range / 128.0;
  
  // Histogram
  std::vector<double> z_names;   // representative x for each bin (center/name)
  std::vector<int> counts;       // histogram counts
  if (width <= 0.0 || !R_finite(width)) width = range / 128.0;
  simple_bin_counts(xnum, width, xnum.front(), z_names, counts);
  const int lz = (int)counts.size();
  if (lz == 0) return Rf_ScalarReal(NA_REAL);
  
  // For fallback paths
  int maxc = 0; for (int c : counts) if (c > maxc) maxc = c;
  
  // ----- Peak detection on SMOOTHED counts (edge-aware 1..3 & concavity) -----
  std::vector<double> cs; smooth_counts_tri7(counts, cs);
  
  // Optional margin above side maxima (in smoothed counts)
  const double MARGIN = 0.0;
  
  std::vector<int> peak_idx; peak_idx.reserve(lz);
  // Require full neighborhoods for offsets 1..3
  for (int i = 3; i <= lz - 4; ++i) {
    double ci = cs[i];
    if (ci <= 0.0) continue;
    
    // Max of neighbors at offsets 1..3 on each side (smoothed series)
    double Ls = std::max(std::max(cs[i-1], cs[i-2]), cs[i-3]);
    double Rs = std::max(std::max(cs[i+1], cs[i+2]), cs[i+3]);
    if (!(ci > Ls + MARGIN && ci > Rs + MARGIN)) continue;
    
    // Negative curvature gate (concavity)
    double curv = cs[i-1] - 2.0*cs[i] + cs[i+1];
    if (!(curv < 0.0)) continue;
    
    peak_idx.push_back(i);
  }
  
  // Non-maximum suppression on smoothed heights: keep peaks >= 4 bins apart
  if (!peak_idx.empty()) {
    std::sort(peak_idx.begin(), peak_idx.end(),
              [&](int a, int b){ return cs[a] > cs[b]; });
    std::vector<int> kept;
    for (int idx : peak_idx) {
      bool too_close = false;
      for (int jdx : kept) if (std::abs(idx - jdx) <= 3) { too_close = true; break; }
      if (!too_close) kept.push_back(idx);
    }
    
    if (!kept.empty()) {
      // Per-peak weighted center over ±3 bins using ORIGINAL counts
      NumericVector centers((int)kept.size());
      for (int t = 0; t < (int)kept.size(); ++t) {
        int zc = kept[t];
        int lo = std::max(0, zc - 3);
        int hi = std::min(lz - 1, zc + 3);
        long double num = 0.0L, den = 0.0L;
        for (int j = lo; j <= hi; ++j) {
          if (std::abs(j - zc) <= 3) {
            num += (long double)z_names[j] * (long double)counts[j];
            den += (long double)counts[j];
          }
        }
        double m = (den > 0.0L) ? (double)(num / den) : z_names[zc];
        centers[t] = m;
      }
      
      if (multi) {
        NumericVector out = clone(centers);
        std::sort(out.begin(), out.end());
        return out;
      } else {
        // GLOBAL-HEIGHT RULE: choose the kept peak with largest smoothed height cs[i]
        int best_t = 0;
        for (int t = 1; t < (int)kept.size(); ++t) {
          if (cs[kept[t]] > cs[kept[best_t]]) best_t = t;
        }
        return Rf_ScalarReal(centers[best_t]);
      }
    }
  }
  
  // Fallback: if multiple global-max bins exist
  int ties = 0; for (int c : counts) if (c == maxc) ++ties;
  if (ties > 1) {
    if (multi) {
      NumericVector out(ties);
      int pos = 0;
      for (int i = 0; i < lz; ++i) if (counts[i] == maxc) out[pos++] = z_names[i];
      std::sort(out.begin(), out.end());
      return out;
    } else {
      // Mean of those bin centers when multi==false
      long double sum = 0.0L;
      int pos = 0;
      for (int i = 0; i < lz; ++i) if (counts[i] == maxc) { sum += (long double)z_names[i]; ++pos; }
      double mean_modes = (pos > 0 ? (double)(sum / (long double)pos) : NA_REAL);
      return Rf_ScalarReal(mean_modes);
    }
  }
  
  // Final fallback: single winning bin -> weighted center around ±1
  int zc = 0; for (int i = 0; i < lz; ++i) if (counts[i] == maxc) { zc = i; break; }
  {
    int lo = std::max(0, zc - 1);
    int hi = std::min(lz - 1, zc + 1);
    long double num = 0.0L, den = 0.0L;
    for (int j = lo; j <= hi; ++j) {
      num += (long double)z_names[j] * (long double)counts[j];
      den += (long double)counts[j];
    }
    double finalv = (den > 0.0L) ? (double)(num / den) : z_names[zc];
    return Rf_ScalarReal(finalv);
  }
}
