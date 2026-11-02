// NNS_dep.cpp
// C++ implementation of NNS.dep and NNS.dep.matrix.
//
// Exported functions (called from R):
//   NNS_dep_pair_cpp    - bivariate dependence given pre-computed partition labels
// (NNS.dep's matrix form loops over NNS_dep_pair_cpp in R/Dependence.R.)
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cstdint>

using namespace Rcpp;
using namespace RcppParallel;

// ============================================================
// INTERNAL HELPERS
// ============================================================

struct DepResult {
  double correlation;
  double dependence;
};

static inline double gravity_pure_cpp(const std::vector<double>& v) {
  size_t n = v.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return v[0];
  if (n == 2) return (v[0] + v[1]) / 2.0;
  
  double sum = 0.0;
  for (double val : v) sum += val;
  return sum / n;
}

static inline int n_unique(const std::vector<double>& v) {
  std::unordered_map<double, int> seen;
  seen.reserve(v.size());
  for (double d : v) seen[d] = 1;
  return static_cast<int>(seen.size());
}

// ------------------------------------------------------------
// Per-cell copula with the independence floor removed (#70).
//
// Inside a partition cell of m points, around the cell means, the copula
// combines two concordance terms:
//   T0 = |2 p0 - 1|, p0 = share of points in the co-lower or co-upper orthant
//   T1 = |2 p1 - 1|, p1 = degree-1 weighted concordant share
//      = |sum dx dy| / sum |dx dy|
// (the discordant terms are exactly 1 - p0 and 1 - p1, so they add nothing).
// Under independence both terms are positive on average, and averaging the
// cells never cancels them: that gave a Dependence floor that decays only like
// n^(-1/4) (0.32 at n = 1000). Each term is now measured against its own
// expectation under the within-cell permutation null,
//   T* = (T - E0[T]) / (1 - E0[T]),
// which keeps perfect dependence at 1 and centres independence near 0.
// E0[T0] is exact (hypergeometric). E0[T1] is exact by enumeration for m <= 6
// and uses the permutation variance of sum dx dy otherwise. Both cost O(m).
// ------------------------------------------------------------

struct CellCopula {
  double signed_value;  // debiased magnitude times the sign of the cell covariance
  double magnitude;     // debiased magnitude; may be negative below chance
};

// E0|2 p0 - 1| for a cell of m points with a points at or below the x mean and
// b points at or below the y mean, under random pairing: the co-lower count K
// is hypergeometric and the concordant count is m - a - b + 2K.
static double null_T0(int m, int a, int b) {
  const int k_lo = std::max(0, a + b - m);
  const int k_hi = std::min(a, b);
  if (k_lo >= k_hi) {  // a single possible count: the term carries no information
    return 1.0;
  }
  const double md = static_cast<double>(m);
  if (m > 400) {
    // Normal approximation with the exact hypergeometric mean and variance.
    const double ad = a, bd = b;
    const double mean_k = ad * bd / md;
    const double var_k = bd * (ad / md) * ((md - ad) / md) * ((md - bd) / (md - 1.0));
    const double mu = (2.0 * (md - ad - bd + 2.0 * mean_k) - md) / md;
    const double sd = 4.0 * std::sqrt(std::max(var_k, 0.0)) / md;
    if (sd <= 0.0) return std::abs(mu);
    const double z = mu / sd;
    return sd * std::sqrt(2.0 / M_PI) * std::exp(-0.5 * z * z) +
      mu * (1.0 - 2.0 * R::pnorm(-z, 0.0, 1.0, 1, 0));
  }
  double e = 0.0;
  for (int k = k_lo; k <= k_hi; ++k) {
    const double pk = R::dhyper(static_cast<double>(k), static_cast<double>(a),
                                static_cast<double>(m - a), static_cast<double>(b), 0);
    const double conc = static_cast<double>(m - a - b + 2 * k);
    e += pk * std::abs(2.0 * conc / md - 1.0);
  }
  return e;
}

// E0[ |sum dx dy| / sum |dx dy| ] under random pairing of dx and dy.
static double null_T1(const std::vector<double>& dx, const std::vector<double>& dy) {
  const int m = static_cast<int>(dx.size());
  if (m <= 6) {
    std::vector<int> perm(m);
    std::iota(perm.begin(), perm.end(), 0);
    double total = 0.0;
    int count = 0;
    do {
      double s = 0.0, d = 0.0;
      for (int i = 0; i < m; ++i) {
        const double p = dx[i] * dy[perm[i]];
        s += p;
        d += std::abs(p);
      }
      total += d > 0.0 ? std::abs(s) / d : 0.0;
      ++count;
    } while (std::next_permutation(perm.begin(), perm.end()));
    return total / static_cast<double>(count);
  }
  double sx2 = 0.0, sy2 = 0.0, sax = 0.0, say = 0.0;
  for (int i = 0; i < m; ++i) {
    sx2 += dx[i] * dx[i];
    sy2 += dy[i] * dy[i];
    sax += std::abs(dx[i]);
    say += std::abs(dy[i]);
  }
  if (sax <= 0.0 || say <= 0.0) return 1.0;
  // sum dx dy has mean 0 and variance sx2 * sy2 / (m - 1) under permutation;
  // sum |dx dy| has mean sax * say / m.
  const double md = static_cast<double>(m);
  const double e = std::sqrt(2.0 / M_PI) * std::sqrt(sx2 * sy2 / (md - 1.0)) * md / (sax * say);
  return std::min(1.0, e);
}

static inline double debias_term(double t, double e0) {
  if (!(e0 < 1.0 - 1e-12)) return 0.0;
  return (t - e0) / (1.0 - e0);
}

static CellCopula copula_signed(const std::vector<double>& xv,
                                const std::vector<double>& yv) {
  const int n = static_cast<int>(xv.size());
  if (n < 2) return {0.0, 0.0};

  double tx = 0.0, ty = 0.0;
  for (int i = 0; i < n; ++i) { tx += xv[i]; ty += yv[i]; }
  tx /= n; ty /= n;

  std::vector<double> dx(n), dy(n);
  int a = 0, b = 0;
  double conc = 0.0, s = 0.0, d = 0.0;
  for (int i = 0; i < n; ++i) {
    dx[i] = xv[i] - tx;
    dy[i] = yv[i] - ty;
    const bool lo_x = dx[i] <= 0.0, lo_y = dy[i] <= 0.0;
    a += lo_x;
    b += lo_y;
    // Concordant = co-lower (<= target) or co-upper (> target); the degree-0
    // counts partition the cell.
    if (lo_x == lo_y) conc += 1.0;
    const double p = dx[i] * dy[i];
    s += p;
    d += std::abs(p);
  }

  // x constant in the cell (discrete x): the cell has no within-cell
  // variation to test, only whether y is determined there. A constant y is
  // perfect dependence, as before; otherwise the cell is uninformative.
  bool x_const = true, y_const = true;
  for (int i = 0; i < n; ++i) {
    if (dx[i] != 0.0) x_const = false;
    if (dy[i] != 0.0) y_const = false;
  }
  if (x_const) return y_const ? CellCopula{1.0, 1.0} : CellCopula{0.0, 0.0};

  const double slope_sign = (s > 0.0) ? 1.0 : (s < 0.0) ? -1.0 : 0.0;
  const double e0 = null_T0(n, a, b);

  // All concordant or all discordant is perfect dependence, unless the cell is
  // too small for that to be informative (a single attainable count).
  // Compare counts, not count * (1/n): 49 * (1/49) != 1.
  if (e0 < 1.0 - 1e-12) {
    if (conc == static_cast<double>(n)) return {1.0, 1.0};
    if (conc == 0.0) return {-1.0, 1.0};
  }

  const double md = static_cast<double>(n);
  const double t0 = std::abs(2.0 * conc / md - 1.0);
  const double t1 = d > 0.0 ? std::abs(s) / d : 0.0;
  const double e1 = null_T1(dx, dy);

  const double inner = 0.5 * (debias_term(t0, e0) + debias_term(t1, e1));
  const double mag = (inner >= 0.0 ? 1.0 : -1.0) * std::sqrt(std::abs(inner));
  if (slope_sign == 0.0) return {0.0, 0.0};
  return {mag * slope_sign, mag};
}

static double copula_degree0_unsigned(const std::vector<double>& xv,
                                      const std::vector<double>& yv) {
  int n = static_cast<int>(xv.size());
  if (n < 2) return 0.0;
  
  double tx = 0.0, ty = 0.0;
  for (int i = 0; i < n; ++i) { tx += xv[i]; ty += yv[i]; }
  tx /= n; ty /= n;
  
  double d0_cupm = 0.0, d0_clpm = 0.0, dpm_d0_count = 0.0;
  for (int i = 0; i < n; ++i) {
    double dx = xv[i] - tx;
    double dy = yv[i] - ty;
    if (dx > 0.0 && dy > 0.0) d0_cupm += 1.0;
    if (dx <= 0.0 && dy <= 0.0) d0_clpm += 1.0;
    // Discordant = neither co-lower (<= target) nor co-upper (> target), so
    // the degree-0 counts partition the cell.
    if (!((dx <= 0.0 && dy <= 0.0) || (dx > 0.0 && dy > 0.0)))
      dpm_d0_count += 1.0;
  }
  
  double inv_n = 1.0 / static_cast<double>(n);
  double d0_Co = (d0_cupm + d0_clpm) * inv_n;
  double dpm_d0 = dpm_d0_count * inv_n;
  
  constexpr double indep_Co = 0.5;
  // Independence null for the discordant partial moment. DpmCountWorker treats
  // a point as concordant when it is all-below OR all-above the target (both
  // fully-aligned orthants), so under independence P(discordant) = 1 - 2*0.5^n,
  // which is 0.5 for the bivariate (n=2) copula -- not 1 - 0.5^n = 0.75. The
  // old 0.75 anchor left a fixed 1/3 residual in each discordant term, giving a
  // ~0.41 dependence floor that never vanished for independent data.
  constexpr double indep_D  = 0.5;
  
  double disc_dep = std::min(1.0, std::max(0.0, std::abs(d0_Co - indep_Co) / indep_Co));
  double nd_disc  = std::abs(dpm_d0 - indep_D) / indep_D;
  
  return std::sqrt((disc_dep + nd_disc) / 2.0);
}

static DepResult NNS_dep_pair_core(const std::vector<double>& xv,
                                   const std::vector<double>& yv,
                                   const std::vector<uint64_t>& quad_xy,
                                   const std::vector<uint64_t>& quad_yx,
                                   bool asym) {
  int n = xv.size();
  
  bool cx = true, cy = true;
  for (int i = 1; i < n; ++i) {
    if (xv[i] != xv[0]) cx = false;
    if (yv[i] != yv[0]) cy = false;
    if (!cx && !cy) break;
  }
  if (cx || cy) return {0.0, 0.0};
  
  std::unordered_map<uint64_t, std::vector<int>> grp_xy;
  grp_xy.reserve(n);
  for (int i = 0; i < n; ++i)
    grp_xy[quad_xy[i]].push_back(i);
  
  std::unordered_map<uint64_t, std::vector<int>> grp_yx;
  grp_yx.reserve(n);
  for (int i = 0; i < n; ++i)
    grp_yx[quad_yx[i]].push_back(i);
  
  CellCopula global_cop = copula_signed(xv, yv);
  if (!std::isfinite(global_cop.signed_value) || !std::isfinite(global_cop.magnitude))
    global_cop = {0.0, 0.0};
  
  double corr_xy = 0.0, dep_xy = 0.0;
  for (auto& kv : grp_xy) {
    const auto& idx = kv.second;
    int nq = static_cast<int>(idx.size());
    if (nq < 1) continue;
    
    std::vector<double> xq(nq), yq(nq);
    for (int k = 0; k < nq; ++k) { xq[k] = xv[idx[k]]; yq[k] = yv[idx[k]]; }
    
    CellCopula cop = copula_signed(xq, yq);
    if (!std::isfinite(cop.signed_value) || !std::isfinite(cop.magnitude)) cop = global_cop;
    
    double w = static_cast<double>(nq) / static_cast<double>(n);
    corr_xy += cop.signed_value * w;
    dep_xy  += cop.magnitude    * w;
  }
  
  double corr_yx = 0.0, dep_yx = 0.0;
  for (auto& kv : grp_yx) {
    const auto& idx = kv.second;
    int nq = static_cast<int>(idx.size());
    if (nq < 1) continue;
    
    std::vector<double> yq(nq), xq(nq);
    for (int k = 0; k < nq; ++k) { yq[k] = yv[idx[k]]; xq[k] = xv[idx[k]]; }
    
    CellCopula cop = copula_signed(yq, xq);
    if (!std::isfinite(cop.signed_value) || !std::isfinite(cop.magnitude)) cop = global_cop;
    
    double w = static_cast<double>(nq) / static_cast<double>(n);
    corr_yx += cop.signed_value * w;
    dep_yx  += cop.magnitude    * w;
  }
  
  // Cells below chance can make the debiased average negative; dependence is a
  // magnitude, so it is floored at 0 after averaging (not per cell, which would
  // reintroduce a floor).
  // The signed Correlation is a signed dependence, so |Correlation| is kept
  // within that direction's Dependence (below-chance cells with a negative
  // slope would otherwise add to it).
  dep_xy = std::max(0.0, dep_xy);
  dep_yx = std::max(0.0, dep_yx);
  corr_xy = std::min(dep_xy, std::max(-dep_xy, corr_xy));
  corr_yx = std::min(dep_yx, std::max(-dep_yx, corr_yx));

  int lx = n_unique(xv);
  int ly = n_unique(yv);
  bool discrete_case = (lx < std::sqrt(static_cast<double>(n))) &&
    (ly < std::sqrt(static_cast<double>(n)));
  
  if (discrete_case) {
    double disc_cop = copula_degree0_unsigned(xv, yv);
    if (!std::isfinite(disc_cop)) disc_cop = std::max(dep_xy, dep_yx);
    
    if (asym) {
      std::vector<double> gv = {dep_xy, disc_cop};
      dep_xy = gravity_pure_cpp(gv);
    } else {
      double dep_sym = std::max(dep_xy, dep_yx);
      std::vector<double> gv = {dep_sym, disc_cop};
      double blended = gravity_pure_cpp(gv);
      dep_xy = blended;
      dep_yx = blended;
    }
  }
  
  if (asym) {
    return {corr_xy, dep_xy};
  }
  
  return {std::max(corr_xy, corr_yx), std::max(dep_xy, dep_yx)};
}

// [[Rcpp::export]]
List NNS_dep_pair_cpp(NumericVector x,
                      NumericVector y,
                      CharacterVector quad_xy,
                      CharacterVector quad_yx,
                      bool asym = false) {
  int n = x.size();
  if (y.size() != n || quad_xy.size() != n || quad_yx.size() != n) {
    stop("x, y, quad_xy, quad_yx must all have the same length");
  }
  
  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> yv(y.begin(), y.end());
  
  std::hash<std::string> hasher;
  std::vector<uint64_t> q_xy(n), q_yx(n);
  for (int i = 0; i < n; ++i) {
    q_xy[i] = hasher(std::string(quad_xy[i]));
    q_yx[i] = hasher(std::string(quad_yx[i]));
  }
  
  DepResult res = NNS_dep_pair_core(xv, yv, q_xy, q_yx, asym);
  
  return List::create(_["Correlation"] = res.correlation,
                      _["Dependence"]  = res.dependence);
}
