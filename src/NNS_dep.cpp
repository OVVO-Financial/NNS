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

static double copula_signed(const std::vector<double>& xv,
                            const std::vector<double>& yv) {
  int n = static_cast<int>(xv.size());
  if (n < 2) return 0.0;
  
  double tx = 0.0, ty = 0.0;
  for (int i = 0; i < n; ++i) { tx += xv[i]; ty += yv[i]; }
  tx /= n; ty /= n;
  
  double d0_cupm = 0.0, d0_clpm = 0.0, dpm_d0_count = 0.0;
  double c1_cupm = 0.0, c1_clpm = 0.0, c1_dpm = 0.0;
  double cov = 0.0, varx = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double dx = xv[i] - tx;
    double dy = yv[i] - ty;
    
    if (dx > 0.0 && dy > 0.0) d0_cupm += 1.0;
    if (dx <= 0.0 && dy <= 0.0) d0_clpm += 1.0;
    if (!((dx < 0.0 && dy < 0.0) || (dx > 0.0 && dy > 0.0)))
      dpm_d0_count += 1.0;
    
    if (dx >= 0.0 && dy >= 0.0) {
      c1_cupm += dx * dy;
    } else if (dx <= 0.0 && dy <= 0.0) {
      c1_clpm += dx * dy;
    } else {
      c1_dpm += std::abs(dx) * std::abs(dy);
    }
    
    cov += dx * dy;
    varx += dx * dx;
  }
  
  double inv_n = 1.0 / static_cast<double>(n);
  double d0_Co = (d0_cupm + d0_clpm) * inv_n;
  if (d0_Co == 1.0 || d0_Co == 0.0) return 1.0;
  
  double c1_total = c1_cupm + c1_clpm + c1_dpm;
  double co_d1 = c1_total > 0.0 ? (c1_cupm + c1_clpm) / c1_total : 0.0;
  double dpm_d0 = dpm_d0_count * inv_n;
  double dpm_d1 = c1_total > 0.0 ? c1_dpm / c1_total : 0.0;
  
  constexpr double indep_Co = 0.5;
  // Independence null for the discordant partial moment. DpmCountWorker treats
  // a point as concordant when it is all-below OR all-above the target (both
  // fully-aligned orthants), so under independence P(discordant) = 1 - 2*0.5^n,
  // which is 0.5 for the bivariate (n=2) copula -- not 1 - 0.5^n = 0.75. The
  // old 0.75 anchor left a fixed 1/3 residual in each discordant term, giving a
  // ~0.41 dependence floor that never vanished for independent data.
  constexpr double indep_D  = 0.5;
  
  double discrete_dep   = std::min(1.0, std::max(0.0, std::abs(d0_Co - indep_Co) / indep_Co));
  double continuous_dep = std::min(1.0, std::max(0.0, std::abs(co_d1 - indep_Co) / indep_Co));
  double nd_disc_dep    = std::abs(dpm_d0 - indep_D) / indep_D;
  double nd_cont_dep    = std::abs(dpm_d1 - indep_D) / indep_D;
  
  double copula_val = std::sqrt(
    (discrete_dep + continuous_dep + nd_disc_dep + nd_cont_dep) / 4.0
  );
  
  double slope_sign = varx == 0.0 ? 0.0 : ((cov > 0.0) ? 1.0 : (cov < 0.0) ? -1.0 : 0.0);
  return copula_val * slope_sign;
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
    if (!((dx < 0.0 && dy < 0.0) || (dx > 0.0 && dy > 0.0)))
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
  
  double global_cop = copula_signed(xv, yv);
  if (!std::isfinite(global_cop)) global_cop = 0.0;
  
  double corr_xy = 0.0, dep_xy = 0.0;
  for (auto& kv : grp_xy) {
    const auto& idx = kv.second;
    int nq = static_cast<int>(idx.size());
    if (nq < 1) continue;
    
    std::vector<double> xq(nq), yq(nq);
    for (int k = 0; k < nq; ++k) { xq[k] = xv[idx[k]]; yq[k] = yv[idx[k]]; }
    
    double cop = copula_signed(xq, yq);
    if (!std::isfinite(cop)) cop = global_cop;
    
    double w = static_cast<double>(nq) / static_cast<double>(n);
    corr_xy += cop           * w;
    dep_xy  += std::abs(cop) * w;
  }
  
  double corr_yx = 0.0, dep_yx = 0.0;
  for (auto& kv : grp_yx) {
    const auto& idx = kv.second;
    int nq = static_cast<int>(idx.size());
    if (nq < 1) continue;
    
    std::vector<double> yq(nq), xq(nq);
    for (int k = 0; k < nq; ++k) { yq[k] = yv[idx[k]]; xq[k] = xv[idx[k]]; }
    
    double cop = copula_signed(yq, xq);
    if (!std::isfinite(cop)) cop = global_cop;
    
    double w = static_cast<double>(nq) / static_cast<double>(n);
    corr_yx += cop           * w;
    dep_yx  += std::abs(cop) * w;
  }
  
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
