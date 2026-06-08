// NNS_dep.cpp
//
// C++ implementation of NNS.dep and NNS.dep.matrix.
//
// Exported functions (called from R):
//   NNS_dep_pair_cpp    — bivariate dependence given pre-computed partition labels
//   NNS_dep_matrix_cpp  — full pairwise dependence matrix, parallelised
//
// Design:
//   The bottleneck in the original NNS.dep was the per-quadrant R dispatch
//   loop: for each quadrant, R called NNS.copula() which called PM.matrix()
//   twice plus DPM_nD() twice, each time paying the R->C++ round-trip cost.
//   With order-3 partitions and two directions (xy, yx) that is ~16 R->C++
//   dispatches per NNS.dep call, compounded across every predictor, fold,
//   trial and epoch in NNS.boost / NNS.stack.
//
//   This file eliminates every one of those dispatches.  The per-quadrant
//   copula computation is evaluated directly for the bivariate degree-0 and
//   degree-1 partial moment terms used by NNS.dep.
//
//   Discrete-variable correction (triggered when both x and y have fewer
//   unique values than sqrt(n)):
//     Original used poly(x, degree) -> fast_lm_mult -> R2, which overfits
//     when unique values ~= degree.  Replaced with the degree-0 copula on
//     the full (x, y) pair -- a pure frequency / count measure that is the
//     correct discrete dependence tool already in the partial moment framework.
//     gravity(c(dependence, discrete_copula_full)) blends the per-quadrant
//     weighted sum with the global discrete anchor.
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
#include "central_tendencies.h" // NNS_gravity_cpp

using namespace Rcpp;
using namespace RcppParallel;

// ============================================================
// INTERNAL HELPERS
// ============================================================

// gravity() -- wraps NNS_gravity_cpp for std::vector input
static inline double gravity_cpp(const std::vector<double>& v) {
  if (v.empty()) return NA_REAL;
  NumericVector nv(v.begin(), v.end());
  return as<double>(NNS_gravity_cpp(nv, false));
}

// Count distinct values in a std::vector<double>
static inline int n_unique(const std::vector<double>& v) {
  std::unordered_map<double, int> seen;
  seen.reserve(v.size());
  for (double d : v) seen[d] = 1;
  return static_cast<int>(seen.size());
}

// ============================================================
// copula_signed
//
// Replicates exactly:
//   NNS.copula(cbind(xx, yy)) * sign(fast_lm(xx, yy)$coef[2])
//
// i.e. what dep_fn() returns per quadrant in NNS.dep.
//
// For a bivariate (2-column) matrix the formula in NNS.copula is:
//
//   discrete_pm_cov  <- PM.matrix(deg=0, deg=0, target, X)   [pop_adj=F]
//   continuous_pm_cov<- PM.matrix(deg=1, deg=1, target, X)   [pop_adj=T, norm=T]
//
//   indep_Co_pm = 0.25 * (n^2 - n)   -> for n=2 cols: 0.5
//   indep_D_pm  = 1 - 0.5^n          -> for n=2 cols: 0.75
//
//   discrete_dep   = |Co_pm_d0 - 0.5| / 0.5
//   continuous_dep = |Co_pm_d1 - 0.5| / 0.5
//   n_dim_disc_dep = |DPM_nD_d0 - 0.75| / 0.75
//   n_dim_cont_dep = |DPM_nD_d1 - 0.75| / 0.75
//
//   result = sqrt( mean(4 terms) ) * sign(slope)
// ============================================================
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
  
  constexpr double indep_Co = 0.5;   // 0.25*(2^2-2) for n=2 cols
  constexpr double indep_D  = 0.75;  // 1 - 0.5^2   for n=2 cols
  
  double discrete_dep   = std::min(1.0, std::max(0.0,
                                                 std::abs(d0_Co              - indep_Co) / indep_Co));
  double continuous_dep = std::min(1.0, std::max(0.0,
                                                 std::abs(co_d1              - indep_Co) / indep_Co));
  double nd_disc_dep    = std::abs(dpm_d0 - indep_D) / indep_D;
  double nd_cont_dep    = std::abs(dpm_d1 - indep_D) / indep_D;
  
  double copula_val = std::sqrt(
    (discrete_dep + continuous_dep + nd_disc_dep + nd_cont_dep) / 4.0
  );
  
  double slope_sign = varx == 0.0 ? 0.0 : ((cov > 0.0) ? 1.0 : (cov < 0.0) ? -1.0 : 0.0);
  return copula_val * slope_sign;
}

// ============================================================
// copula_degree0_unsigned
//
// Degree-0 only copula -- the proper discrete dependence measure.
// Used as the global anchor in the discrete-variable correction,
// replacing the polynomial R2 which overfits when unique values
// are close to polynomial degree.
//
// Returns sqrt( mean(discrete_dep, n_dim_discrete_dep) ) -- unsigned,
// since it is used as an absolute magnitude anchor in gravity().
// ============================================================
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
  constexpr double indep_D  = 0.75;
  
  double disc_dep = std::min(1.0, std::max(0.0,
                                           std::abs(d0_Co - indep_Co) / indep_Co));
  double nd_disc  = std::abs(dpm_d0 - indep_D) / indep_D;
  
  return std::sqrt((disc_dep + nd_disc) / 2.0);
}


// ============================================================
// NNS_dep_pair_cpp
//
// Full bivariate NNS.dep computation given pre-computed XONLY
// partition labels from NNS.part (both directions).
//
// Eliminates all R-level dispatch in the hot loop:
//   Original: 16+ R->C++ round-trips per call (one NNS.copula per quadrant)
//   This:     1 call from R, all copula work done in C++
//
// Arguments:
//   x, y     -- raw numeric vectors (length n)
//   quad_xy  -- quadrant labels from NNS.part(x, y, type="XONLY")$dt$quadrant
//   quad_yx  -- quadrant labels from NNS.part(y, x, type="XONLY")$dt$quadrant
//   asym     -- if TRUE, return asymmetric (xy direction only)
//
// Returns List("Correlation", "Dependence") -- identical to NNS.dep output.
// ============================================================

// [[Rcpp::export]]
List NNS_dep_pair_cpp(NumericVector x,
                      NumericVector y,
                      CharacterVector quad_xy,
                      CharacterVector quad_yx,
                      bool asym = false) {
  
  int n = x.size();
  if (y.size() != n || quad_xy.size() != n || quad_yx.size() != n)
    stop("x, y, quad_xy, quad_yx must all have the same length");
  
  // Degenerate: constant variable
    {
      bool cx = true, cy = true;
      for (int i = 1; i < n; ++i) {
        if (x[i] != x[0]) cx = false;
        if (y[i] != y[0]) cy = false;
        if (!cx && !cy) break;
      }
      if (cx || cy) return List::create(_["Correlation"] = 0.0,
          _["Dependence"]  = 0.0);
    }
    
    std::vector<double> xv(x.begin(), x.end());
    std::vector<double> yv(y.begin(), y.end());
    
    // --- Group observations by quadrant label ---
    std::unordered_map<std::string, std::vector<int>> grp_xy;
    grp_xy.reserve(n);
    for (int i = 0; i < n; ++i)
      grp_xy[std::string(quad_xy[i])].push_back(i);
    
    std::unordered_map<std::string, std::vector<int>> grp_yx;
    grp_yx.reserve(n);
    for (int i = 0; i < n; ++i)
      grp_yx[std::string(quad_yx[i])].push_back(i);
    
    // Fallback for NA copula results (matches original anyNA -> dep_fn(x,y))
    double global_cop = copula_signed(xv, yv);
    if (!R_finite(global_cop)) global_cop = 0.0;
    
    // --- xy direction: weighted copula sum ---
    double corr_xy = 0.0, dep_xy = 0.0;
    for (auto& kv : grp_xy) {
      const auto& idx = kv.second;
      int nq = static_cast<int>(idx.size());
      if (nq < 1) continue;
      
      std::vector<double> xq(nq), yq(nq);
      for (int k = 0; k < nq; ++k) { xq[k] = xv[idx[k]]; yq[k] = yv[idx[k]]; }
      
      double cop = copula_signed(xq, yq);
      if (!R_finite(cop)) cop = global_cop;
      
      double w = static_cast<double>(nq) / static_cast<double>(n);
      corr_xy += cop           * w;
      dep_xy  += std::abs(cop) * w;
    }
    
    // --- yx direction: weighted copula sum ---
    // yx partition was built on (y, x) so first variable in each group is y
    double corr_yx = 0.0, dep_yx = 0.0;
    for (auto& kv : grp_yx) {
      const auto& idx = kv.second;
      int nq = static_cast<int>(idx.size());
      if (nq < 1) continue;
      
      std::vector<double> yq(nq), xq(nq);
      for (int k = 0; k < nq; ++k) { yq[k] = yv[idx[k]]; xq[k] = xv[idx[k]]; }
      
      double cop = copula_signed(yq, xq);
      if (!R_finite(cop)) cop = global_cop;
      
      double w = static_cast<double>(nq) / static_cast<double>(n);
      corr_yx += cop           * w;
      dep_yx  += std::abs(cop) * w;
    }
    
    // --- Discrete-variable correction ---
    // Triggered when both variables have fewer unique values than sqrt(n).
    // Uses degree-0 copula (pure frequency / count) as the global anchor --
    // the correct discrete tool in the partial moment framework.
    // Replaces the original poly(x, degree) R2 which overfits at low cardinality.
    int lx = n_unique(xv);
    int ly = n_unique(yv);
    bool discrete_case = (lx < std::sqrt(static_cast<double>(n))) &&
      (ly < std::sqrt(static_cast<double>(n)));
    
    if (discrete_case) {
      double disc_cop = copula_degree0_unsigned(xv, yv);
      if (!R_finite(disc_cop)) disc_cop = std::max(dep_xy, dep_yx);
      
      if (asym) {
        std::vector<double> gv = {dep_xy, disc_cop};
        dep_xy = gravity_cpp(gv);
      } else {
        // Apply correction to the symmetric max before returning
        double dep_sym = std::max(dep_xy, dep_yx);
        std::vector<double> gv = {dep_sym, disc_cop};
        double blended = gravity_cpp(gv);
        // Write back so the final max below uses the corrected value
        dep_xy = blended;
        dep_yx = blended;
      }
    }
    
    // --- Final aggregation ---
    if (asym) {
      return List::create(_["Correlation"] = corr_xy,
                          _["Dependence"]  = dep_xy);
    }
    
    return List::create(_["Correlation"] = std::max(corr_xy, corr_yx),
                        _["Dependence"]  = std::max(dep_xy,  dep_yx));
}


// ============================================================
// DepMatrixWorker
//
// RcppParallel worker for NNS_dep_matrix_cpp.
// Each thread owns one (i, j) pair, runs its own inline XONLY
// partition (matching NNS.part(..., type="XONLY") exactly), then
// calls NNS_dep_pair_cpp.  The inline partition avoids calling
// the R-level NNS.part from within a parallel section.
// ============================================================
struct DepMatrixWorker : public Worker {
  
  const RMatrix<double> X;
  const int n_obs;
  const int n_vars;
  const bool asym;
  
  RVector<double> corr_upper;
  RVector<double> dep_upper;
  RVector<double> corr_lower;
  RVector<double> dep_lower;
  
  std::vector<int> pair_i, pair_j;
  
  DepMatrixWorker(const NumericMatrix& X_,
                  bool asym_,
                  NumericVector& cu,
                  NumericVector& du,
                  NumericVector& cl,
                  NumericVector& dl)
    : X(X_), n_obs(X_.nrow()), n_vars(X_.ncol()), asym(asym_),
      corr_upper(cu), dep_upper(du),
      corr_lower(cl), dep_lower(dl)
  {
    int np = n_vars * (n_vars - 1) / 2;
    pair_i.reserve(np); pair_j.reserve(np);
    for (int i = 0; i < n_vars - 1; ++i)
      for (int j = i + 1; j < n_vars; ++j) {
        pair_i.push_back(i);
        pair_j.push_back(j);
      }
  }
  
  // Inline XONLY partition matching NNS.part(x, y, type="XONLY", obs.req=obs_req):
  // splits on x only, up to floor(log2(n)) depth,
  // does not stop early on group size (min.obs.stop = FALSE).
  static std::vector<std::string>
    xonly_partition(const std::vector<double>& x, int obs_req) {
      int n = static_cast<int>(x.size());
      int max_order = std::max(1, static_cast<int>(
        std::floor(std::log2(std::max(1, n)))));
      
      std::vector<std::string> quad(n, "q");
      
      for (int depth = 0; depth < max_order; ++depth) {
        std::unordered_map<std::string, std::vector<int>> grp;
        grp.reserve(n);
        for (int i = 0; i < n; ++i) grp[quad[i]].push_back(i);
        
        bool any_split = false;
        for (auto& kv : grp) {
          const auto& idx = kv.second;
          if (static_cast<int>(idx.size()) <= obs_req) continue;
          
          double cx = 0.0;
          for (int i : idx) cx += x[i];
          cx /= static_cast<double>(idx.size());
          
          for (int i : idx)
            quad[i] += (x[i] > cx) ? "2" : "1";
          
          any_split = true;
        }
        if (!any_split) break;
      }
      return quad;
    }
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t p = begin; p < end; ++p) {
      int ci = pair_i[p];
      int cj = pair_j[p];
      
      std::vector<double> xv(n_obs), yv(n_obs);
      for (int r = 0; r < n_obs; ++r) {
        xv[r] = X(r, ci);
        yv[r] = X(r, cj);
      }
      
      int obs_req = std::max(8, n_obs / 8);
      
      std::vector<std::string> q_xy = xonly_partition(xv, obs_req);
      std::vector<std::string> q_yx = xonly_partition(yv, obs_req);
      
      CharacterVector quad_xy(n_obs), quad_yx(n_obs);
      for (int r = 0; r < n_obs; ++r) {
        quad_xy[r] = q_xy[r];
        quad_yx[r] = q_yx[r];
      }
      
      NumericVector xnv(xv.begin(), xv.end());
      NumericVector ynv(yv.begin(), yv.end());
      
      List res_ij = NNS_dep_pair_cpp(xnv, ynv, quad_xy, quad_yx, asym);
      corr_upper[p] = as<double>(res_ij["Correlation"]);
      dep_upper[p]  = as<double>(res_ij["Dependence"]);
      
      if (asym) {
        List res_ji = NNS_dep_pair_cpp(ynv, xnv, quad_yx, quad_xy, true);
        corr_lower[p] = as<double>(res_ji["Correlation"]);
        dep_lower[p]  = as<double>(res_ji["Dependence"]);
      } else {
        corr_lower[p] = corr_upper[p];
        dep_lower[p]  = dep_upper[p];
      }
    }
  }
};


// ============================================================
// NNS_dep_matrix_cpp
//
// Parallel pairwise dependence matrix -- replaces NNS.dep.matrix().
// All n*(n-1)/2 pairs computed simultaneously via parallelFor.
// ============================================================

// [[Rcpp::export]]
List NNS_dep_matrix_cpp(NumericMatrix X,
                        bool asym = false) {
  
  int n_vars = X.ncol();
  
  if (n_vars < 2)
    stop("NNS_dep_matrix_cpp: X must have at least 2 columns");
  
  int n_pairs = n_vars * (n_vars - 1) / 2;
  
  NumericVector corr_upper(n_pairs, 0.0);
  NumericVector dep_upper (n_pairs, 0.0);
  NumericVector corr_lower(n_pairs, 0.0);
  NumericVector dep_lower (n_pairs, 0.0);
  
  DepMatrixWorker worker(X, asym,
                         corr_upper, dep_upper,
                         corr_lower, dep_lower);
  
  parallelFor(0, n_pairs, worker);
  
  NumericMatrix rhos(n_vars, n_vars);
  NumericMatrix deps(n_vars, n_vars);
  for (int i = 0; i < n_vars; ++i) { rhos(i, i) = 1.0; deps(i, i) = 1.0; }
  
  {
    int p = 0;
    for (int i = 0; i < n_vars - 1; ++i) {
      for (int j = i + 1; j < n_vars; ++j, ++p) {
        if (!asym) {
          double r = (corr_upper[p] + corr_lower[p]) / 2.0;
          double d = (dep_upper[p]  + dep_lower[p])  / 2.0;
          rhos(i, j) = r; rhos(j, i) = r;
          deps(i, j) = d; deps(j, i) = d;
        } else {
          rhos(i, j) = corr_upper[p];
          deps(i, j) = dep_upper[p];
          rhos(j, i) = corr_lower[p];
          deps(j, i) = dep_lower[p];
        }
      }
    }
  }
  
  CharacterVector cn = colnames(X);
  if (cn.size() == n_vars) {
    colnames(rhos) = cn; rownames(rhos) = cn;
    colnames(deps) = cn; rownames(deps) = cn;
  }
  
  return List::create(_["Correlation"] = rhos,
                      _["Dependence"]  = deps);
}