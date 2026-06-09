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

static inline double gravity_cpp(const std::vector<double>& v) {
  if (v.empty()) return NA_REAL;
  NumericVector nv(v.begin(), v.end());
  return as<double>(NNS_gravity_cpp(nv, false));
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
  constexpr double indep_D  = 0.75;
  
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
  constexpr double indep_D  = 0.75;
  
  double disc_dep = std::min(1.0, std::max(0.0, std::abs(d0_Co - indep_Co) / indep_Co));
  double nd_disc  = std::abs(dpm_d0 - indep_D) / indep_D;
  
  return std::sqrt((disc_dep + nd_disc) / 2.0);
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
  
  {
    bool cx = true, cy = true;
    for (int i = 1; i < n; ++i) {
      if (x[i] != x[0]) cx = false;
      if (y[i] != y[0]) cy = false;
      if (!cx && !cy) break;
    }
    if (cx || cy) return List::create(_["Correlation"] = 0.0, _["Dependence"]  = 0.0);
  }
  
  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> yv(y.begin(), y.end());
  
  std::unordered_map<std::string, std::vector<int>> grp_xy;
  grp_xy.reserve(n);
  for (int i = 0; i < n; ++i)
    grp_xy[std::string(quad_xy[i])].push_back(i);
  
  std::unordered_map<std::string, std::vector<int>> grp_yx;
  grp_yx.reserve(n);
  for (int i = 0; i < n; ++i)
    grp_yx[std::string(quad_yx[i])].push_back(i);
  
  double global_cop = copula_signed(xv, yv);
  if (!R_finite(global_cop)) global_cop = 0.0;
  
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
      double dep_sym = std::max(dep_xy, dep_yx);
      std::vector<double> gv = {dep_sym, disc_cop};
      double blended = gravity_cpp(gv);
      dep_xy = blended;
      dep_yx = blended;
    }
  }
  
  if (asym) {
    return List::create(_["Correlation"] = corr_xy, _["Dependence"]  = dep_xy);
  }
  
  return List::create(_["Correlation"] = std::max(corr_xy, corr_yx),
                      _["Dependence"]  = std::max(dep_xy,  dep_yx));
}

// Worker 1: Parallel Precomputation of Space Partition Labels
struct PrecomputePartitionsWorker : public Worker {
  const RMatrix<double> X;
  const int n_obs;
  const int obs_req;
  std::vector<std::vector<std::string>>& all_quads;
  
  PrecomputePartitionsWorker(const NumericMatrix& X_, int obs_req_, std::vector<std::vector<std::string>>& all_quads_)
    : X(X_), n_obs(X_.nrow()), obs_req(obs_req_), all_quads(all_quads_) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t j = begin; j < end; ++j) {
      int max_order = std::max(1, static_cast<int>(std::floor(std::log2(std::max(1, n_obs)))));
      std::vector<std::string> quad(n_obs, "q");
      
      for (int depth = 0; depth < max_order; ++depth) {
        std::unordered_map<std::string, std::vector<int>> grp;
        grp.reserve(n_obs);
        for (int i = 0; i < n_obs; ++i) grp[quad[i]].push_back(i);
        
        bool any_split = false;
        for (auto& kv : grp) {
          const auto& idx = kv.second;
          if (static_cast<int>(idx.size()) <= obs_req) continue;
          
          double cx = 0.0;
          for (int i : idx) cx += X(i, j);
          cx /= static_cast<double>(idx.size());
          
          for (int i : idx)
            quad[i] += (X(i, j) > cx) ? "2" : "1";
          
          any_split = true;
        }
        if (!any_split) break;
      }
      all_quads[j] = std::move(quad);
    }
  }
};

// Worker 2: Concurrent Evaluation of Cross-Product Asset Dependencies
struct DepMatrixWorker : public Worker {
  const RMatrix<double> X;
  const int n_obs;
  const int n_vars;
  const bool asym;
  const std::vector<std::vector<std::string>>& all_quads;
  
  RVector<double> corr_upper;
  RVector<double> dep_upper;
  RVector<double> corr_lower;
  RVector<double> dep_lower;
  
  std::vector<int> pair_i, pair_j;
  
  DepMatrixWorker(const NumericMatrix& X_,
                  bool asym_,
                  const std::vector<std::vector<std::string>>& all_quads_,
                  NumericVector& cu,
                  NumericVector& du,
                  NumericVector& cl,
                  NumericVector& dl)
    : X(X_), n_obs(X_.nrow()), n_vars(X_.ncol()), asym(asym_), all_quads(all_quads_),
      corr_upper(cu), dep_upper(du), corr_lower(cl), dep_lower(dl)
  {
    int np = n_vars * (n_vars - 1) / 2;
    pair_i.reserve(np); pair_j.reserve(np);
    for (int i = 0; i < n_vars - 1; ++i) {
      for (int j = i + 1; j < n_vars; ++j) {
        pair_i.push_back(i);
        pair_j.push_back(j);
      }
    }
  }
  
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t p = begin; p < end; ++p) {
      int ci = pair_i[p];
      int cj = pair_j[p];
      
      const std::vector<std::string>& q_xy = all_quads[ci];
      const std::vector<std::string>& q_yx = all_quads[cj];
      
      CharacterVector quad_xy(n_obs), quad_yx(n_obs);
      NumericVector xnv(n_obs), ynv(n_obs);
      
      for (int r = 0; r < n_obs; ++r) {
        quad_xy[r] = q_xy[r];
        quad_yx[r] = q_yx[r];
        xnv[r] = X(r, ci);
        ynv[r] = X(r, cj);
      }
      
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

// [[Rcpp::export]]
List NNS_dep_matrix_cpp(NumericMatrix X, bool asym = false) {
  int n_vars = X.ncol();
  int n_obs  = X.nrow();
  if (n_vars < 2)
    stop("NNS_dep_matrix_cpp: X must have at least 2 columns");
  
  int n_pairs = n_vars * (n_vars - 1) / 2;
  
  int obs_req = std::max(8, n_obs / 8);
  std::vector<std::vector<std::string>> all_quads(n_vars);
  PrecomputePartitionsWorker partitioner(X, obs_req, all_quads);
  parallelFor(0, n_vars, partitioner);
  
  NumericVector corr_upper(n_pairs, 0.0);
  NumericVector dep_upper (n_pairs, 0.0);
  NumericVector corr_lower(n_pairs, 0.0);
  NumericVector dep_lower (n_pairs, 0.0);
  
  DepMatrixWorker worker(X, asym, all_quads, corr_upper, dep_upper, corr_lower, dep_lower);
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
  
  return List::create(_["Correlation"] = rhos, _["Dependence"]  = deps);
}
