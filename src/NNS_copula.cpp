// NNS.copula's numeric core in one call.
//
// R/Copula.R used to reach the partial-moment routines through about ten
// separate R calls with R arithmetic in between. NNS.reg's default order
// calls NNS.copula twice per regressor, so that bookkeeping, not the
// arithmetic, was most of a short series' regression time.
//
// The kernel calls the same routines, in the same order, on the same values,
// and repeats R's own arithmetic for the steps in between: sum(), prod() and
// mean() accumulate in long double as R's do (mean() with R's second
// refinement pass), and `^` is R_pow(). The result is bit-identical to the R
// implementation.

// Keep a * b + c * d as R evaluates it, unfused, on every compiler.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include "partial_moments_rcpp.h"

using namespace Rcpp;

namespace {

// sum() of the strict upper triangle in column order, as sum(m[upper.tri(m)]).
double upper_sum(const NumericMatrix& m) {
  const int n = m.ncol();
  long double s = 0.0L;
  for (int j = 1; j < n; ++j)
    for (int i = 0; i < j; ++i) s += m(i, j);
  return static_cast<double>(s);
}

struct Anchors { double Co; double D; };

// list(Co = sum((outer(p, p) + outer(1 - p, 1 - p))[utr]),
//      D  = 1 - prod(p) - prod(1 - p))
Anchors indep_anchors(const std::vector<double>& p) {
  const int n = static_cast<int>(p.size());
  std::vector<double> q(n);
  for (int i = 0; i < n; ++i) q[i] = 1 - p[i];
  long double s = 0.0L;
  for (int j = 1; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      const double pp = p[i] * p[j];
      const double qq = q[i] * q[j];
      s += pp + qq;
    }
  }
  long double prod_p = 1.0L, prod_q = 1.0L;
  for (int i = 0; i < n; ++i) { prod_p *= p[i]; prod_q *= q[i]; }
  Anchors a;
  a.Co = static_cast<double>(s);
  a.D = (1 - static_cast<double>(prod_p)) - static_cast<double>(prod_q);
  return a;
}

// if (obs == indep) 0 else abs(obs - indep) / indep
double rel_dev(double obs, double indep) {
  if (ISNAN(obs) || ISNAN(indep)) stop("missing value where TRUE/FALSE needed");
  if (obs == indep) return 0.0;
  return std::fabs(obs - indep) / indep;
}

// R's max() and min() of two numbers: a missing value wins.
double r_max(double a, double b) {
  if (ISNAN(a)) return a;
  if (ISNAN(b)) return b;
  return a > b ? a : b;
}
double r_min(double a, double b) {
  if (ISNAN(a)) return a;
  if (ISNAN(b)) return b;
  return a < b ? a : b;
}

// R's mean() of a short vector.
double r_mean(const double* v, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += v[i];
  if (!R_FINITE(static_cast<double>(s))) {
    s = 0.0L;
    for (int i = 0; i < n; ++i) s += v[i] / n;
    return static_cast<double>(s);
  }
  s /= n;
  long double t = 0.0L;
  for (int i = 0; i < n; ++i) t += (v[i] - s);
  s += t / n;
  return static_cast<double>(s);
}

}  // namespace

// [[Rcpp::export]]
double NNS_copula_cpp(NumericMatrix X, NumericVector target, bool continuous) {
  int n = X.ncol();
  const int rows = X.nrow();

  // A constant margin is independent of the rest and has no degree-1 lower share.
  std::vector<int> keep;
  for (int j = 0; j < n; ++j) {
    bool vary = false;
    for (int i = 1; i < rows; ++i) {
      if (X(i, j) != X(0, j)) { vary = true; break; }
    }
    if (vary) keep.push_back(j);
  }
  if (static_cast<int>(keep.size()) != n) {
    if (keep.size() < 2) return 0.0;
    NumericMatrix Xv(rows, static_cast<int>(keep.size()));
    NumericVector tv(keep.size());
    for (std::size_t a = 0; a < keep.size(); ++a) {
      for (int i = 0; i < rows; ++i) Xv(i, a) = X(i, keep[a]);
      tv[a] = target[keep[a]];
    }
    X = Xv;
    target = tv;
    n = X.ncol();
  }

  // Pairwise
  List discrete = PMMatrix_RCPP(0.0, 0.0, target, X, false, false);
  NumericMatrix d_cupm = discrete["cupm"], d_clpm = discrete["clpm"];
  const double discrete_Co_pm = upper_sum(d_cupm) + upper_sum(d_clpm);

  // Perfect dependence: every pair fully concordant (1) or fully discordant (0).
  // A margin with every observation on one side of its target makes that trivially true.
  std::vector<NumericVector> column(n);
  std::vector<double> lower_share(n);
  for (int j = 0; j < n; ++j) {
    column[j] = X.column(j);
    lower_share[j] = LPM_RCPP(0.0, NumericVector::create(target[j]), column[j], false)[0];
  }
  // Normalized region shares have no informative contrast at this target.
  bool all_edge = true, all_inside = true;
  for (int j = 0; j < n; ++j) {
    const double p = lower_share[j];
    if (!(p == 0 || p == 1)) all_edge = false;
    if (!(p > 0 && p < 1)) all_inside = false;
  }
  if (all_edge) return 0.0;
  if (all_inside) {
    bool extreme = true;
    for (int j = 1; j < n && extreme; ++j) {
      for (int i = 0; i < j; ++i) {
        const double pair = d_cupm(i, j) + d_clpm(i, j);
        if (!(std::fabs(pair - 1) < 1e-9 || std::fabs(pair) < 1e-9)) { extreme = false; break; }
      }
    }
    if (extreme) return 1.0;
  }

  double continuous_Co_pm = discrete_Co_pm;
  if (continuous) {
    List cont = PMMatrix_RCPP(1.0, 1.0, target, X, true, true);
    NumericMatrix c_cupm = cont["cupm"], c_clpm = cont["clpm"];
    continuous_Co_pm = upper_sum(c_cupm) + upper_sum(c_clpm);
  }

  // Independence anchors from each margin's own lower-tail share, so they hold at any target
  const Anchors discrete_indep = indep_anchors(lower_share);
  Anchors continuous_indep = discrete_indep;
  if (continuous) {
    std::vector<double> ratio(n);
    for (int j = 0; j < n; ++j) {
      ratio[j] = LPM_ratio_RCPP(1.0, NumericVector::create(target[j]), column[j])[0];
    }
    continuous_indep = indep_anchors(ratio);
  }

  // A boundary target can zero an anchor; an observed 0 then matches it exactly
  double discrete_dep = rel_dev(discrete_Co_pm, discrete_indep.Co);
  double continuous_dep = rel_dev(continuous_Co_pm, continuous_indep.Co);
  discrete_dep = r_min(r_max(discrete_dep, 0), 1);
  continuous_dep = r_min(r_max(continuous_dep, 0), 1);

  // n-dimensional
  const double discrete_D_pm = DPM_nD_RCPP(X, target, 0.0, true);
  const double continuous_D_pm = continuous ? DPM_nD_RCPP(X, target, 1.0, true) : discrete_D_pm;
  double n_dim_discrete_dep = rel_dev(discrete_D_pm, discrete_indep.D);
  double n_dim_continuous_dep = rel_dev(continuous_D_pm, continuous_indep.D);
  n_dim_discrete_dep = r_min(n_dim_discrete_dep, 1);
  n_dim_continuous_dep = r_min(n_dim_continuous_dep, 1);

  const double parts[4] = {discrete_dep, continuous_dep, n_dim_discrete_dep, n_dim_continuous_dep};
  return R_pow(r_mean(parts, 4), 0.5);
}
