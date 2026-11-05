// table_ops.cpp
// Rcpp helpers that replace the former data.table grouped operations, so the
// R routines stay fully backed by compiled code (no per-call R grouping loops).

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include "central_tendencies.h"

using namespace Rcpp;

// gravity() of a vector, filtering non-finite values first (matches
// NNS_gravity_cpp / NNS.gravity exactly).
static inline double grav_finite(const std::vector<double>& v) {
  std::vector<double> f;
  f.reserve(v.size());
  for (double d : v) if (R_finite(d)) f.push_back(d);
  return gravity_value(f, false);
}

// Grouped regression-point reduction for NNS.M.reg.
//
// Replaces the (former data.table) "group by NNS.ID, reduce each IV column and
// the DV by the chosen central tendency" stage plus the residual-bias step.
//
//   X          : n x m matrix of independent variables (finite numeric)
//   dv         : length-n dependent variable
//   ids        : length-n group labels (NNS.ID)
//   reducer    : 0 gravity (off), 1 mean, 2 median, 3 mode, 4 mode_class,
//                5 identity (first row of group; order is non-numeric / "max")
//   class_mode : round fitted y.hat toward the nearest class (type = "CLASS")
//
// Returns list(rpm, yhat, residuals):
//   rpm        : G x (m+1); one row per unique NNS.ID in ascending (C-locale)
//                id order; columns = reduced IVs then reduced DV ("y.hat").
//   yhat       : length-n fitted values in ORIGINAL observation order, after the
//                grouped residual-bias correction.
//   residuals  : length-n (yhat - dv) after the bias correction.
//
// [[Rcpp::export]]
List NNS_mreg_reduce_cpp(NumericMatrix X, NumericVector dv,
                         CharacterVector ids, int reducer, bool class_mode) {
  const int n = X.nrow();
  const int m = X.ncol();

  std::vector<std::string> id(n);
  for (int i = 0; i < n; ++i) id[i] = as<std::string>(ids[i]);

  // unique ids in ascending byte (C-locale) order -- matches order(method = "radix")
  std::vector<std::string> uid(id);
  std::sort(uid.begin(), uid.end());
  uid.erase(std::unique(uid.begin(), uid.end()), uid.end());
  const int G = (int) uid.size();

  std::unordered_map<std::string, int> pos;
  pos.reserve(G * 2 + 1);
  for (int g = 0; g < G; ++g) pos[uid[g]] = g;

  std::vector<int> grp(n);
  std::vector< std::vector<int> > rows(G);
  for (int i = 0; i < n; ++i) { const int g = pos[id[i]]; grp[i] = g; rows[g].push_back(i); }

  NumericMatrix rpm(G, m + 1);
  for (int g = 0; g < G; ++g) {
    const std::vector<int>& r = rows[g];
    for (int j = 0; j <= m; ++j) {
      std::vector<double> col;
      col.reserve(r.size());
      for (size_t t = 0; t < r.size(); ++t) col.push_back(j < m ? X(r[t], j) : dv[r[t]]);

      double val;
      if (reducer == 5) {                              // identity: first row
        val = col.front();
      } else if (reducer == 1) {                       // mean
        long double s = 0.0L; for (double d : col) s += d;
        val = (double)(s / (long double) col.size());
      } else if (reducer == 2) {                       // median
        std::vector<double> t2(col);
        std::sort(t2.begin(), t2.end());
        const int L = (int) t2.size();
        val = (L % 2) ? t2[L / 2] : 0.5 * (t2[L / 2 - 1] + t2[L / 2]);
      } else if (reducer == 3 || reducer == 4) {       // mode / mode_class
        NumericVector mr(NNS_mode_cpp(wrap(col), reducer == 4, false));
        val = mr.size() ? mr[0] : NA_REAL;
      } else {                                         // gravity (off / default)
        val = grav_finite(col);
      }
      rpm(g, j) = val;
    }
  }

  // initial fitted y.hat in original order = each obs's group reduced DV
  NumericVector yhat(n);
  for (int i = 0; i < n; ++i) yhat[i] = rpm(grp[i], m);
  if (class_mode) {
    for (int i = 0; i < n; ++i) {
      const double v = yhat[i], f = v - std::floor(v);
      yhat[i] = (f < 0.5) ? std::floor(v) : std::ceil(v);
    }
  }

  // residual-bias correction by NNS.ID: bias = gravity(residuals); y.hat -= bias
  std::vector<double> resid(n);
  for (int i = 0; i < n; ++i) resid[i] = yhat[i] - dv[i];
  for (int g = 0; g < G; ++g) {
    const std::vector<int>& r = rows[g];
    std::vector<double> rr;
    rr.reserve(r.size());
    for (size_t t = 0; t < r.size(); ++t) rr.push_back(resid[r[t]]);
    const double b = grav_finite(rr);
    for (size_t t = 0; t < r.size(); ++t) yhat[r[t]] -= b;
  }

  NumericVector residuals(n);
  for (int i = 0; i < n; ++i) residuals[i] = yhat[i] - dv[i];

  return List::create(_["rpm"] = rpm, _["yhat"] = yhat, _["residuals"] = residuals);
}
