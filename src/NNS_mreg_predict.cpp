// Native kernel for the repaired NNS.M.reg prediction rule.
// Replicates R/Multivariate_Regression.R's .nns_mreg_predict_one exactly:
//   - range-normalized L2 / L1 or Hamming (FACTOR) distances,
//   - stable ascending ordering with index tie-break,
//   - k = 1 aggregates all exact nearest-distance ties (gravity / class mode),
//   - k > 1 uses the eight-component ensemble weights built from Rmath
//     densities, matching the R helper term for term,
//   - exact double-precision weighted class mode (max total, tie -> min value).
// Serial by design: no global RcppParallel thread state is touched.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// central_tendencies.cpp
double gravity_value(std::vector<double> x, bool discrete);

namespace {

// .nns_mreg_normalize_weights: nonfinite/negative -> 0; sum<=0 -> uniform.
inline void normalize_weights(std::vector<double>& w) {
  const std::size_t k = w.size();
  double s = 0.0;
  for (double& v : w) {
    if (!R_finite(v) || v < 0.0) v = 0.0;
    s += v;
  }
  if (s <= 0.0) {
    for (double& v : w) v = 1.0 / static_cast<double>(k);
  } else {
    for (double& v : w) v /= s;
  }
}

inline double sample_sd(const std::vector<double>& v) {
  const std::size_t n = v.size();
  if (n < 2) return NA_REAL;
  double mu = 0.0;
  for (double x : v) mu += x;
  mu /= static_cast<double>(n);
  double acc = 0.0;
  for (double x : v) { const double d = x - mu; acc += d * d; }
  return std::sqrt(acc / static_cast<double>(n - 1));
}

// .nns_mreg_weighted_mode: totals per distinct value, max total, tie -> min.
inline double weighted_mode(const std::vector<double>& vals,
                            const std::vector<double>& w) {
  std::vector<std::pair<double, double> > items;
  items.reserve(vals.size());
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if (R_finite(vals[i]) && R_finite(w[i]) && w[i] >= 0.0) {
      items.push_back(std::make_pair(vals[i], w[i]));
    }
  }
  if (items.empty()) return NA_REAL;
  std::sort(items.begin(), items.end());
  double best_val = items[0].first;
  double best_tot = -1.0;
  std::size_t i = 0;
  while (i < items.size()) {
    const double v = items[i].first;
    double tot = 0.0;
    while (i < items.size() && items[i].first == v) { tot += items[i].second; ++i; }
    if (tot > best_tot) { best_tot = tot; best_val = v; }
  }
  return best_val;
}

// .nns_mreg_ensemble_weights for the k selected (ascending) distances.
inline std::vector<double> ensemble_weights(const std::vector<double>& d) {
  const int k = static_cast<int>(d.size());
  std::vector<double> total(k, 0.0);
  const double kk = static_cast<double>(k);

  std::vector<double> comp(k);

  // uniform
  for (int i = 0; i < k; ++i) total[i] += 1.0 / kk;

  // student t density of distances, df = k
  for (int i = 0; i < k; ++i) comp[i] = ::Rf_dt(d[i], kk, 0);
  normalize_weights(comp);
  for (int i = 0; i < k; ++i) total[i] += comp[i];

  // inverse distance
  for (int i = 0; i < k; ++i) comp[i] = 1.0 / std::max(d[i], 1e-12);
  normalize_weights(comp);
  for (int i = 0; i < k; ++i) total[i] += comp[i];

  // exponential density of ranks, rate = 1/k (Rf_dexp takes the scale = k)
  for (int i = 0; i < k; ++i) comp[i] = ::Rf_dexp(static_cast<double>(i + 1), kk, 0);
  normalize_weights(comp);
  for (int i = 0; i < k; ++i) total[i] += comp[i];

  // reversed |log lognormal density| of ranks, sdlog = sample sd of ranks
  {
    std::vector<double> ranks(k);
    for (int i = 0; i < k; ++i) ranks[i] = static_cast<double>(i + 1);
    const double rank_sd = sample_sd(ranks);
    if (R_finite(rank_sd) && rank_sd > 0.0) {
      for (int i = 0; i < k; ++i)
        comp[i] = std::fabs(::Rf_dlnorm(ranks[i], 0.0, rank_sd, 1));
      normalize_weights(comp);
      std::reverse(comp.begin(), comp.end());
      for (int i = 0; i < k; ++i) total[i] += comp[i];
    }
  }

  // power-law on ranks
  for (int i = 0; i < k; ++i) {
    const double r = static_cast<double>(i + 1);
    comp[i] = 1.0 / (r * r);
  }
  normalize_weights(comp);
  for (int i = 0; i < k; ++i) total[i] += comp[i];

  // normal density of distances, sd = sample sd of distances
  {
    const double dist_sd = sample_sd(d);
    if (R_finite(dist_sd) && dist_sd > 0.0) {
      for (int i = 0; i < k; ++i) comp[i] = ::Rf_dnorm4(d[i], 0.0, dist_sd, 0);
      normalize_weights(comp);
      for (int i = 0; i < k; ++i) total[i] += comp[i];
    }
  }

  // RBF on distances, bandwidth = 2 * sample variance
  {
    const double dist_sd = sample_sd(d);
    const double dist_var = R_finite(dist_sd) ? dist_sd * dist_sd : NA_REAL;
    if (R_finite(dist_var) && dist_var > 0.0) {
      for (int i = 0; i < k; ++i) comp[i] = std::exp(-d[i] / (2.0 * dist_var));
      normalize_weights(comp);
      for (int i = 0; i < k; ++i) total[i] += comp[i];
    }
  }

  normalize_weights(total);
  return total;
}

}  // namespace

namespace {

inline void row_distances(const NumericMatrix& rpm_x, const NumericMatrix& Xtest,
                          int r, int dist_code,
                          const std::vector<int>& active,
                          const std::vector<double>& inv_range,
                          std::vector<double>& d) {
  const int n = rpm_x.nrow(), p = rpm_x.ncol();
  if (dist_code == 2) {
    for (int i = 0; i < n; ++i) {
      double mismatches = 0.0;
      for (int j = 0; j < p; ++j)
        if (rpm_x(i, j) != Xtest(r, j)) mismatches += 1.0;
      d[i] = mismatches / static_cast<double>(p);
    }
  } else if (active.empty()) {
    std::fill(d.begin(), d.end(), 0.0);
  } else if (dist_code == 1) {
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      for (std::size_t a = 0; a < active.size(); ++a) {
        const int j = active[a];
        acc += std::fabs((rpm_x(i, j) - Xtest(r, j)) * inv_range[a]);
      }
      d[i] = acc;
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      for (std::size_t a = 0; a < active.size(); ++a) {
        const int j = active[a];
        const double z = (rpm_x(i, j) - Xtest(r, j)) * inv_range[a];
        acc += z * z;
      }
      d[i] = std::sqrt(acc);
    }
  }
}

inline double aggregate_min_ties(const std::vector<double>& d,
                                 const NumericVector& yhat,
                                 bool is_class) {
  const int n = static_cast<int>(d.size());
  double dmin = d[0];
  for (int i = 1; i < n; ++i) if (d[i] < dmin) dmin = d[i];
  std::vector<double> tied;
  for (int i = 0; i < n; ++i) if (d[i] == dmin) tied.push_back(yhat[i]);
  if (is_class) {
    std::vector<double> ones(tied.size(), 1.0);
    return weighted_mode(tied, ones);
  }
  std::vector<double> finite_tied;
  finite_tied.reserve(tied.size());
  for (double v : tied) if (R_finite(v)) finite_tied.push_back(v);
  return gravity_value(finite_tied, false);
}

inline void active_columns(const NumericVector& mins, const NumericVector& maxs,
                           std::vector<int>& active, std::vector<double>& inv_range) {
  const int p = mins.size();
  for (int j = 0; j < p; ++j) {
    const double range = maxs[j] - mins[j];
    if (R_finite(range) && range > 0.0) {
      active.push_back(j);
      inv_range.push_back(1.0 / range);
    }
  }
}

}  // namespace

// Prediction path for every k = 1..kmax under the same rule as
// NNS_mreg_predict_cpp, so cross-validated n.best selection scores exactly
// the estimator the final NNS.reg fit will use.
// [[Rcpp::export]]
NumericMatrix NNS_mreg_predict_path_cpp(const NumericMatrix& rpm_x,
                                        const NumericVector& yhat,
                                        const NumericMatrix& Xtest,
                                        int kmax,
                                        int dist_code,
                                        const NumericVector& mins,
                                        const NumericVector& maxs,
                                        bool is_class) {
  const int n = rpm_x.nrow(), p = rpm_x.ncol(), m = Xtest.nrow();
  if (yhat.size() != n) stop("yhat length must equal nrow(rpm_x)");
  if (Xtest.ncol() != p) stop("Xtest and rpm_x must have the same columns");
  if (mins.size() != p || maxs.size() != p) stop("mins/maxs must have one value per column");
  if (kmax < 1) stop("kmax must be >= 1");
  if (kmax > n) kmax = n;

  std::vector<int> active;
  std::vector<double> inv_range;
  active_columns(mins, maxs, active, inv_range);

  NumericMatrix out(m, kmax);
  std::vector<double> d(n);
  std::vector<int> idx(n);

  for (int r = 0; r < m; ++r) {
    row_distances(rpm_x, Xtest, r, dist_code, active, inv_range, d);

    for (int i = 0; i < n; ++i) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(),
                     [&d](int a, int b) { return d[a] < d[b]; });

    std::vector<double> dk(kmax), yk(kmax);
    for (int i = 0; i < kmax; ++i) {
      dk[i] = d[idx[i]];
      yk[i] = yhat[idx[i]];
    }

    out(r, 0) = aggregate_min_ties(d, yhat, is_class);
    for (int k = 2; k <= kmax; ++k) {
      const std::vector<double> dsub(dk.begin(), dk.begin() + k);
      const std::vector<double> w = ensemble_weights(dsub);
      if (is_class) {
        const std::vector<double> ysub(yk.begin(), yk.begin() + k);
        out(r, k - 1) = weighted_mode(ysub, w);
      } else {
        double dot = 0.0;
        for (int i = 0; i < k; ++i) dot += yk[i] * w[i];
        out(r, k - 1) = dot;
      }
    }
  }
  return out;
}

// dist_code: 0 = L2, 1 = L1, 2 = FACTOR (Hamming over encoded columns).
// [[Rcpp::export]]
NumericVector NNS_mreg_predict_cpp(const NumericMatrix& rpm_x,
                                   const NumericVector& yhat,
                                   const NumericMatrix& Xtest,
                                   int k,
                                   int dist_code,
                                   const NumericVector& mins,
                                   const NumericVector& maxs,
                                   bool is_class) {
  const int n = rpm_x.nrow(), p = rpm_x.ncol(), m = Xtest.nrow();
  if (yhat.size() != n) stop("yhat length must equal nrow(rpm_x)");
  if (Xtest.ncol() != p) stop("Xtest and rpm_x must have the same columns");
  if (mins.size() != p || maxs.size() != p) stop("mins/maxs must have one value per column");
  if (k < 1) stop("k must be >= 1");

  // Active columns and reciprocal ranges for the normalized metrics.
  std::vector<int> active;
  std::vector<double> inv_range;
  for (int j = 0; j < p; ++j) {
    const double range = maxs[j] - mins[j];
    if (R_finite(range) && range > 0.0) {
      active.push_back(j);
      inv_range.push_back(1.0 / range);
    }
  }

  NumericVector out(m);
  std::vector<double> d(n);
  std::vector<int> idx(n);

  for (int r = 0; r < m; ++r) {
    // distances
    if (dist_code == 2) {
      for (int i = 0; i < n; ++i) {
        double mismatches = 0.0;
        for (int j = 0; j < p; ++j)
          if (rpm_x(i, j) != Xtest(r, j)) mismatches += 1.0;
        d[i] = mismatches / static_cast<double>(p);
      }
    } else if (active.empty()) {
      std::fill(d.begin(), d.end(), 0.0);
    } else if (dist_code == 1) {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a];
          acc += std::fabs((rpm_x(i, j) - Xtest(r, j)) * inv_range[a]);
        }
        d[i] = acc;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a];
          const double z = (rpm_x(i, j) - Xtest(r, j)) * inv_range[a];
          acc += z * z;
        }
        d[i] = std::sqrt(acc);
      }
    }

    const int kk = std::min(k, n);

    if (kk == 1) {
      // Aggregate all exact nearest-distance ties deterministically.
      double dmin = d[0];
      for (int i = 1; i < n; ++i) if (d[i] < dmin) dmin = d[i];
      std::vector<double> tied;
      for (int i = 0; i < n; ++i) if (d[i] == dmin) tied.push_back(yhat[i]);
      if (is_class) {
        std::vector<double> ones(tied.size(), 1.0);
        out[r] = weighted_mode(tied, ones);
      } else {
        // R's gravity() drops nonfinite values before aggregating.
        std::vector<double> finite_tied;
        finite_tied.reserve(tied.size());
        for (double v : tied) if (R_finite(v)) finite_tied.push_back(v);
        out[r] = gravity_value(finite_tied, false);
      }
      continue;
    }

    // stable ascending order: distance, then original index
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(),
                     [&d](int a, int b) { return d[a] < d[b]; });

    std::vector<double> dk(kk), yk(kk);
    for (int i = 0; i < kk; ++i) {
      dk[i] = d[idx[i]];
      yk[i] = yhat[idx[i]];
    }

    const std::vector<double> w = ensemble_weights(dk);
    if (is_class) {
      out[r] = weighted_mode(yk, w);
    } else {
      double dot = 0.0;
      for (int i = 0; i < kk; ++i) dot += yk[i] * w[i];
      out[r] = dot;
    }
  }

  return out;
}
