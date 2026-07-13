// Native kernel for the repaired NNS.M.reg prediction rule.
// Replicates R/Multivariate_Regression.R's .nns_mreg_predict_one exactly:
//   - range-normalized L2 / L1 or Hamming (FACTOR) distances,
//   - stable ascending ordering with index tie-break,
//   - k = 1 aggregates all exact nearest-distance ties (gravity / class mode),
//   - k > 1 uses the eight-component ensemble weights built from Rmath
//     densities, matching the R helper term for term,
//   - exact double-precision weighted class mode (max total, tie -> min value).
// Row-parallel with a per-call thread count; no global RcppParallel thread state is touched.
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;
using namespace RcppParallel;

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


inline void rank_components(int k, std::vector<double>& uniform,
                            std::vector<double>& exponential,
                            std::vector<double>& lognormal,
                            std::vector<double>& power) {
  uniform.assign(k, 1.0 / static_cast<double>(k));
  exponential.resize(k); lognormal.assign(k, 0.0); power.resize(k);
  const double kk = static_cast<double>(k);
  for (int i = 0; i < k; ++i) exponential[i] = ::Rf_dexp(static_cast<double>(i + 1), kk, 0);
  normalize_weights(exponential);
  std::vector<double> ranks(k);
  for (int i = 0; i < k; ++i) ranks[i] = static_cast<double>(i + 1);
  const double rank_sd = sample_sd(ranks);
  if (R_finite(rank_sd) && rank_sd > 0.0) {
    for (int i = 0; i < k; ++i) lognormal[i] = std::fabs(::Rf_dlnorm(ranks[i], 0.0, rank_sd, 1));
    normalize_weights(lognormal);
    std::reverse(lognormal.begin(), lognormal.end());
  }
  for (int i = 0; i < k; ++i) { const double r = static_cast<double>(i + 1); power[i] = 1.0 / (r * r); }
  normalize_weights(power);
}

struct RankComponents {
  std::vector< std::vector<double> > uniform, exponential, lognormal, power;
  explicit RankComponents(int kmax) : uniform(kmax + 1), exponential(kmax + 1),
    lognormal(kmax + 1), power(kmax + 1) {
    for (int k = 2; k <= kmax; ++k) rank_components(k, uniform[k], exponential[k], lognormal[k], power[k]);
  }
};

inline void normalize_prefix(std::vector<double>& values, int k) {
  double s = 0.0;
  for (int i = 0; i < k; ++i) {
    if (!R_finite(values[i]) || values[i] < 0.0) values[i] = 0.0;
    s += values[i];
  }
  if (s <= 0.0) {
    const double u = 1.0 / static_cast<double>(k);
    for (int i = 0; i < k; ++i) values[i] = u;
  } else {
    for (int i = 0; i < k; ++i) values[i] /= s;
  }
}

inline double weighted_mode_prefix(const std::vector<double>& yk,
                                   const std::vector<double>& weights,
                                   int k,
                                   std::vector<double>& cls,
                                   std::vector<double>& totals) {
  int ncls = 0;
  for (int i = 0; i < k; ++i) {
    const double v = yk[i], w = weights[i];
    if (!R_finite(v) || !R_finite(w) || w < 0.0) continue;
    int pos = -1;
    for (int c = 0; c < ncls; ++c) if (cls[c] == v) { pos = c; break; }
    if (pos < 0) { pos = ncls++; cls[pos] = v; totals[pos] = 0.0; }
    totals[pos] += w;
  }
  if (ncls == 0) return NA_REAL;
  double best_val = cls[0], best_tot = totals[0];
  for (int c = 1; c < ncls; ++c) {
    if (totals[c] > best_tot || (totals[c] == best_tot && cls[c] < best_val)) {
      best_tot = totals[c]; best_val = cls[c];
    }
  }
  return best_val;
}

inline double predict_sorted_k_noalloc(const std::vector<double>& dk,
                                       const std::vector<double>& yk,
                                       int k, bool is_class,
                                       const RankComponents& rank,
                                       std::vector<double>& student,
                                       std::vector<double>& inverse,
                                       std::vector<double>& normal,
                                       std::vector<double>& rbf,
                                       std::vector<double>& total,
                                       std::vector<double>& cls,
                                       std::vector<double>& cls_totals) {
  for (int i = 0; i < k; ++i) student[i] = ::Rf_dt(dk[i], static_cast<double>(k), 0);
  normalize_prefix(student, k);
  for (int i = 0; i < k; ++i) inverse[i] = 1.0 / std::max(dk[i], 1e-12);
  normalize_prefix(inverse, k);

  double mu = 0.0;
  for (int i = 0; i < k; ++i) mu += dk[i];
  mu /= static_cast<double>(k);
  double acc = 0.0;
  for (int i = 0; i < k; ++i) { const double z = dk[i] - mu; acc += z * z; }
  const double sd = (k > 1) ? std::sqrt(acc / static_cast<double>(k - 1)) : NA_REAL;
  const double var = R_finite(sd) ? sd * sd : NA_REAL;
  if (R_finite(sd) && sd > 0.0) {
    for (int i = 0; i < k; ++i) normal[i] = ::Rf_dnorm4(dk[i], 0.0, sd, 0);
    normalize_prefix(normal, k);
  } else {
    for (int i = 0; i < k; ++i) normal[i] = 0.0;
  }
  if (R_finite(var) && var > 0.0) {
    for (int i = 0; i < k; ++i) rbf[i] = std::exp(-dk[i] / (2.0 * var));
    normalize_prefix(rbf, k);
  } else {
    for (int i = 0; i < k; ++i) rbf[i] = 0.0;
  }

  for (int i = 0; i < k; ++i) {
    total[i] = rank.uniform[k][i] + student[i] + inverse[i] +
      rank.exponential[k][i] + rank.lognormal[k][i] + rank.power[k][i] +
      normal[i] + rbf[i];
  }
  normalize_prefix(total, k);

  if (is_class) return weighted_mode_prefix(yk, total, k, cls, cls_totals);
  double dot = 0.0;
  for (int i = 0; i < k; ++i) dot += yk[i] * total[i];
  return dot;
}

struct PredictPathWorker : public Worker {
  RMatrix<double> rpm_x;
  RVector<double> yhat;
  RMatrix<double> Xtest;
  RVector<double> mins;
  RVector<double> maxs;
  RMatrix<double> out;
  const int kmax, dist_code, n, p;
  const bool is_class;
  const RankComponents& rank;
  std::vector<int> active;
  std::vector<double> inv_range;

  PredictPathWorker(const NumericMatrix& rpm_x_, const NumericVector& yhat_,
                    const NumericMatrix& Xtest_, int kmax_, int dist_code_,
                    const NumericVector& mins_, const NumericVector& maxs_,
                    bool is_class_, const RankComponents& rank_, NumericMatrix& out_)
      : rpm_x(rpm_x_), yhat(yhat_), Xtest(Xtest_), mins(mins_), maxs(maxs_),
        out(out_), kmax(kmax_), dist_code(dist_code_), n(rpm_x_.nrow()),
        p(rpm_x_.ncol()), is_class(is_class_), rank(rank_) {
    for (int j = 0; j < p; ++j) {
      const double range = maxs[j] - mins[j];
      if (R_finite(range) && range > 0.0) { active.push_back(j); inv_range.push_back(1.0 / range); }
    }
  }

  void distances_for_row(std::size_t r, std::vector<double>& d) const {
    if (dist_code == 2) {
      for (int i = 0; i < n; ++i) {
        double mismatches = 0.0;
        for (int j = 0; j < p; ++j) if (rpm_x(i, j) != Xtest(r, j)) mismatches += 1.0;
        d[i] = mismatches / static_cast<double>(p);
      }
    } else if (active.empty()) {
      std::fill(d.begin(), d.end(), 0.0);
    } else if (dist_code == 1) {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a]; acc += std::fabs((rpm_x(i, j) - Xtest(r, j)) * inv_range[a]);
        }
        d[i] = acc;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a]; const double z = (rpm_x(i, j) - Xtest(r, j)) * inv_range[a]; acc += z * z;
        }
        d[i] = std::sqrt(acc);
      }
    }
  }

  double min_ties(const std::vector<double>& d, std::vector<double>& tied,
                  std::vector<double>& weights) const {
    double dmin = d[0];
    for (int i = 1; i < n; ++i) if (d[i] < dmin) dmin = d[i];
    tied.clear();
    for (int i = 0; i < n; ++i) if (d[i] == dmin && R_finite(yhat[i])) tied.push_back(yhat[i]);
    if (is_class) {
      weights.assign(tied.size(), 1.0);
      return weighted_mode(tied, weights);
    }
    return gravity_value(tied, false);
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> d(n), dk(kmax), yk(kmax), tied, one_weights;
    std::vector<int> idx(n);
    std::vector<double> student(kmax), inverse(kmax), normal(kmax), rbf(kmax), total(kmax);
    std::vector<double> cls(kmax), cls_totals(kmax);
    for (std::size_t r = begin; r < end; ++r) {
      distances_for_row(r, d);
      for (int i = 0; i < n; ++i) idx[i] = i;
      std::stable_sort(idx.begin(), idx.end(), [&d](int a, int b) { return d[a] < d[b]; });
      for (int i = 0; i < kmax; ++i) { dk[i] = d[idx[i]]; yk[i] = yhat[idx[i]]; }
      out(r, 0) = min_ties(d, tied, one_weights);
      for (int k = 2; k <= kmax; ++k) {
        out(r, k - 1) = predict_sorted_k_noalloc(dk, yk, k, is_class, rank,
                                                 student, inverse, normal, rbf,
                                                 total, cls, cls_totals);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix NNS_mreg_predict_path_v2_cpp(const NumericMatrix& rpm_x,
                                           const NumericVector& yhat,
                                           const NumericMatrix& Xtest,
                                           int kmax,
                                           int dist_code,
                                           const NumericVector& mins,
                                           const NumericVector& maxs,
                                           bool is_class,
                                           int nthreads) {
  const int n = rpm_x.nrow(), p = rpm_x.ncol(), m = Xtest.nrow();
  if (yhat.size() != n) stop("yhat length must equal nrow(rpm_x)");
  if (Xtest.ncol() != p) stop("Xtest and rpm_x must have the same columns");
  if (mins.size() != p || maxs.size() != p) stop("mins/maxs must have one value per column");
  if (kmax < 1) stop("kmax must be >= 1");
  if (kmax > n) kmax = n;
  RankComponents rank(kmax);
  NumericMatrix out(m, kmax);
  PredictPathWorker worker(rpm_x, yhat, Xtest, kmax, dist_code, mins, maxs,
                           is_class, rank, out);
  const int threads = std::max(1, nthreads);
  if (threads == 1 || m < 2) worker(0, m); else parallelFor(0, m, worker, 1, threads);
  return out;
}

// [[Rcpp::export]]
NumericMatrix NNS_mreg_predict_path_cpp(const NumericMatrix& rpm_x,
                                        const NumericVector& yhat,
                                        const NumericMatrix& Xtest,
                                        int kmax,
                                        int dist_code,
                                        const NumericVector& mins,
                                        const NumericVector& maxs,
                                        bool is_class) {
  return NNS_mreg_predict_path_v2_cpp(rpm_x, yhat, Xtest, kmax, dist_code, mins, maxs, is_class, 1);
}

// [[Rcpp::export]]
NumericVector NNS_mreg_predict_v2_cpp(const NumericMatrix& rpm_x,
                                      const NumericVector& yhat,
                                      const NumericMatrix& Xtest,
                                      int k,
                                      int dist_code,
                                      const NumericVector& mins,
                                      const NumericVector& maxs,
                                      bool is_class,
                                      int nthreads) {
  NumericMatrix path = NNS_mreg_predict_path_v2_cpp(rpm_x, yhat, Xtest, k, dist_code, mins, maxs, is_class, nthreads);
  return path(_, k - 1);
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
