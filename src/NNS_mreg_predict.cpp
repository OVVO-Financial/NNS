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
#include <utility>
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

inline void validate_dist_code(int dist_code) {
  if (dist_code < 0 || dist_code > 3) {
    stop("Unknown NNS multivariate distance code");
  }
}

// A neighbour under the kernel's total order: distance, then original row.
struct Near { double d; int idx; };
inline bool near_less(const Near& a, const Near& b) {
  return a.d < b.d || (a.d == b.d && a.idx < b.idx);
}

// Largest sum of squares whose square root does not exceed s. sqrt() is
// monotone, so a running sum above it has a strictly larger distance.
inline double l2_raw_limit(double s) {
  if (!R_finite(s)) return R_PosInf;
  double c = s * s;
  if (!R_finite(c)) return R_PosInf;
  for (int it = 0; it < 64 && std::sqrt(c) > s; ++it) c = std::nextafter(c, 0.0);
  if (std::sqrt(c) > s) return R_PosInf;
  for (int it = 0; it < 64; ++it) {
    const double up = std::nextafter(c, R_PosInf);
    if (std::sqrt(up) <= s) c = up; else return c;
  }
  return R_PosInf;  // bound not established: no pruning
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
  const int kmax, dist_code, n, p, only_k;
  const bool is_class;
  const RankComponents& rank;
  std::vector<int> active;
  std::vector<double> col_range;
  // 1-D fast path (rank-one duplicated coordinate, e.g. dy.d_'s cbind(X*, X*)):
  // the range-normalised distance is non-decreasing in |coord - t|, so
  // neighbours are found in O(log n + k) per query.
  bool is_1d = false;
  int axis = -1;
  std::vector<double> sc;   // RPM coordinate sorted ascending
  std::vector<int> so;      // original indices for the sorted order
  // Neighbour search for the range-normalised metrics: the active coordinates
  // of every RPM row, stored row by row. With enough queries to repay the
  // sort, the rows are in ascending order of one coordinate (the walk axis),
  // so a query visits them outward from its own position.
  bool search = false;
  int pa = 0;                    // active columns
  int walk = -1;                 // walk axis as a position in `active`, or -1
  std::vector<double> rows;      // n x pa, row-major, in visiting order
  std::vector<int> row_id;       // original row of each stored row
  std::vector<double> axis_val;  // walk-axis coordinate of each stored row

  struct Search {
    std::vector<Near> best;
    std::vector<int> ties;
    std::vector<double> q;
  };

  PredictPathWorker(const NumericMatrix& rpm_x_, const NumericVector& yhat_,
                    const NumericMatrix& Xtest_, int kmax_, int dist_code_,
                    const NumericVector& mins_, const NumericVector& maxs_,
                    bool is_class_, const RankComponents& rank_, NumericMatrix& out_,
                    int only_k_ = 0)
      : rpm_x(rpm_x_), yhat(yhat_), Xtest(Xtest_), mins(mins_), maxs(maxs_),
        out(out_), kmax(kmax_), dist_code(dist_code_), n(rpm_x_.nrow()),
        p(rpm_x_.ncol()), only_k(only_k_), is_class(is_class_), rank(rank_) {
    for (int j = 0; j < p; ++j) {
      const double range = maxs[j] - mins[j];
      if (R_finite(range) && range > 0.0) { active.push_back(j); col_range.push_back(range); }
    }
    // Detect an effectively 1-D design: all active columns element-wise
    // identical in both rpm_x and Xtest (FACTOR/Hamming metric excluded).
    if (dist_code != 3 && !active.empty()) {
      axis = active[0];
      bool identical = true;
      const int mq = Xtest.nrow();
      for (std::size_t a = 1; a < active.size() && identical; ++a) {
        const int j = active[a];
        for (int i = 0; i < n; ++i) if (rpm_x(i, j) != rpm_x(i, axis)) { identical = false; break; }
        for (int i = 0; identical && i < mq; ++i) if (Xtest(i, j) != Xtest(i, axis)) { identical = false; break; }
      }
      if (identical) {
        is_1d = true;
        std::vector<double> coordv(n);
        for (int i = 0; i < n; ++i) coordv[i] = rpm_x(i, axis);
        so.resize(n);
        for (int i = 0; i < n; ++i) so[i] = i;
        std::sort(so.begin(), so.end(), [&coordv](int a, int b) {
          return coordv[a] < coordv[b] || (coordv[a] == coordv[b] && a < b);
        });
        sc.resize(n);
        for (int i = 0; i < n; ++i) sc[i] = coordv[so[i]];
      }
    }
    if (!is_1d && dist_code != 3 && !active.empty()) {
      search = true;
      pa = static_cast<int>(active.size());
      std::vector<int> ord(n);
      for (int i = 0; i < n; ++i) ord[i] = i;
      if (Xtest.nrow() >= 4 && n >= 64) {
        // Walk along the coordinate most spread out relative to its range:
        // its term rules out the most rows.
        double best_spread = 0.0;
        for (int a = 0; a < pa; ++a) {
          const int j = active[a];
          bool finite = true;
          double mu = 0.0;
          for (int i = 0; i < n && finite; ++i) {
            finite = R_finite(rpm_x(i, j));
            mu += rpm_x(i, j) / n;
          }
          if (!finite) continue;
          double spread = 0.0;
          for (int i = 0; i < n; ++i) {
            const double u = (rpm_x(i, j) - mu) / col_range[a];
            spread += u * u;
          }
          if (spread > best_spread) { best_spread = spread; walk = a; }
        }
        if (walk >= 0) {
          const int j = active[walk];
          std::sort(ord.begin(), ord.end(), [this, j](int x, int y) {
            const double vx = rpm_x(x, j), vy = rpm_x(y, j);
            return vx < vy || (vx == vy && x < y);
          });
        }
      }
      rows.resize(static_cast<std::size_t>(n) * pa);
      for (int s = 0; s < n; ++s) {
        double* x = rows.data() + static_cast<std::size_t>(s) * pa;
        for (int a = 0; a < pa; ++a) x[a] = rpm_x(ord[s], active[a]);
      }
      row_id.swap(ord);
      if (walk >= 0) {
        axis_val.resize(n);
        for (int s = 0; s < n; ++s) axis_val[s] = rows[static_cast<std::size_t>(s) * pa + walk];
      }
    }
  }

  // Query r's `keep` nearest RPM rows under (distance, row), nearest first in
  // S.best, and, when `ties` is set, every row at the smallest distance in
  // S.ties, ascending. Each visited row's distance is accumulated over the
  // columns exactly as distance() does. Candidates collect in a buffer that
  // is cut back to the `keep` best whenever it doubles; the worst of those
  // bounds every later row. A row strictly beyond that bound cannot enter,
  // and along the walk axis neither can any row further out than `reach`:
  // that coordinate's own term already exceeds the bound, and the terms only
  // add up.
  void neighbours(std::size_t r, int keep, bool ties, Search& S) const {
    const int P = pa;
    double* q = S.q.data();
    for (int a = 0; a < P; ++a) q[a] = Xtest(r, active[a]);
    const double* range = col_range.data();
    const double* data = rows.data();
    std::vector<Near>& best = S.best;
    best.clear();
    S.ties.clear();
    const bool all = keep >= n;  // every row is kept: nothing to rule out
    const std::size_t cap = static_cast<std::size_t>(keep) * 2 + 32;
    bool seen = false;
    double limit = R_PosInf;     // a sum above this cannot matter
    double bound = R_PosInf;     // the distance `limit` was derived from
    double reach = R_PosInf;     // walk-axis gap beyond which no row matters
    double dmin = R_PosInf;
    const double rg = walk >= 0 ? range[walk] : 1.0;

    auto set_bound = [&](double b) {
      if (b == bound) return;
      bound = b;
      limit = dist_code == 1 ? l2_raw_limit(b) : b;
      if (walk >= 0 && limit < R_PosInf) {
        // The largest normalised gap whose own term stays within the limit,
        // widened well past rounding; the floor keeps gaps whose term would
        // underflow.
        const double zmax = dist_code == 0 ? (std::sqrt(1.0 + 4.0 * limit) - 1.0) / 2.0
                          : dist_code == 2 ? limit : std::sqrt(limit);
        reach = std::max(zmax * rg * (1.0 + 1e-9), rg * 1e-150);
      }
    };
    auto cut = [&]() {
      std::nth_element(best.begin(), best.begin() + (keep - 1), best.end(), near_less);
      best.resize(keep);
    };
    auto visit = [&](int s) {
      const double* x = data + static_cast<std::size_t>(s) * P;
      double acc = 0.0;
      if (dist_code == 0) {
        for (int a = 0; a < P; ++a) {
          const double z = (x[a] - q[a]) / range[a];
          acc += std::fabs(z) + z * z;
        }
      } else if (dist_code == 2) {
        for (int a = 0; a < P; ++a) {
          const double z = (x[a] - q[a]) / range[a];
          acc += std::fabs(z);
        }
      } else {
        for (int a = 0; a < P; ++a) {
          const double z = (x[a] - q[a]) / range[a];
          acc += z * z;
        }
      }
      if (acc > limit) return;
      const double d = dist_code == 1 ? std::sqrt(acc) : acc;
      const int id = row_id[s];
      if (ties) {
        if (!seen || d < dmin) { dmin = d; S.ties.clear(); seen = true; }
        if (d == dmin) S.ties.push_back(id);
        // With ties alone, only rows at the smallest distance matter.
        if (keep == 0) set_bound(dmin);
      }
      if (keep > 0) {
        best.push_back(Near{d, id});
        if (!all && best.size() >= cap) {
          cut();
          set_bound(best[keep - 1].d);
        }
      }
    };

    if (walk < 0 || (keep > 0 && cap >= static_cast<std::size_t>(n))) {
      // No walk axis, or a buffer that never fills: every row is visited.
      for (int s = 0; s < n; ++s) visit(s);
    } else {
      const double t = q[walk];
      int hi = static_cast<int>(std::lower_bound(axis_val.begin(), axis_val.end(), t) -
                                axis_val.begin());
      int lo = hi - 1;
      // The rows nearest along the axis first, until they set a bound.
      while (!(limit < R_PosInf) && (lo >= 0 || hi < n)) {
        const bool go_left = hi >= n || (lo >= 0 && (t - axis_val[lo]) <= (axis_val[hi] - t));
        if (go_left) visit(lo--); else visit(hi++);
      }
      for (; lo >= 0 && (t - axis_val[lo]) <= reach; --lo) visit(lo);
      for (; hi < n && (axis_val[hi] - t) <= reach; ++hi) visit(hi);
    }
    if (keep > 0) {
      if (!all && static_cast<int>(best.size()) > keep) cut();
      std::sort(best.begin(), best.end(), near_less);
    }
    if (ties) std::sort(S.ties.begin(), S.ties.end());
  }

  // The k = 1 estimate from the rows tied at the smallest distance.
  double ties_value(const std::vector<int>& ids, std::vector<double>& tied,
                    std::vector<double>& weights) const {
    tied.clear();
    for (std::size_t i = 0; i < ids.size(); ++i)
      if (R_finite(yhat[ids[i]])) tied.push_back(yhat[ids[i]]);
    if (is_class) {
      weights.assign(tied.size(), 1.0);
      return weighted_mode(tied, weights);
    }
    return gravity_value(tied, false);
  }

  // Range-normalised distance of training row i to query r. The 1-D path uses
  // it too, so it selects neighbours and ties exactly as the generic path does.
  double distance(int i, std::size_t r) const {
    double acc = 0.0;
    for (std::size_t a = 0; a < active.size(); ++a) {
      const int j = active[a];
      const double z = (rpm_x(i, j) - Xtest(r, j)) / col_range[a];
      if (dist_code == 2) acc += std::fabs(z);
      else if (dist_code == 0) acc += std::fabs(z) + z * z;
      else acc += z * z;
    }
    return dist_code == 1 ? std::sqrt(acc) : acc;
  }

  // Walk outward in |coord - t| order, then extend over distance() ties at the boundary.
  void topk_1d(std::size_t row, double t, int k, std::vector<int>& out_idx,
               std::vector<double>& out_d) const {
    const int nn = static_cast<int>(sc.size());
    k = std::min(k, nn);
    int r = static_cast<int>(std::lower_bound(sc.begin(), sc.end(), t) - sc.begin());
    int lo = r - 1;
    std::vector< std::pair<double, int> > win;
    while (static_cast<int>(win.size()) < k && (lo >= 0 || r < nn)) {
      const double dl = (lo >= 0) ? std::fabs(sc[lo] - t) : R_PosInf;
      const double dr = (r < nn) ? std::fabs(sc[r] - t) : R_PosInf;
      const int pos = (dl <= dr) ? lo-- : r++;
      win.emplace_back(distance(so[pos], row), so[pos]);
    }
    double thr = 0.0;
    for (const auto& w : win) thr = std::max(thr, w.first);
    for (double dd; lo >= 0 && (dd = distance(so[lo], row)) <= thr; --lo) win.emplace_back(dd, so[lo]);
    for (double dd; r < nn && (dd = distance(so[r], row)) <= thr; ++r) win.emplace_back(dd, so[r]);
    std::sort(win.begin(), win.end());  // (distance, original index)
    out_idx.clear();
    for (int i = 0; i < k; ++i) { out_idx.push_back(win[i].second); out_d[i] = win[i].first; }
  }

  double min_ties_1d(std::size_t row, double t, std::vector<double>& tied, std::vector<double>& weights) const {
    const int nn = static_cast<int>(sc.size());
    const int r = static_cast<int>(std::lower_bound(sc.begin(), sc.end(), t) - sc.begin());
    const int lo = r - 1;
    const double dl = (lo >= 0) ? distance(so[lo], row) : R_PosInf;
    const double dr = (r < nn) ? distance(so[r], row) : R_PosInf;
    const double mn = std::min(dl, dr);
    tied.clear();
    for (int li = lo; li >= 0 && distance(so[li], row) == mn; --li)
      if (R_finite(yhat[so[li]])) tied.push_back(yhat[so[li]]);
    for (int ri = r; ri < nn && distance(so[ri], row) == mn; ++ri)
      if (R_finite(yhat[so[ri]])) tied.push_back(yhat[so[ri]]);
    if (is_class) { weights.assign(tied.size(), 1.0); return weighted_mode(tied, weights); }
    return gravity_value(tied, false);
  }

  void distances_for_row(std::size_t r, std::vector<double>& d) const {
    if (dist_code == 3) {
      for (int i = 0; i < n; ++i) {
        double mismatches = 0.0;
        for (int j = 0; j < p; ++j) if (rpm_x(i, j) != Xtest(r, j)) mismatches += 1.0;
        d[i] = mismatches / static_cast<double>(p);
      }
    } else if (active.empty()) {
      std::fill(d.begin(), d.end(), 0.0);
    } else if (dist_code <= 2) {
      for (int i = 0; i < n; ++i) d[i] = distance(i, r);
    } else {
      stop("Unknown NNS multivariate distance code");
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
    // The full distance vector is only needed without a neighbour search
    // (Hamming distances, or no active column).
    const bool scan = !is_1d && !search;
    std::vector<double> d(scan ? n : 0), dk(kmax), yk(kmax), tied, one_weights;
    std::vector<int> idx(scan ? n : 0), tk;
    Search S;
    if (search) {
      S.q.resize(pa);
      S.best.reserve(std::min(static_cast<std::size_t>(kmax) * 2 + 33, static_cast<std::size_t>(n) + 1));
    }
    std::vector<double> student(kmax), inverse(kmax), normal(kmax), rbf(kmax), total(kmax);
    std::vector<double> cls(kmax), cls_totals(kmax);
    // Total-order comparator (distance, then original index) so partial
    // selection reproduces std::stable_sort with `d[a] < d[b]` byte-for-byte.
    auto cmp = [&d](int a, int b) { return d[a] < d[b] || (d[a] == d[b] && a < b); };
    for (std::size_t r = begin; r < end; ++r) {
      if (is_1d) {
        // Rank-one duplicated coordinate: O(log n + k) neighbour search.
        const double t = Xtest(r, axis);
        if (only_k > 0) {
          if (only_k == 1) {
            out(r, 0) = min_ties_1d(r, t, tied, one_weights);
          } else {
            topk_1d(r, t, only_k, tk, dk);
            for (int i = 0; i < only_k; ++i) yk[i] = yhat[tk[i]];
            out(r, 0) = predict_sorted_k_noalloc(dk, yk, only_k, is_class, rank,
                                                 student, inverse, normal, rbf,
                                                 total, cls, cls_totals);
          }
        } else {
          out(r, 0) = min_ties_1d(r, t, tied, one_weights);
          if (kmax >= 2) {
            topk_1d(r, t, kmax, tk, dk);
            for (int i = 0; i < kmax; ++i) yk[i] = yhat[tk[i]];
            for (int k = 2; k <= kmax; ++k) {
              out(r, k - 1) = predict_sorted_k_noalloc(dk, yk, k, is_class, rank,
                                                       student, inverse, normal, rbf,
                                                       total, cls, cls_totals);
            }
          }
        }
        continue;
      }
      if (search) {
        if (only_k == 1 || kmax == 1) {
          neighbours(r, 0, true, S);
          out(r, 0) = ties_value(S.ties, tied, one_weights);
        } else if (only_k > 0) {
          neighbours(r, only_k, false, S);
          for (int i = 0; i < only_k; ++i) { dk[i] = S.best[i].d; yk[i] = yhat[S.best[i].idx]; }
          out(r, 0) = predict_sorted_k_noalloc(dk, yk, only_k, is_class, rank,
                                               student, inverse, normal, rbf,
                                               total, cls, cls_totals);
        } else {
          neighbours(r, kmax, true, S);
          out(r, 0) = ties_value(S.ties, tied, one_weights);
          for (int i = 0; i < kmax; ++i) { dk[i] = S.best[i].d; yk[i] = yhat[S.best[i].idx]; }
          for (int k = 2; k <= kmax; ++k) {
            out(r, k - 1) = predict_sorted_k_noalloc(dk, yk, k, is_class, rank,
                                                     student, inverse, normal, rbf,
                                                     total, cls, cls_totals);
          }
        }
        continue;
      }
      distances_for_row(r, d);
      if (only_k == 1 || kmax == 1) {
        // One neighbour is every row tied at the smallest distance.
        out(r, 0) = min_ties(d, tied, one_weights);
        continue;
      }
      for (int i = 0; i < n; ++i) idx[i] = i;
      // Retain only the kmax nearest rows: O(n + kmax log kmax) instead of a
      // full O(n log n) sort.
      if (kmax < n) {
        std::nth_element(idx.begin(), idx.begin() + kmax, idx.end(), cmp);
        std::sort(idx.begin(), idx.begin() + kmax, cmp);
      } else {
        std::sort(idx.begin(), idx.end(), cmp);
      }
      for (int i = 0; i < kmax; ++i) { dk[i] = d[idx[i]]; yk[i] = yhat[idx[i]]; }
      if (only_k > 0) {
        // Single selected n.best: emit just the k-th prediction, no 1..k path.
        out(r, 0) = (only_k == 1)
            ? min_ties(d, tied, one_weights)
            : predict_sorted_k_noalloc(dk, yk, only_k, is_class, rank,
                                       student, inverse, normal, rbf,
                                       total, cls, cls_totals);
      } else {
        out(r, 0) = min_ties(d, tied, one_weights);
        for (int k = 2; k <= kmax; ++k) {
          out(r, k - 1) = predict_sorted_k_noalloc(dk, yk, k, is_class, rank,
                                                   student, inverse, normal, rbf,
                                                   total, cls, cls_totals);
        }
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
  validate_dist_code(dist_code);
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
  const int n = rpm_x.nrow(), p = rpm_x.ncol(), m = Xtest.nrow();
  if (yhat.size() != n) stop("yhat length must equal nrow(rpm_x)");
  if (Xtest.ncol() != p) stop("Xtest and rpm_x must have the same columns");
  if (mins.size() != p || maxs.size() != p) stop("mins/maxs must have one value per column");
  if (k < 1) stop("k must be >= 1");
  validate_dist_code(dist_code);
  if (k > n) k = n;
  // Genuine single-k kernel: select the k nearest rows and emit only the k-th
  // prediction, instead of computing the whole 1..k path and discarding it.
  RankComponents rank(k);
  NumericMatrix out(m, 1);
  PredictPathWorker worker(rpm_x, yhat, Xtest, k, dist_code, mins, maxs,
                           is_class, rank, out, /*only_k=*/k);
  const int threads = std::max(1, nthreads);
  if (threads == 1 || m < 2) worker(0, m); else parallelFor(0, m, worker, 1, threads);
  return out(_, 0);
}

// dist_code: 0 = native NNS, 1 = L2, 2 = L1, 3 = FACTOR (Hamming over encoded columns).
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
  validate_dist_code(dist_code);

  // Active columns and reciprocal ranges for the normalized metrics.
  std::vector<int> active;
  std::vector<double> col_range;
  for (int j = 0; j < p; ++j) {
    const double range = maxs[j] - mins[j];
    if (R_finite(range) && range > 0.0) {
      active.push_back(j);
      col_range.push_back(range);
    }
  }

  NumericVector out(m);
  std::vector<double> d(n);
  std::vector<int> idx(n);

  for (int r = 0; r < m; ++r) {
    // distances
    if (dist_code == 3) {
      for (int i = 0; i < n; ++i) {
        double mismatches = 0.0;
        for (int j = 0; j < p; ++j)
          if (rpm_x(i, j) != Xtest(r, j)) mismatches += 1.0;
        d[i] = mismatches / static_cast<double>(p);
      }
    } else if (active.empty()) {
      std::fill(d.begin(), d.end(), 0.0);
    } else if (dist_code == 2) {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a];
          acc += std::fabs((rpm_x(i, j) - Xtest(r, j)) / col_range[a]);
        }
        d[i] = acc;
      }
    } else if (dist_code == 0) {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a];
          const double z = (rpm_x(i, j) - Xtest(r, j)) / col_range[a];
          acc += std::fabs(z) + z * z;
        }
        d[i] = acc;
      }
    } else if (dist_code == 1) {
      for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t a = 0; a < active.size(); ++a) {
          const int j = active[a];
          const double z = (rpm_x(i, j) - Xtest(r, j)) / col_range[a];
          acc += z * z;
        }
        d[i] = std::sqrt(acc);
      }
    } else {
      stop("Unknown NNS multivariate distance code");
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

    // Retain only the kk nearest rows in stable (distance, then original
    // index) order via partial selection - O(n + kk log kk) per query.
    for (int i = 0; i < n; ++i) idx[i] = i;
    auto cmp = [&d](int a, int b) { return d[a] < d[b] || (d[a] == d[b] && a < b); };
    if (kk < n) {
      std::nth_element(idx.begin(), idx.begin() + kk, idx.end(), cmp);
      std::sort(idx.begin(), idx.begin() + kk, cmp);
    } else {
      std::sort(idx.begin(), idx.end(), cmp);
    }

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
