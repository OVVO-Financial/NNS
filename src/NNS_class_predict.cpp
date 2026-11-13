// The class-probability prediction loop of R/Classification.R in one call.
//
// For every point the R code measured the distance to every cell centre,
// ordered the cells, and for each neighbour count built the ensemble weights
// and the weighted class shares, all in interpreted R. This kernel does the
// same arithmetic in the same order:
//   - sum(), mean(), var(), colMeans() and rowSums() accumulate in long
//     double exactly as R's own C code does (var() and mean() with R's second
//     pass), and `^` is R_pow();
//   - the densities are R's own (Rmath);
//   - the weighted shares are the dgemv call crossprod() makes, so they match
//     on whichever BLAS R is running.
// The result is bit-identical to the R loop it replaces.

// Keep a * b + c as R evaluates it, unfused, on every compiler.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rcpp.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <map>
#include <vector>

using namespace Rcpp;

namespace {

// sum(x)
inline double r_sum(const double* x, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  if (s > DBL_MAX) return R_PosInf;
  if (s < -DBL_MAX) return R_NegInf;
  return static_cast<double>(s);
}

// mean(x)
inline double r_mean(const double* x, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  const bool finite_s = R_FINITE(static_cast<double>(s));
  if (finite_s) {
    s /= n;
  } else {
    s = 0.0L;
    for (int i = 0; i < n; ++i) s += x[i] / n;
  }
  if (finite_s && R_FINITE(static_cast<double>(s))) {
    long double t = 0.0L;
    for (int i = 0; i < n; ++i) t += (x[i] - s);
    s += t / n;
  } else if (R_FINITE(static_cast<double>(s))) {
    long double t = 0.0L;
    for (int i = 0; i < n; ++i) t += (x[i] - s) / n;
    s += t;
  }
  return static_cast<double>(s);
}

// var(x) of a vector with at least two values
inline double r_var(const double* x, int n) {
  for (int k = 0; k < n; ++k) if (ISNAN(x[k])) return NA_REAL;
  long double sum = 0.0L;
  for (int k = 0; k < n; ++k) sum += x[k];
  long double tmp = sum / n;
  if (R_FINITE(static_cast<double>(tmp))) {
    sum = 0.0L;
    for (int k = 0; k < n; ++k) sum += (x[k] - tmp);
    tmp = tmp + sum / n;
  }
  const long double xm = static_cast<double>(tmp);
  sum = 0.0L;
  for (int k = 0; k < n; ++k) sum += (x[k] - xm) * (x[k] - xm);
  return static_cast<double>(sum / (n - 1));
}

// .nns_mreg_normalize_weights()
inline void normalize(double* w, int k) {
  for (int i = 0; i < k; ++i) if (!R_FINITE(w[i]) || w[i] < 0.0) w[i] = 0.0;
  const double s = r_sum(w, k);
  if (s <= 0.0) {
    const double u = 1.0 / k;
    for (int i = 0; i < k; ++i) w[i] = u;
  } else {
    for (int i = 0; i < k; ++i) w[i] /= s;
  }
}

// .nns_mreg_rank_weights(k), k >= 2
struct RankTerms {
  double uniform;
  std::vector<double> exponential, lognormal, power;
};

void rank_terms(int k, RankTerms& out) {
  std::vector<double> ranks(k);
  for (int i = 0; i < k; ++i) ranks[i] = i + 1;
  const double rank_sd = std::sqrt(r_var(ranks.data(), k));
  out.uniform = 1.0 / k;
  out.exponential.resize(k);
  const double rate = 1.0 / k;        // dexp(ranks, rate = 1 / k)
  const double scale = 1.0 / rate;    // which R passes on as 1 / rate
  for (int i = 0; i < k; ++i) out.exponential[i] = ::Rf_dexp(ranks[i], scale, 0);
  normalize(out.exponential.data(), k);
  out.lognormal.assign(k, 0.0);
  if (R_FINITE(rank_sd) && rank_sd > 0.0) {
    for (int i = 0; i < k; ++i)
      out.lognormal[i] = std::fabs(::Rf_dlnorm(ranks[i], 0.0, rank_sd, 1));
    normalize(out.lognormal.data(), k);
    std::reverse(out.lognormal.begin(), out.lognormal.end());
  }
  out.power.resize(k);
  for (int i = 0; i < k; ++i) out.power[i] = R_pow(ranks[i], -2.0);
  normalize(out.power.data(), k);
}

struct Scratch {
  std::vector<double> student, inverse, normal, rbf, weights, rows, z;
  void size(int m, int K) {
    if (static_cast<int>(student.size()) < m) {
      student.resize(m); inverse.resize(m); normal.resize(m); rbf.resize(m);
      weights.resize(m);
    }
    if (rows.size() < static_cast<std::size_t>(m) * K) rows.resize(static_cast<std::size_t>(m) * K);
    if (static_cast<int>(z.size()) < K) z.resize(K);
  }
};

// .nns_mreg_ensemble_weights(near, rank), then equal distances share their
// mean weight. `near` is ascending.
void ensemble(const double* near, int m, const RankTerms* rank, Scratch& S) {
  double* w = S.weights.data();
  if (m == 1) { w[0] = 1.0; return; }
  double* student = S.student.data();
  double* inverse = S.inverse.data();
  double* normal = S.normal.data();
  double* rbf = S.rbf.data();
  for (int i = 0; i < m; ++i) student[i] = ::Rf_dt(near[i], static_cast<double>(m), 0);
  normalize(student, m);
  for (int i = 0; i < m; ++i) inverse[i] = 1.0 / std::max(near[i], 1e-12);
  normalize(inverse, m);
  const double var = r_var(near, m);
  const double sd = std::sqrt(var);
  if (R_FINITE(sd) && sd > 0.0) {
    for (int i = 0; i < m; ++i) normal[i] = ::Rf_dnorm4(near[i], 0.0, sd, 0);
    normalize(normal, m);
  } else {
    for (int i = 0; i < m; ++i) normal[i] = 0.0;
  }
  if (R_FINITE(var) && var > 0.0) {
    const double width = 2.0 * var;
    for (int i = 0; i < m; ++i) rbf[i] = std::exp(-near[i] / width);
    normalize(rbf, m);
  } else {
    for (int i = 0; i < m; ++i) rbf[i] = 0.0;
  }
  for (int i = 0; i < m; ++i) {
    w[i] = rank->uniform + student[i] + inverse[i] + rank->exponential[i] +
      rank->lognormal[i] + rank->power[i] + normal[i] + rbf[i];
  }
  normalize(w, m);
  for (int i = 0; i < m;) {
    int j = i + 1;
    while (j < m && near[j] == near[i]) ++j;
    if (j - i > 1) {
      const double shared = r_mean(w + i, j - i);
      for (int l = i; l < j; ++l) w[l] = shared;
    }
    i = j;
  }
}

// crossprod(weights, rows): rows is m x K, column-major.
void weighted_shares(Scratch& S, int m, int K) {
  double* w = S.weights.data();
  double* y = S.rows.data();
  double* z = S.z.data();
  bool finite = true;
  for (int i = 0; i < m && finite; ++i) finite = R_FINITE(w[i]);
  for (std::size_t i = 0; i < static_cast<std::size_t>(m) * K && finite; ++i) finite = R_FINITE(y[i]);
  if (!finite) {
    // R keeps values that are not finite out of the BLAS.
    for (int k = 0; k < K; ++k) {
      double sum = 0.0;
      for (int j = 0; j < m; ++j) sum += w[j] * y[j + static_cast<std::size_t>(k) * m];
      z[k] = sum;
    }
    return;
  }
  const double one = 1.0, zero = 0.0;
  const int ione = 1;
  if (K == 1) {
    F77_CALL(dgemv)("T", &m, &ione, &one, w, &m, y, &ione, &zero, z, &ione FCONE);
  } else {
    F77_CALL(dgemv)("T", &m, &K, &one, y, &m, w, &ione, &zero, z, &ione FCONE);
  }
}

}  // namespace

// One matrix of unnormalised class shares (points x classes) per neighbour
// count in ks. dist_code: 0 = NNS, 1 = L2, 2 = L1, 3 = FACTOR.
// [[Rcpp::export]]
List NNS_class_predict_cpp(const NumericMatrix& centres, const NumericMatrix& shares,
                           const NumericMatrix& points, const IntegerVector& ks,
                           int dist_code, const NumericVector& minimums,
                           const NumericVector& maximums) {
  const int n = centres.nrow(), p = centres.ncol(), K = shares.ncol();
  const int m_points = points.nrow(), nk = ks.size();
  if (shares.nrow() != n) stop("shares and centres must have the same rows");
  if (points.ncol() != p) stop("points and centres must have the same columns");
  if (minimums.size() != p || maximums.size() != p) stop("one minimum and maximum per column");
  if (dist_code < 0 || dist_code > 3) stop("Unknown NNS multivariate distance code");
  if (n < 1) stop("no cells");
  int k_max = 0;
  for (int a = 0; a < nk; ++a) {
    if (ks[a] < 1 || ks[a] > n) stop("neighbour counts must lie between 1 and the number of cells");
    k_max = std::max(k_max, ks[a]);
  }

  List out(nk);
  std::vector<double*> out_ptr(nk);
  for (int a = 0; a < nk; ++a) {
    NumericMatrix m(m_points, K);
    out[a] = m;
    out_ptr[a] = REAL(m);
  }
  if (m_points == 0 || nk == 0) return out;

  // Range-scaled coordinates use the columns with a positive finite range;
  // the Hamming distance uses them all.
  std::vector<int> active;
  std::vector<double> range;
  if (dist_code == 3) {
    for (int j = 0; j < p; ++j) active.push_back(j);
  } else {
    for (int j = 0; j < p; ++j) {
      const double r = maximums[j] - minimums[j];
      if (R_FINITE(r) && r > 0.0) { active.push_back(j); range.push_back(r); }
    }
  }
  const int pa = static_cast<int>(active.size());
  std::vector<double> cell(static_cast<std::size_t>(n) * pa);
  for (int c = 0; c < n; ++c)
    for (int a = 0; a < pa; ++a) cell[static_cast<std::size_t>(c) * pa + a] = centres(c, active[a]);

  std::vector<double> d(n), q(pa), ds(n), exact_share(K);
  std::vector<int> idx(n), exact;
  std::map<int, RankTerms> ranks;
  Scratch S;
  auto nearer = [&d](int a, int b) { return d[a] < d[b] || (d[a] == d[b] && a < b); };

  for (int i = 0; i < m_points; ++i) {
    if (i % 256 == 0) checkUserInterrupt();
    for (int a = 0; a < pa; ++a) q[a] = points(i, active[a]);

    // Distances, each row summed in long double as rowSums() and rowMeans() do.
    for (int c = 0; c < n; ++c) {
      const double* x = cell.data() + static_cast<std::size_t>(c) * pa;
      long double acc = 0.0L;
      if (dist_code == 0) {
        for (int a = 0; a < pa; ++a) {
          const double z = (x[a] - q[a]) / range[a];
          const double term = std::fabs(z) + z * z;
          acc += term;
        }
        d[c] = static_cast<double>(acc);
      } else if (dist_code == 2) {
        for (int a = 0; a < pa; ++a) {
          const double z = (x[a] - q[a]) / range[a];
          acc += std::fabs(z);
        }
        d[c] = static_cast<double>(acc);
      } else if (dist_code == 1) {
        for (int a = 0; a < pa; ++a) {
          const double z = (x[a] - q[a]) / range[a];
          const double term = z * z;
          acc += term;
        }
        d[c] = std::sqrt(static_cast<double>(acc));
      } else {
        for (int a = 0; a < pa; ++a) acc += (x[a] != q[a]) ? 1.0 : 0.0;
        acc /= p;
        d[c] = static_cast<double>(acc);
      }
    }

    // Exact matches keep their own class frequencies.
    exact.clear();
    for (int c = 0; c < n; ++c) if (d[c] == 0.0) exact.push_back(c);
    if (!exact.empty()) {
      const int cnt = static_cast<int>(exact.size());
      for (int k = 0; k < K; ++k) {
        long double sum = 0.0L;
        for (int e = 0; e < cnt; ++e) sum += shares(exact[e], k);
        sum /= cnt;
        exact_share[k] = static_cast<double>(sum);
      }
      for (int a = 0; a < nk; ++a)
        for (int k = 0; k < K; ++k) out_ptr[a][i + static_cast<std::size_t>(k) * m_points] = exact_share[k];
      continue;
    }

    // The cells in order(d), as far as the largest count and its ties need.
    for (int c = 0; c < n; ++c) idx[c] = c;
    int have = n;
    if (k_max < n) {
      std::nth_element(idx.begin(), idx.begin() + (k_max - 1), idx.end(), nearer);
      const double edge = d[idx[k_max - 1]];
      have = static_cast<int>(std::partition(idx.begin() + k_max, idx.end(),
                                             [&](int c) { return d[c] <= edge; }) - idx.begin());
    }
    std::sort(idx.begin(), idx.begin() + have, nearer);
    for (int l = 0; l < have; ++l) ds[l] = d[idx[l]];

    for (int a = 0; a < nk; ++a) {
      // The k nearest cells and every cell tied at the k-th distance.
      const int m = static_cast<int>(std::upper_bound(ds.begin(), ds.begin() + have, ds[ks[a] - 1]) -
                                     ds.begin());
      S.size(m, K);
      const RankTerms* rank = NULL;
      if (m > 1) {
        std::map<int, RankTerms>::iterator it = ranks.find(m);
        if (it == ranks.end()) {
          it = ranks.insert(std::make_pair(m, RankTerms())).first;
          rank_terms(m, it->second);
        }
        rank = &it->second;
      }
      ensemble(ds.data(), m, rank, S);
      for (int k = 0; k < K; ++k)
        for (int l = 0; l < m; ++l) S.rows[l + static_cast<std::size_t>(k) * m] = shares(idx[l], k);
      weighted_shares(S, m, K);
      for (int k = 0; k < K; ++k) out_ptr[a][i + static_cast<std::size_t>(k) * m_points] = S.z[k];
    }
  }
  return out;
}

// colMeans() of each group's rows, the groups numbered from 1 in `codes`:
// one long double sum per column over the group's rows in order, as R does.
// [[Rcpp::export]]
NumericMatrix NNS_group_means_cpp(const NumericMatrix& X, const IntegerVector& codes,
                                  int groups) {
  const int n = X.nrow(), p = X.ncol();
  if (codes.size() != n) stop("one group code per row");
  if (groups < 0) stop("the number of groups cannot be negative");
  std::vector<int> count(groups, 0);
  for (int i = 0; i < n; ++i) {
    if (codes[i] < 1 || codes[i] > groups) stop("group codes must lie between 1 and the number of groups");
    ++count[codes[i] - 1];
  }
  NumericMatrix out(groups, p);
  std::vector<long double> sum(groups);
  for (int j = 0; j < p; ++j) {
    std::fill(sum.begin(), sum.end(), 0.0L);
    for (int i = 0; i < n; ++i) sum[codes[i] - 1] += X(i, j);
    for (int g = 0; g < groups; ++g) {
      sum[g] /= count[g];
      out(g, j) = static_cast<double>(sum[g]);
    }
  }
  return out;
}
