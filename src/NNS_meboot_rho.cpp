// The rho-targeting loop of NNS.meboot in one compiled call.
//
// For every replicate, R/NNS_meboot.R sorted the replicate's residuals, laid
// them out aligned (e) and anti-aligned (m) with the original residuals' ranks
// and ran optimize() over the weight t of t * m + (1 - t) * e, each objective
// value ranking a length-n vector in R. This kernel does the same arithmetic:
//   - rank(ties.method = "average") as R's do_rank() assigns it;
//   - mean() and sum() accumulate in long double as R's do, mean() with R's
//     second pass;
//   - the minimizer is R's own Brent_fmin (nns_fmin.c).
// Each replicate is bit-identical to the R loop. A replicate is left for the
// R loop (and flagged) when optimize() would warn there: a value that is not
// finite, or a candidate with no spread.

// Keep a * b + c as R evaluates it, unfused, on every compiler.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <utility>
#include <vector>

using namespace Rcpp;

extern "C" double nns_brent_fmin(double ax, double bx, double (*f)(double, void*),
                                 void* info, double tol);

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

struct Replicate {
  int n;
  bool spearman;
  const double* e;         // aligned with the original residuals' ranks
  const double* m;         // anti-aligned
  const double* target;    // centred target
  double target_norm;
  double rho;
  bool finite;             // every objective value was finite
  std::vector<double> comb, rank, centred, product;
  std::vector<std::pair<double, int> > order;

  // rank(comb, ties.method = "average"): equal values share the mean of the
  // positions they occupy, so their order within a tie does not matter.
  void average_ranks() {
    for (int i = 0; i < n; ++i) order[i] = std::make_pair(comb[i], i);
    std::sort(order.begin(), order.end());
    for (int i = 0; i < n;) {
      int j = i;
      while (j < n - 1 && order[j].first == order[j + 1].first) ++j;
      const double shared = (i + j + 2) / 2.;
      for (int k = i; k <= j; ++k) rank[order[k].second] = shared;
      i = j + 1;
    }
  }

  // abs(fast_corr(t * m + (1 - t) * e) - rho)
  double objective(double t) {
    const double s = 1 - t;
    for (int i = 0; i < n; ++i) {
      const double left = t * m[i];
      const double right = s * e[i];
      comb[i] = left + right;
    }
    const double* xv = comb.data();
    if (spearman) {
      average_ranks();
      xv = rank.data();
    }
    const double mu = r_mean(xv, n);
    for (int i = 0; i < n; ++i) centred[i] = xv[i] - mu;
    for (int i = 0; i < n; ++i) product[i] = centred[i] * centred[i];
    const double xn = std::sqrt(r_sum(product.data(), n));
    if (xn == 0) return NA_REAL;
    for (int i = 0; i < n; ++i) product[i] = centred[i] * target[i];
    const double scale = xn * target_norm;
    const double corr = r_sum(product.data(), n) / scale;
    if (!R_FINITE(corr)) return NA_REAL;
    return std::fabs(corr - rho);
  }
};

double replicate_objective(double t, void* info) {
  Replicate* rep = static_cast<Replicate*>(info);
  const double value = rep->objective(t);
  if (!R_FINITE(value)) {
    // optimize() warns and substitutes here; the R loop does that.
    rep->finite = false;
    return DBL_MAX;
  }
  return value;
}

}  // namespace

// The replicates after rho targeting, and for each whether the R loop must
// redo it. aligned / anti index the sorted replicate (1-based).
// [[Rcpp::export]]
List NNS_meboot_rho_cpp(const NumericMatrix& replicates, const IntegerVector& aligned,
                        const IntegerVector& anti, const NumericVector& target_centered,
                        double target_norm, const NumericVector& rho, bool spearman) {
  const int n = replicates.nrow(), reps = replicates.ncol();
  if (aligned.size() != n || anti.size() != n || target_centered.size() != n)
    stop("the rank indices and the target need one value per observation");
  if (rho.size() != reps) stop("one rho per replicate");
  for (int i = 0; i < n; ++i) {
    if (aligned[i] < 1 || aligned[i] > n || anti[i] < 1 || anti[i] > n)
      stop("rank indices must lie between 1 and the number of observations");
  }

  NumericMatrix out = clone(replicates);
  LogicalVector redo(reps);
  std::vector<double> sorted(n), e(n), m(n);
  Replicate rep;
  rep.n = n;
  rep.spearman = spearman;
  rep.target = target_centered.begin();
  rep.target_norm = target_norm;
  rep.comb.resize(n); rep.rank.resize(n); rep.centred.resize(n); rep.product.resize(n);
  rep.order.resize(n);

  for (int c = 0; c < reps; ++c) {
    if (c % 64 == 0) checkUserInterrupt();
    const double* column = &replicates(0, c);
    bool finite = R_FINITE(rho[c]);
    for (int i = 0; i < n && finite; ++i) finite = R_FINITE(column[i]);
    if (!finite) { redo[c] = TRUE; continue; }

    // sort(res_mat[, i])
    std::copy(column, column + n, sorted.begin());
    std::stable_sort(sorted.begin(), sorted.end());
    for (int i = 0; i < n; ++i) {
      e[i] = sorted[aligned[i] - 1];
      m[i] = sorted[anti[i] - 1];
    }
    rep.e = e.data();
    rep.m = m.data();
    rep.rho = rho[c];
    rep.finite = true;

    // optimize(obj, interval = c(0, 1), tol = 0.01)$minimum
    const double t_hat = nns_brent_fmin(0.0, 1.0, replicate_objective, &rep, 0.01);
    if (!rep.finite) { redo[c] = TRUE; continue; }

    const double s_hat = 1 - t_hat;
    for (int i = 0; i < n; ++i) {
      const double left = t_hat * m[i];
      const double right = s_hat * e[i];
      out(i, c) = left + right;
    }
  }
  return List::create(_["replicates"] = out, _["redo"] = redo);
}

// apply(x, 2, sort) for a matrix of finite values, or NULL when a value is
// not finite and sort() would drop it.
// [[Rcpp::export]]
SEXP NNS_sort_columns_cpp(const NumericMatrix& x) {
  const int n = x.nrow(), p = x.ncol();
  const double* in = x.begin();
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(n) * p; ++i) {
    if (!R_FINITE(in[i])) return R_NilValue;
  }
  NumericMatrix out(n, p);
  for (int j = 0; j < p; ++j) {
    double* column = &out(0, j);
    std::copy(in + static_cast<R_xlen_t>(j) * n, in + static_cast<R_xlen_t>(j + 1) * n, column);
    std::stable_sort(column, column + n);
  }
  return out;
}
