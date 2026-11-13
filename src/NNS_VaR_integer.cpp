// LPM.VaR / UPM.VaR for integer degrees 1 to 4 in one compiled call.
//
// R/LPM_UPM_VaR.R's .NNS_LPM_VaR_integer sorted the data, built prefix power
// sums and the ratio at every distinct value, located each percentile's
// interval, and then solved one polynomial per percentile through
// stats::uniroot() with an R closure. This kernel is the same computation:
//   - the power sums accumulate in long double as cumsum() and sum() do, and
//     powers go through R_pow() as R's `^` does;
//   - each polynomial term is evaluated in R's order, unfused;
//   - the root search is R's own zeroin (nns_zeroin.c).
// The result is bit-identical to the R implementation, which is kept as
// .NNS_LPM_VaR_integer_R and takes over whenever this kernel declines: a
// value that is not finite, a percentile that is missing, or a search that
// would warn or stop in uniroot().

// Keep a * b + c as R evaluates it, unfused, on every compiler.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

using namespace Rcpp;

extern "C" double nns_zeroin2(double ax, double bx, double fa, double fb,
                              double (*f)(double x, void* info), void* info,
                              double* Tol, int* Maxit);

namespace {

// One interval's polynomial: the prefix and suffix power sums either side of
// it, and the coefficient choose(degree, j) * (-1)^j of each term.
struct Interval {
  int degree;
  double percentile;
  const double* prefix;     // degree + 1 prefix power sums
  double suffix[5];         // total - prefix
  double lower_sign[5];     // choose(degree, j) * (-1)^j
  double upper_sign[5];     // choose(degree, j) * (-1)^(degree - j)
  bool finite;              // every value of f was finite

  // .NNS_LPM_raw_integer
  double lpm(double t) const {
    double value = 0.0;
    for (int j = 0; j <= degree; ++j) {
      const double term = lower_sign[j] * R_pow(t, static_cast<double>(degree - j)) * prefix[j];
      value = value + term;
    }
    return value;
  }
  // .NNS_UPM_raw_integer
  double upm(double t) const {
    double value = 0.0;
    for (int j = 0; j <= degree; ++j) {
      const double term = upper_sign[j] * R_pow(t, static_cast<double>(degree - j)) * suffix[j];
      value = value + term;
    }
    return value;
  }
  // .NNS_pm_root_value_integer
  double f(double t) const {
    const double l = lpm(t), u = upm(t);
    const double left = (1 - percentile) * l;
    const double right = percentile * u;
    return left - right;
  }
  double ratio_gap(double t) const {
    const double l = lpm(t), u = upm(t);
    return l / (l + u) - percentile;
  }
};

double root_value(double t, void* info) {
  Interval* cell = static_cast<Interval*>(info);
  const double value = cell->f(t);
  if (!R_FINITE(value)) {
    // uniroot() warns and substitutes here; the R implementation does that.
    cell->finite = false;
    return value == R_NegInf ? -DBL_MAX : DBL_MAX;
  }
  return value;
}

}  // namespace

// LPM.VaR(percentile, degree, x) for degree 1 to 4 and x without missing
// values, or NULL when the R implementation has to handle the call.
// [[Rcpp::export]]
SEXP NNS_LPM_VaR_integer_cpp(const NumericVector& percentile, int degree,
                             const NumericVector& x) {
  const int n = x.size(), np = percentile.size();
  if (degree < 1 || degree > 4 || n < 1) return R_NilValue;
  for (int i = 0; i < n; ++i) if (!R_FINITE(x[i])) return R_NilValue;
  for (int i = 0; i < np; ++i) if (ISNAN(percentile[i])) return R_NilValue;

  // sort(x), centred at its median order statistic.
  std::vector<double> xs(x.begin(), x.end());
  std::stable_sort(xs.begin(), xs.end());
  const double shift = xs[(n + 1) / 2 - 1];
  for (int i = 0; i < n; ++i) xs[i] = xs[i] - shift;
  const double x_min = xs[0], x_max = xs[n - 1];
  NumericVector out(np);
  if (n == 1 || x_min == x_max) {
    const double value = x_min + shift;
    for (int i = 0; i < np; ++i) out[i] = value;
    return out;
  }

  // rle(): each run's last value and the count up to it.
  std::vector<double> unique_x;
  std::vector<int> k_break;
  for (int i = 0; i < n; ++i) {
    if (i == n - 1 || xs[i + 1] != xs[i]) {
      unique_x.push_back(xs[i]);
      k_break.push_back(i + 1);
    }
  }
  const int U = static_cast<int>(unique_x.size());

  // Prefix power sums at the breaks and their totals, as cumsum() and sum().
  const int rows = degree + 1;
  std::vector<double> prefix(static_cast<std::size_t>(rows) * U);
  double total[5];
  for (int u = 0; u < U; ++u) prefix[static_cast<std::size_t>(u) * rows] = k_break[u];
  total[0] = n;
  for (int j = 1; j <= degree; ++j) {
    long double sum = 0.0L;
    int u = 0;
    for (int i = 0; i < n; ++i) {
      sum += R_pow(xs[i], static_cast<double>(j));
      if (i + 1 == k_break[u]) {
        prefix[static_cast<std::size_t>(u) * rows + j] = static_cast<double>(sum);
        ++u;
      }
    }
    total[j] = sum > DBL_MAX ? R_PosInf : sum < -DBL_MAX ? R_NegInf : static_cast<double>(sum);
  }

  double lower_sign[5], upper_sign[5];
  for (int j = 0; j <= degree; ++j) {
    const double c = ::Rf_choose(static_cast<double>(degree), static_cast<double>(j));
    lower_sign[j] = c * R_pow(-1.0, static_cast<double>(j));
    upper_sign[j] = c * R_pow(-1.0, static_cast<double>(degree - j));
  }

  // .NNS_LPM_ratio_at_breaks, clamped to [0, 1] and made monotone.
  std::vector<double> ratio(U);
  for (int u = 0; u < U; ++u) {
    const double t = unique_x[u];
    const double* pre = prefix.data() + static_cast<std::size_t>(u) * rows;
    double lpm = 0.0, upm = 0.0;
    for (int j = 0; j <= degree; ++j) {
      const double power = R_pow(t, static_cast<double>(degree - j));
      const double below = lower_sign[j] * power * pre[j];
      lpm = lpm + below;
      const double suffix = total[j] - pre[j];
      const double above = upper_sign[j] * power * suffix;
      upm = upm + above;
    }
    double r = lpm / (lpm + upm);
    if (ISNAN(r)) r = 0.0;
    if (0.0 > r) r = 0.0;
    if (1.0 < r) r = 1.0;
    ratio[u] = r;
  }
  for (int u = 1; u < U; ++u) if (ratio[u - 1] > ratio[u]) ratio[u] = ratio[u - 1];
  ratio[0] = 0.0;
  ratio[U - 1] = 1.0;

  const double gap_tol = R_pow(DBL_EPSILON, 0.75);
  const double eps2 = DBL_EPSILON * DBL_EPSILON;

  for (int i = 0; i < np; ++i) {
    const double p = percentile[i];
    if (p <= 0.0) { out[i] = x_min + shift; continue; }
    if (p >= 1.0) { out[i] = x_max + shift; continue; }

    // findInterval(p, ratio, rightmost.closed = TRUE), kept inside the breaks.
    int interval = static_cast<int>(std::upper_bound(ratio.begin(), ratio.end(), p) - ratio.begin());
    if (interval < 1) interval = 1;
    if (interval > U - 1) interval = U - 1;

    const double lower = unique_x[interval - 1], upper = unique_x[interval];
    if (lower == upper) { out[i] = lower + shift; continue; }

    Interval cell;
    cell.degree = degree;
    cell.percentile = p;
    cell.prefix = prefix.data() + static_cast<std::size_t>(interval - 1) * rows;
    cell.finite = true;
    for (int j = 0; j <= degree; ++j) {
      cell.suffix[j] = total[j] - cell.prefix[j];
      cell.lower_sign[j] = lower_sign[j];
      cell.upper_sign[j] = upper_sign[j];
    }

    const double f_lower = cell.f(lower), f_upper = cell.f(upper);
    const double g_lower = cell.ratio_gap(lower), g_upper = cell.ratio_gap(upper);

    double value;
    if (R_FINITE(g_lower) && std::fabs(g_lower) <= gap_tol) {
      value = lower;
    } else if (R_FINITE(g_upper) && std::fabs(g_upper) <= gap_tol) {
      value = upper;
    } else if (R_FINITE(f_lower) && R_FINITE(f_upper) && f_lower * f_upper <= 0.0) {
      // uniroot(f, c(lower, upper), tol = .Machine$double.eps^2 * (upper - lower))
      double tol = eps2 * (upper - lower);
      int iterations = 1000;
      if (!R_FINITE(tol) || tol <= 0.0) return R_NilValue;
      const double root = nns_zeroin2(lower, upper, f_lower, f_upper, root_value, &cell,
                                      &tol, &iterations);
      if (!cell.finite || iterations < 0) return R_NilValue;
      value = root;
      if (lower > value) value = lower;   // pmin(pmax(root, lower), upper)
      if (upper < value) value = upper;
    } else {
      // No sign change, from rounding at a break: the closer endpoint.
      if (ISNAN(g_lower) || ISNAN(g_upper)) return R_NilValue;
      value = std::fabs(g_lower) <= std::fabs(g_upper) ? lower : upper;
    }
    out[i] = value + shift;
  }
  return out;
}
