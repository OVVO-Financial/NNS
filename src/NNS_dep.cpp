// NNS_dep.cpp
// C++ implementation of NNS.dep and NNS.dep.matrix.
//
// Exported functions (called from R):
//   NNS_dep_pair_cpp    - bivariate dependence given pre-computed partition labels
// (NNS.dep's matrix form loops over NNS_dep_pair_cpp in R/Dependence.R.)
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <cfloat>

using namespace Rcpp;
using namespace RcppParallel;

// ============================================================
// INTERNAL HELPERS
// ============================================================

struct DepResult {
  double correlation;
  double dependence;
};

// ------------------------------------------------------------
// Per-cell dependence with the independence floor removed (#70).
//
// Inside a partition cell of m points, around the cell means, two concordance
// terms are measured:
//   T0 = |2 p0 - 1|, p0 = share of points in the co-lower or co-upper orthant
//        (degree 0: frequency)
//   T1 = |sum dx dy| / sum |dx dy|  (degree 1: area)
// (the discordant terms are exactly 1 - p0 and 1 - p1, so they add nothing).
// Each term is measured against its own expectation under the within-cell
// permutation null,
//   T* = (T - E0[T]) / (1 - E0[T]),
// which keeps perfect dependence at 1 and centres independence near 0. The
// cell's dependence is (T0* + T1*) / 2: no square root, which would turn the
// sampling noise of an unbiased estimate of 0 back into a floor of order
// n^(-1/4).
// E0[T0] uses the hypergeometric distribution (a normal approximation for
// cells larger than 400). E0[T1] is exact whenever the distinct
// arrangements of the cell can be enumerated cheaply (every cell of up to 6
// points, and larger cells with heavily tied values) and uses the permutation
// variance of sum dx dy otherwise.
// ------------------------------------------------------------

// A plain aggregate, so `return {a, b};` compiles as C++11 too: the package
// sets no C++ standard, and older R (4.2 on Windows, for one) builds as C++11.
struct CellCopula {
  double signed_value;  // debiased dependence times the sign of the cell covariance
  double magnitude;     // debiased dependence; may be negative below chance
};

// E0|2 p0 - 1| for a cell of m points with a points at or below the x mean and
// b points at or below the y mean, under random pairing: the co-lower count K
// is hypergeometric and the concordant count is m - a - b + 2K.
static double null_T0(int m, int a, int b) {
  const int k_lo = std::max(0, a + b - m);
  const int k_hi = std::min(a, b);
  if (k_lo >= k_hi) {  // a single possible count: the term carries no information
    return 1.0;
  }
  const double md = static_cast<double>(m);
  if (m > 400) {
    // Normal approximation with the exact hypergeometric mean and variance.
    const double ad = a, bd = b;
    const double mean_k = ad * bd / md;
    const double var_k = bd * (ad / md) * ((md - ad) / md) * ((md - bd) / (md - 1.0));
    const double mu = (2.0 * (md - ad - bd + 2.0 * mean_k) - md) / md;
    const double sd = 4.0 * std::sqrt(std::max(var_k, 0.0)) / md;
    if (sd <= 0.0) return std::abs(mu);
    const double z = mu / sd;
    return sd * std::sqrt(2.0 / M_PI) * std::exp(-0.5 * z * z) +
      mu * (1.0 - 2.0 * R::pnorm(-z, 0.0, 1.0, 1, 0));
  }
  // Hypergeometric probabilities by the ratio recurrence from the lowest count.
  double pk = R::dhyper(static_cast<double>(k_lo), static_cast<double>(a),
                        static_cast<double>(m - a), static_cast<double>(b), 0);
  double e = 0.0;
  for (int k = k_lo; k <= k_hi; ++k) {
    const double conc = static_cast<double>(m - a - b + 2 * k);
    e += pk * std::abs(2.0 * conc / md - 1.0);
    if (k < k_hi) {
      pk *= (static_cast<double>(a - k) * static_cast<double>(b - k)) /
            (static_cast<double>(k + 1) * static_cast<double>(m - a - b + k + 1));
    }
  }
  return e;
}

// log of the number of distinct arrangements of v: m! / prod(multiplicity!).
static double log_arrangements(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const int m = static_cast<int>(v.size());
  double out = std::lgamma(static_cast<double>(m) + 1.0);
  for (int i = 0; i < m;) {
    int j = i;
    while (j < m && v[j] == v[i]) ++j;
    out -= std::lgamma(static_cast<double>(j - i) + 1.0);
    i = j;
  }
  return out;
}

// Enumerate up to this many distinct arrangements exactly (6! = 720, so every
// cell of up to 6 points is exact, and so are larger cells with heavy ties,
// such as 7 points of which 6 are tied).
static const double kLogExactArrangements = std::log(720.0) + 1e-9;
static const int kExactCellSize = 6;

// E0[ |sum dx dy| / sum |dx dy| ] under random pairing of dx and dy. The
// arrangement counts are only worth computing when a variable has ties.
// legacy = true keeps the rule NNS.dep used before #95 (exact only for cells of
// at most 6 points), for the local measure behind NNS.reg's default order.
static double null_T1(const std::vector<double>& dx, const std::vector<double>& dy,
                      bool x_ties, bool y_ties, bool legacy = false) {
  const int m = static_cast<int>(dx.size());
  // Pairing dy with a random permutation of dx is the same null as the reverse,
  // so permute whichever has fewer distinct arrangements. Every distinct
  // arrangement of a multiset is equally likely under a uniform permutation.
  double lx = 1e300, ly = 1e300;
  if (m <= kExactCellSize || (x_ties && !legacy)) lx = log_arrangements(dx);
  if (m <= kExactCellSize || (y_ties && !legacy)) ly = log_arrangements(dy);
  if (std::min(lx, ly) <= kLogExactArrangements) {
    std::vector<double> v = (lx <= ly) ? dx : dy;
    const std::vector<double>& w = (lx <= ly) ? dy : dx;
    std::sort(v.begin(), v.end());
    double total = 0.0;
    double count = 0.0;
    do {
      double s = 0.0, d = 0.0;
      for (int i = 0; i < m; ++i) {
        const double p = w[i] * v[i];
        s += p;
        d += std::abs(p);
      }
      total += d > 0.0 ? std::abs(s) / d : 0.0;
      count += 1.0;
    } while (std::next_permutation(v.begin(), v.end()));
    return total / count;
  }
  double sx2 = 0.0, sy2 = 0.0, sax = 0.0, say = 0.0;
  for (int i = 0; i < m; ++i) {
    sx2 += dx[i] * dx[i];
    sy2 += dy[i] * dy[i];
    sax += std::abs(dx[i]);
    say += std::abs(dy[i]);
  }
  if (sax <= 0.0 || say <= 0.0) return 1.0;
  // sum dx dy has mean 0 and variance sx2 * sy2 / (m - 1) under permutation;
  // sum |dx dy| has mean sax * say / m.
  const double md = static_cast<double>(m);
  const double e = std::sqrt(2.0 / M_PI) * std::sqrt(sx2 * sy2 / (md - 1.0)) * md / (sax * say);
  return std::min(1.0, e);
}

static inline double debias_term(double t, double e0) {
  if (!(e0 < 1.0 - 1e-12)) return 0.0;
  return (t - e0) / (1.0 - e0);
}

// Constancy of a cell, for the rules applied by the partition (best_level).
enum CellKind { CELL_VARYING = 0, CELL_X_CONST = 1, CELL_BOTH_CONST = 2 };

static CellCopula cell_dependence(const std::vector<double>& xv,
                                  const std::vector<double>& yv,
                                  bool x_ties, bool y_ties, CellKind& kind,
                                  bool legacy_null = false) {
  const int n = static_cast<int>(xv.size());
  kind = CELL_VARYING;
  if (n < 2) return {0.0, 0.0};

  // Measure from the cell minimum, so deviations round at the cell's spread rather than its level.
  double rx = 0.0, ry = 0.0;
  if (!legacy_null) {
    rx = *std::min_element(xv.begin(), xv.end());
    ry = *std::min_element(yv.begin(), yv.end());
  }
  double tx = 0.0, ty = 0.0, ax = 0.0, ay = 0.0;
  for (int i = 0; i < n; ++i) {
    tx += xv[i] - rx; ty += yv[i] - ry;
    ax = std::max(ax, xv[i] - rx); ay = std::max(ay, yv[i] - ry);
  }
  tx /= n; ty /= n;
  // Bound on the mean's rounding, which depends on row order: a value equal to the mean can't fall either side.
  const double tie_x = legacy_null ? 0.0 : (n + 1) * DBL_EPSILON * ax;
  const double tie_y = legacy_null ? 0.0 : (n + 1) * DBL_EPSILON * ay;

  std::vector<double> dx(n), dy(n);
  int a = 0, b = 0;
  double conc = 0.0, s = 0.0, d = 0.0;
  for (int i = 0; i < n; ++i) {
    dx[i] = (xv[i] - rx) - tx;
    dy[i] = (yv[i] - ry) - ty;
    const bool lo_x = dx[i] <= tie_x, lo_y = dy[i] <= tie_y;
    a += lo_x;
    b += lo_y;
    // Concordant = co-lower (<= target) or co-upper (> target); the degree-0
    // counts partition the cell.
    if (lo_x == lo_y) conc += 1.0;
    const double p = dx[i] * dy[i];
    s += p;
    d += std::abs(p);
  }

  // x constant in the cell (discrete x): there is no within-cell variation to
  // test, only whether y is determined there. The partition scores that (see
  // best_level). Compare the raw values: a repeated 0.1 need not centre to
  // exactly 0.
  bool x_const = true, y_const = true;
  for (int i = 1; i < n; ++i) {
    if (xv[i] != xv[0]) x_const = false;
    if (yv[i] != yv[0]) y_const = false;
  }
  if (x_const) {
    kind = y_const ? CELL_BOTH_CONST : CELL_X_CONST;
    return {0.0, 0.0};
  }

  // s is summed in partition order, so a zero covariance can round either way.
  const double s_tol = legacy_null ? 0.0 : 1e-12 * d;
  const double slope_sign = (s > s_tol) ? 1.0 : (s < -s_tol) ? -1.0 : 0.0;
  const double e0 = null_T0(n, a, b);

  // All concordant or all discordant is perfect dependence, unless the cell is
  // too small for that to be informative (a single attainable count).
  // Compare counts, not count * (1/n): 49 * (1/49) != 1.
  if (e0 < 1.0 - 1e-12) {
    if (conc == static_cast<double>(n)) return {1.0, 1.0};
    if (conc == 0.0) return {-1.0, 1.0};
  }

  const double md = static_cast<double>(n);
  const double t0 = std::abs(2.0 * conc / md - 1.0);
  const double t1 = d > 0.0 ? std::abs(s) / d : 0.0;
  const double e1 = null_T1(dx, dy, x_ties, y_ties, legacy_null);

  const double dep = 0.5 * (debias_term(t0, e0) + debias_term(t1, e1));
  if (slope_sign == 0.0) return {0.0, 0.0};
  return {dep * slope_sign, dep};
}

struct LevelBest {
  double dependence;
  double correlation;
};

enum PartitionKind { CONDITIONAL, JOINT };

// Probability that m values drawn without replacement from v are all equal,
// i.e. that y is constant in a cell of m points under the permutation null.
// counts holds the multiplicities of v's distinct values, largest first.
static double p_constant(int m, const std::vector<int>& counts, int n,
                         std::vector<double>& cache) {
  if (cache[m] >= 0.0) return cache[m];
  const double lchoose_n = std::lgamma(n + 1.0) - std::lgamma(m + 1.0) - std::lgamma(n - m + 1.0);
  double p = 0.0;
  for (int c : counts) {
    if (c < m) break;
    p += std::exp(std::lgamma(c + 1.0) - std::lgamma(m + 1.0) - std::lgamma(c - m + 1.0) - lchoose_n);
  }
  cache[m] = std::min(1.0, p);
  return cache[m];
}

static std::vector<int> value_counts(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  std::vector<int> counts;
  for (size_t i = 0; i < v.size();) {
    size_t j = i;
    while (j < v.size() && v[j] == v[i]) ++j;
    counts.push_back(static_cast<int>(j - i));
    i = j;
  }
  std::sort(counts.begin(), counts.end(), std::greater<int>());
  return counts;
}

// A partition label: a view of the R string, which stays alive with the
// CharacterVector it was read from. (std::string_view needs C++17.)
// Comparisons are bytewise, as std::string_view's are.
struct Label {
  const char* p;
  std::size_t n;
  std::size_t size() const { return n; }
};

static inline bool operator<(const Label& a, const Label& b) {
  const int c = std::memcmp(a.p, b.p, std::min(a.n, b.n));
  return c < 0 || (c == 0 && a.n < b.n);
}

static inline bool operator==(const Label& a, const Label& b) {
  return a.n == b.n && std::memcmp(a.p, b.p, a.n) == 0;
}

// The first len characters of a and b agree (a shorter label counts whole).
static inline bool same_prefix(const Label& a, const Label& b, std::size_t len) {
  const std::size_t la = std::min(a.n, len), lb = std::min(b.n, len);
  return la == lb && std::memcmp(a.p, b.p, la) == 0;
}

// The rows in label order with ties in row order, as a stable sort by label
// gives them, and each distinct label with the start of its run.
struct LabelOrder {
  std::vector<int> ord;
  std::vector<Label> distinct;  // ascending
  std::vector<int> start;       // run starts in ord, then the row count
  std::size_t longest;          // longest label, at least 1
};

// A partition has few distinct labels, and equal R strings share one cached
// CHARSXP. So the rows are grouped by string pointer, only the distinct labels
// are compared as strings, and a counting sort lays the rows out: the order a
// stable sort of all the labels would give, without comparing every pair.
static LabelOrder order_labels(const CharacterVector& q) {
  const int n = static_cast<int>(q.size());
  LabelOrder out;
  std::vector<SEXP> seen;
  std::vector<int> code(n);
  // Open addressing on the pointer; the table stays at most half full.
  std::size_t mask = 63;
  std::vector<int> slot(mask + 1, -1);
  for (int i = 0; i < n; ++i) {
    const SEXP s = STRING_ELT(q, i);
    std::size_t h = (reinterpret_cast<std::uintptr_t>(s) >> 4) * static_cast<std::size_t>(0x9E3779B1u);
    h &= mask;
    while (slot[h] >= 0 && seen[slot[h]] != s) h = (h + 1) & mask;
    if (slot[h] < 0) {
      slot[h] = static_cast<int>(seen.size());
      seen.push_back(s);
      code[i] = slot[h];
      if (seen.size() * 2 > mask) {
        mask = mask * 2 + 1;
        slot.assign(mask + 1, -1);
        for (std::size_t d = 0; d < seen.size(); ++d) {
          std::size_t g = (reinterpret_cast<std::uintptr_t>(seen[d]) >> 4) * static_cast<std::size_t>(0x9E3779B1u);
          g &= mask;
          while (slot[g] >= 0) g = (g + 1) & mask;
          slot[g] = static_cast<int>(d);
        }
      }
    } else {
      code[i] = slot[h];
    }
  }

  const int found = static_cast<int>(seen.size());
  std::vector<Label> label(found);
  out.longest = 1;
  for (int d = 0; d < found; ++d) {
    label[d].p = CHAR(seen[d]);
    label[d].n = static_cast<std::size_t>(LENGTH(seen[d]));
    out.longest = std::max(out.longest, label[d].n);
  }
  std::vector<int> by(found);
  std::iota(by.begin(), by.end(), 0);
  std::sort(by.begin(), by.end(), [&](int a, int b) { return label[a] < label[b]; });
  // Two strings with the same bytes would be one label.
  std::vector<int> rank(found);
  for (int r = 0; r < found; ++r) {
    if (r == 0 || !(label[by[r]] == label[by[r - 1]])) out.distinct.push_back(label[by[r]]);
    rank[by[r]] = static_cast<int>(out.distinct.size()) - 1;
  }
  const int groups = static_cast<int>(out.distinct.size());
  out.start.assign(groups + 1, 0);
  for (int i = 0; i < n; ++i) ++out.start[rank[code[i]] + 1];
  for (int g = 0; g < groups; ++g) out.start[g + 1] += out.start[g];
  std::vector<int> next(out.start.begin(), out.start.end() - 1);
  out.ord.resize(n);
  for (int i = 0; i < n; ++i) out.ord[next[rank[code[i]]]++] = i;
  return out;
}

// Dependence of v on u over every level of a recursive mean-split partition
// whose cell labels are prefix codes ("q", "q1", "q12", ...): the cells at level
// l are the label prefixes of length l + 1. Each level's dependence is the
// cell-size weighted average of the cell dependences; the strongest level is
// returned (ties go to the finer level), with the signed average at that level
// as the correlation.
//
// Cells constant in u: in a conditional partition (on u alone) such a cell
// asks whether v is determined at that value of u. Its score is 1 when v is
// constant there and 0 otherwise, measured against p0, the chance that v is
// constant in a cell of that size under the permutation null:
// (score - p0) / (1 - p0). In the joint partition a constant cell is constant
// by construction, so it carries no information and scores 0.
static LevelBest best_level(const std::vector<double>& u,
                            const std::vector<double>& v,
                            const LabelOrder& lab,
                            double between_sign, PartitionKind kind,
                            bool u_ties, bool v_ties,
                            std::vector<double>* levels = nullptr) {
  const int n = static_cast<int>(u.size());
  const std::vector<int>& ord = lab.ord;
  const int groups = static_cast<int>(lab.distinct.size());
  const size_t L = lab.longest;

  std::vector<int> v_counts;
  std::vector<double> p0_cache;
  if (kind == CONDITIONAL && v_ties) {
    v_counts = value_counts(v);
    p0_cache.assign(n + 1, -1.0);
  }

  LevelBest best{-1e300, 0.0};
  int prev_cells = -1;
  std::vector<double> uq, vq;
  std::vector<int> starts;
  for (size_t len = 1; len <= L; ++len) {
    // Cells at this level are runs of equal prefixes in label order, so they
    // begin where a distinct label's prefix differs from the one before it.
    starts.clear();
    starts.push_back(0);
    for (int g = 1; g < groups; ++g) {
      if (!same_prefix(lab.distinct[g], lab.distinct[g - 1], len)) starts.push_back(lab.start[g]);
    }
    const int cells = static_cast<int>(starts.size());
    if (cells == prev_cells) continue;  // nested partition unchanged: same value
    prev_cells = cells;
    starts.push_back(n);

    double dep = 0.0, corr = 0.0;
    for (int c = 0; c < cells; ++c) {
      const int lo = starts[c], hi = starts[c + 1], m = hi - lo;
      uq.resize(m); vq.resize(m);
      for (int k = 0; k < m; ++k) { uq[k] = u[ord[lo + k]]; vq[k] = v[ord[lo + k]]; }
      CellKind ck;
      CellCopula cop = cell_dependence(uq, vq, u_ties, v_ties, ck);
      if (ck != CELL_VARYING) {
        double mag = 0.0;
        if (kind == CONDITIONAL) {
          const double p0 = v_ties ? p_constant(m, v_counts, n, p0_cache) : 0.0;
          const double score = (ck == CELL_BOTH_CONST) ? 1.0 : 0.0;
          if (p0 < 1.0 - 1e-12) mag = (score - p0) / (1.0 - p0);
        }
        // A constant cell has no within-cell sign; it takes the sign of the
        // overall covariance, read as 0 when that is within rounding of 0.
        cop = {between_sign * mag, mag};
      }
      if (!std::isfinite(cop.signed_value) || !std::isfinite(cop.magnitude)) cop = {0.0, 0.0};
      const double w = static_cast<double>(m) / static_cast<double>(n);
      dep += w * cop.magnitude;
      corr += w * cop.signed_value;
    }
    if (levels) levels->push_back(dep);
    if (dep >= best.dependence) best = {dep, corr};
  }
  // Below-chance levels can make every average negative; dependence is a
  // magnitude, so it is floored at 0, and |Correlation| is kept within it.
  best.dependence = std::max(0.0, best.dependence);
  best.correlation = std::min(best.dependence, std::max(-best.dependence, best.correlation));
  return best;
}

static bool has_ties(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return std::adjacent_find(v.begin(), v.end()) != v.end();
}

// x_ties / y_ties: whether the variable has repeated values, when the caller
// already knows (NA to have it found here).
// [[Rcpp::export]]
List NNS_dep_pair_cpp(NumericVector x,
                      NumericVector y,
                      CharacterVector quad_xy,
                      CharacterVector quad_yx,
                      CharacterVector quad_joint,
                      bool asym = false,
                      int x_ties_known = NA_INTEGER,
                      int y_ties_known = NA_INTEGER) {
  const int n = x.size();
  if (y.size() != n || quad_xy.size() != n) {
    stop("x, y and quad_xy must all have the same length");
  }
  if (!asym && (quad_yx.size() != n || quad_joint.size() != n)) {
    stop("quad_yx and quad_joint must have the length of x when asym = FALSE");
  }

  std::vector<double> xv(x.begin(), x.end());
  std::vector<double> yv(y.begin(), y.end());

  bool cx = true, cy = true;
  for (int i = 1; i < n; ++i) {
    if (xv[i] != xv[0]) cx = false;
    if (yv[i] != yv[0]) cy = false;
    if (!cx && !cy) break;
  }
  if (n < 2 || cx || cy) {
    return List::create(_["Correlation"] = 0.0, _["Dependence"] = 0.0);
  }

  // Sign of the overall covariance, accumulated from the first observation so
  // a large common level doesn't round away.
  const double x0 = xv[0], y0 = yv[0];
  double mx = 0.0, my = 0.0, cov = 0.0, cov_abs = 0.0;
  for (int i = 0; i < n; ++i) { mx += xv[i] - x0; my += yv[i] - y0; }
  mx /= n; my /= n;
  for (int i = 0; i < n; ++i) {
    const double p = ((xv[i] - x0) - mx) * ((yv[i] - y0) - my);
    cov += p;
    cov_abs += std::abs(p);
  }
  const double cov_tol = 1e-12 * cov_abs;
  const double between_sign = (cov > cov_tol) ? 1.0 : (cov < -cov_tol) ? -1.0 : 0.0;

  const bool x_ties = x_ties_known == NA_INTEGER ? has_ties(xv) : x_ties_known != 0;
  const bool y_ties = y_ties_known == NA_INTEGER ? has_ties(yv) : y_ties_known != 0;

  // y on x: the partition on x.
  LevelBest best = best_level(xv, yv, order_labels(quad_xy), between_sign, CONDITIONAL,
                              x_ties, y_ties);
  if (!asym) {
    // Symmetric: the strongest of the joint quadrants of (x, y), which carry
    // relations such as a circle or y = x * noise, and the two conditional
    // partitions, which carry functions of either variable. The joint
    // partition wins ties.
    LevelBest joint = best_level(xv, yv, order_labels(quad_joint), between_sign, JOINT,
                                 x_ties, y_ties);
    LevelBest on_y = best_level(yv, xv, order_labels(quad_yx), between_sign, CONDITIONAL,
                                y_ties, x_ties);
    LevelBest out = joint;
    if (best.dependence > out.dependence) out = best;
    if (on_y.dependence > out.dependence) out = on_y;
    best = out;
  }

  return List::create(_["Correlation"] = best.correlation,
                      _["Dependence"]  = best.dependence);
}

// Unfloored dependence of y on x at each distinct level of the partition on x
// (internal; NNS.caus compares the two directions on their level average).
// [[Rcpp::export]]
NumericVector NNS_dep_levels_cpp(NumericVector x, NumericVector y, CharacterVector quad_xy) {
  const int n = x.size();
  if (y.size() != n || quad_xy.size() != n) stop("x, y and quad_xy must all have the same length");
  std::vector<double> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  bool cx = true, cy = true;
  for (int i = 1; i < n; ++i) { if (xv[i] != xv[0]) cx = false; if (yv[i] != yv[0]) cy = false; }
  if (n < 2 || cx || cy) return NumericVector::create(0.0);
  std::vector<double> levels;
  best_level(xv, yv, order_labels(quad_xy), 0.0, CONDITIONAL, has_ties(xv), has_ties(yv), &levels);
  return wrap(levels);
}

// ------------------------------------------------------------
// Local (finest-cell) dependence: NNS.dep's measure before the multi-level
// redesign (#95). Each finest cell's debiased concordance enters as a signed
// square root, the cells are averaged and floored at 0, and discrete pairs
// are blended with the degree-0 copula. It is high only when structure
// survives inside the finest cells, which is how NNS.reg's default order
// (order ~ 10 x dependence) and dy.d_'s step rule read it: as a fine-scale
// signal-to-noise measure. NNS.dep reports the multi-level measure instead,
// which also counts dependence between cells (a noisy linear relationship)
// and would push those heuristics to over-deep partitions.
// ------------------------------------------------------------

static double copula_degree0_unsigned(const std::vector<double>& xv,
                                      const std::vector<double>& yv) {
  const int n = static_cast<int>(xv.size());
  if (n < 2) return 0.0;
  double tx = 0.0, ty = 0.0;
  for (int i = 0; i < n; ++i) { tx += xv[i]; ty += yv[i]; }
  tx /= n; ty /= n;
  double conc = 0.0, disc = 0.0;
  for (int i = 0; i < n; ++i) {
    const bool lo_x = (xv[i] - tx) <= 0.0, lo_y = (yv[i] - ty) <= 0.0;
    if (lo_x == lo_y) conc += 1.0; else disc += 1.0;
  }
  const double d0_Co = conc / n, dpm_d0 = disc / n;
  const double disc_dep = std::min(1.0, std::max(0.0, std::abs(d0_Co - 0.5) / 0.5));
  const double nd_disc = std::abs(dpm_d0 - 0.5) / 0.5;
  return std::sqrt((disc_dep + nd_disc) / 2.0);
}

static int count_unique(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return static_cast<int>(std::unique(v.begin(), v.end()) - v.begin());
}

static double local_dep_dir(const std::vector<double>& u, const std::vector<double>& v,
                            const LabelOrder& lab,
                            bool u_ties, bool v_ties) {
  const int n = static_cast<int>(u.size());
  auto root = [](double inner) { return (inner >= 0.0 ? 1.0 : -1.0) * std::sqrt(std::abs(inner)); };
  auto cell_value = [&](const std::vector<double>& uq, const std::vector<double>& vq) -> double {
    CellKind ck;
    CellCopula c = cell_dependence(uq, vq, u_ties, v_ties, ck, true);
    if (ck == CELL_BOTH_CONST) return 1.0;
    if (ck == CELL_X_CONST) return 0.0;
    return root(c.magnitude);
  };
  double global = cell_value(u, v);
  if (!std::isfinite(global)) global = 0.0;

  // Each finest cell is one distinct label's run of rows.
  const std::vector<int>& ord = lab.ord;
  const int groups = static_cast<int>(lab.distinct.size());
  double dep = 0.0;
  std::vector<double> uq, vq;
  for (int g = 0; g < groups; ++g) {
    const int lo = lab.start[g], hi = lab.start[g + 1];
    const int m = hi - lo;
    uq.resize(m); vq.resize(m);
    for (int k = 0; k < m; ++k) { uq[k] = u[ord[lo + k]]; vq[k] = v[ord[lo + k]]; }
    double val = cell_value(uq, vq);
    if (!std::isfinite(val)) val = global;
    dep += val * static_cast<double>(m) / static_cast<double>(n);
  }
  return std::max(0.0, dep);
}

// [[Rcpp::export]]
double NNS_dep_local_cpp(NumericVector x, NumericVector y,
                         CharacterVector quad_xy, CharacterVector quad_yx,
                         bool asym = true) {
  const int n = x.size();
  if (y.size() != n || quad_xy.size() != n) stop("x, y and quad_xy must all have the same length");
  if (!asym && quad_yx.size() != n) stop("quad_yx must have the length of x when asym = FALSE");
  std::vector<double> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  bool cx = true, cy = true;
  for (int i = 1; i < n; ++i) { if (xv[i] != xv[0]) cx = false; if (yv[i] != yv[0]) cy = false; }
  if (n < 2 || cx || cy) return 0.0;
  // One sort per variable gives its distinct count, and so whether it has ties.
  const int x_unique = count_unique(xv), y_unique = count_unique(yv);
  const bool x_ties = x_unique < n, y_ties = y_unique < n;
  double dep = local_dep_dir(xv, yv, order_labels(quad_xy), x_ties, y_ties);
  if (!asym) dep = std::max(dep, local_dep_dir(yv, xv, order_labels(quad_yx), y_ties, x_ties));
  const double rn = std::sqrt(static_cast<double>(n));
  if (x_unique < rn && y_unique < rn) {
    double disc = copula_degree0_unsigned(xv, yv);
    if (!std::isfinite(disc)) disc = dep;
    dep = 0.5 * (dep + disc);
  }
  return dep;
}
