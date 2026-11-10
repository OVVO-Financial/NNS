// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include "central_tendencies.h"

using namespace Rcpp;

// Recursive mean-split partition.
//
// A cell is a contiguous range of an index array, and a split rearranges that
// range into its children in digit order. The cells are therefore always in
// label order ("q1" < "q12" < "q2"), and the labels themselves are only built
// once, for the final cells. Observations are never relabelled or regrouped
// level by level.
//
// The default reducer (gravity) only depends on a cell's sorted values, so the
// values are sorted once and each split keeps its children sorted: a level
// costs one pass over the data. The other reducers gather a cell's values in
// row order, as they depend on it (a long double mean) or sort for themselves.

namespace {

inline bool not_nan(double x) { return !ISNAN(x); }

enum Noise { NOISE_MEAN, NOISE_MEDIAN, NOISE_MODE, NOISE_MODE_CLASS, NOISE_OFF };

inline Noise noise_code(const std::string& s) {
  if (s == "mean") return NOISE_MEAN;
  if (s == "median") return NOISE_MEDIAN;
  if (s == "mode") return NOISE_MODE;
  if (s == "mode_class") return NOISE_MODE_CLASS;
  return NOISE_OFF;
}

// ---- reducers of one cell's values, given in row order ----

double mean_no_na(const std::vector<double>& v) {
  long double s = 0.0L; std::size_t m = 0;
  for (double xi : v) if (not_nan(xi)) { s += xi; ++m; }
  return m ? static_cast<double>(s / m) : NA_REAL;
}

double median_no_na(const std::vector<double>& v, std::vector<double>& a) {
  a.clear();
  for (double xi : v) if (not_nan(xi)) a.push_back(xi);
  if (a.empty()) return NA_REAL;
  std::size_t n = a.size();
  std::nth_element(a.begin(), a.begin() + n / 2, a.end());
  double hi = a[n / 2];
  if (n & 1u) return hi;
  auto lm = std::max_element(a.begin(), a.begin() + n / 2);
  return (*lm + hi) * 0.5;
}

double gravity_finite(const std::vector<double>& v, std::vector<double>& a) {
  a.clear();
  for (double xi : v) if (R_finite(xi)) a.push_back(xi);
  return gravity_inplace(a, false);
}

double mode_of(const std::vector<double>& v, bool discrete) {
  NumericVector r(v.begin(), v.end());
  return as<double>(NNS_mode_cpp(r, discrete, false));
}

struct Agg {
  Noise noise;
  std::vector<double> scratch;

  double for_x(const std::vector<double>& v) {
    switch (noise) {
      case NOISE_MEAN:   return mean_no_na(v);
      case NOISE_MEDIAN: return median_no_na(v, scratch);
      // Regression cells take the continuous mode: the discrete one rounds to
      // integers, which collapses data on a finer scale.
      case NOISE_MODE:   return mode_of(v, false);
      default:           return gravity_finite(v, scratch);
    }
  }

  double for_y(const std::vector<double>& v) {
    switch (noise) {
      case NOISE_MEAN:       return mean_no_na(v);
      case NOISE_MEDIAN:     return median_no_na(v, scratch);
      case NOISE_MODE:       return mode_of(v, false);
      case NOISE_MODE_CLASS: return mode_of(v, true);
      default:               return gravity_finite(v, scratch);
    }
  }
};

struct Node {
  int lo, hi;      // range of the index array(s)
  int parent;      // -1 for the root
  char digit;      // appended to the parent's label
};

// A data.frame assembled here instead of by DataFrame::create(), which calls
// as.data.frame() back in R and took about 95% of a small partition's time.
// Same columns, names, class and automatic (compact) row names.
List plain_data_frame(List cols) {
  const R_xlen_t nr = cols.size() ? Rf_xlength(cols[0]) : 0;
  cols.attr("row.names") = nr > 0 ? IntegerVector::create(NA_INTEGER, -(int)nr)
                                  : IntegerVector(0);
  cols.attr("class") = "data.frame";
  return cols;
}

// Indices ordered by value, ties by index.
void sorted_index(const NumericVector& v, std::vector<int>& idx) {
  const int n = v.size();
  std::vector<std::pair<double, int> > key(n);
  for (int i = 0; i < n; ++i) key[i] = std::make_pair(v[i], i);
  std::sort(key.begin(), key.end());
  idx.resize(n);
  for (int i = 0; i < n; ++i) idx[i] = key[i].second;
}

// Rearrange order[lo, hi) by the digit of each observation, keeping the order
// within a digit.
inline void stable_by_digit(std::vector<int>& order, int lo, int hi,
                            const std::vector<unsigned char>& digit,
                            const int* start, std::vector<int>& tmp) {
  int at[5];
  for (int d = 1; d <= 4; ++d) at[d] = start[d];
  for (int k = lo; k < hi; ++k) {
    const int i = order[k];
    tmp[at[digit[i]]++] = i;
  }
  std::copy(tmp.begin(), tmp.begin() + (hi - lo), order.begin() + lo);
}

}  // namespace

// [[Rcpp::export]]
List NNS_part_cpp(NumericVector x,
                  NumericVector y,
                  Nullable<std::string> type,
                  Nullable<int> order_in,
                  int obs_req,
                  bool min_obs_stop,
                  std::string noise_reduction,
                  bool quadrants_only = false) {

  const int n = x.size();
  if (y.size() != n) stop("x and y must have same length");

  int default_order = std::max((int)std::ceil(std::log2(std::max(1, n))), 1);
  int max_order = order_in.isNotNull() ? as<int>(order_in) : default_order;
  if (max_order == 0) max_order = 1;
  const int level_cap = (int)std::floor(std::log2(std::max(1, n)));
  const bool xonly = type.isNotNull();
  std::transform(noise_reduction.begin(), noise_reduction.end(),
                 noise_reduction.begin(), ::tolower);
  Agg agg;
  agg.noise = noise_code(noise_reduction);

  bool x_finite = true, y_finite = true;
  for (int i = 0; i < n; ++i) {
    if (!R_finite(x[i])) x_finite = false;
    if (!R_finite(y[i])) y_finite = false;
  }
  // Splits on gravity alone can work from sorted values.
  const bool sorted_x = xonly && x_finite &&
    (agg.noise == NOISE_OFF || agg.noise == NOISE_MODE_CLASS);
  const bool sorted_xy = !xonly && x_finite && y_finite && agg.noise == NOISE_OFF;

  // `order` holds every cell as a contiguous range. Sorted splits keep each
  // cell ordered by x (and `order_y` by y); the others keep row order.
  std::vector<int> order, order_y;
  std::vector<double> xs;                 // x in sorted order (x-only splits)
  if (sorted_x) {
    sorted_index(x, order);
    xs.resize(n);
    for (int k = 0; k < n; ++k) xs[k] = x[order[k]];
  } else if (sorted_xy) {
    sorted_index(x, order);
    sorted_index(y, order_y);
  } else {
    order.resize(n);
    for (int i = 0; i < n; ++i) order[i] = i;
  }

  std::vector<Node> nodes;
  nodes.push_back(Node{0, n, -1, 0});
  std::vector<int> cells(1, 0), next;      // current cells, in label order

  std::vector<unsigned char> digit;
  std::vector<int> tmp;
  if (!sorted_x) { digit.assign(n, 0); tmp.resize(n); }
  std::vector<double> xv, yv;

  std::vector<double> H_x0, H_x1, H_y;
  std::vector<double> V_x,  V_y0, V_y1;
  std::vector<double> V_lines;
  const bool want_lines = !quadrants_only;

  int depth = 0;
  while (true) {
    if (depth >= max_order) break;
    if (depth >= level_cap) break;

    bool any_split = false;
    for (int id : cells) {
      if (nodes[id].hi - nodes[id].lo > obs_req) { any_split = true; break; }
    }
    if (!any_split) break;

    if (xonly && want_lines) {
      for (int id : cells) {
        const int lo = nodes[id].lo, hi = nodes[id].hi;
        double minx = R_PosInf, maxx = R_NegInf;
        if (sorted_x) {
          if (hi > lo) { minx = xs[lo]; maxx = xs[hi - 1]; }
        } else {
          for (int k = lo; k < hi; ++k) {
            const double xi = x[order[k]];
            if (not_nan(xi)) { if (xi < minx) minx = xi; if (xi > maxx) maxx = xi; }
          }
        }
        if (not_nan(minx)) V_lines.push_back(minx);
        if (not_nan(maxx)) V_lines.push_back(maxx);
      }
    }

    next.clear();
    for (int id : cells) {
      const int lo = nodes[id].lo, hi = nodes[id].hi, m = hi - lo;
      if (m <= obs_req) { next.push_back(id); continue; }

      int count[5] = {0, 0, 0, 0, 0};

      if (sorted_x) {
        // The cell is a run of sorted x: the values at or below its gravity
        // are its first `left` values.
        const double* b = xs.data() + lo;
        const double cx = gravity_sorted(b, m, false);
        const int left = not_nan(cx) ? (int)(std::upper_bound(b, b + m, cx) - b) : m;
        count[1] = left;
        count[2] = m - left;
      } else if (sorted_xy) {
        xv.resize(m); yv.resize(m);
        for (int k = 0; k < m; ++k) { xv[k] = x[order[lo + k]]; yv[k] = y[order_y[lo + k]]; }
        const double cx = gravity_sorted(xv.data(), m, false);
        const double cy = gravity_sorted(yv.data(), m, false);
        if (want_lines) {
          if (not_nan(cy)) { H_x0.push_back(xv[0]); H_x1.push_back(xv[m - 1]); H_y.push_back(cy); }
          if (not_nan(cx)) { V_x.push_back(cx); V_y0.push_back(yv[0]); V_y1.push_back(yv[m - 1]); }
        }
        const bool okx = not_nan(cx), oky = not_nan(cy);
        for (int k = lo; k < hi; ++k) {
          const int i = order[k];
          const int lox = okx ? (x[i] <= cx) : 0;
          const int loy = oky ? (y[i] <= cy) : 0;
          const int q = 1 + lox + 2 * loy;
          digit[i] = (unsigned char)q;
          ++count[q];
        }
      } else {
        xv.resize(m);
        for (int k = 0; k < m; ++k) xv[k] = x[order[lo + k]];
        const double cx = agg.for_x(xv);
        double cy = NA_REAL;
        if (!xonly) {
          yv.resize(m);
          for (int k = 0; k < m; ++k) yv[k] = y[order[lo + k]];
          cy = agg.for_y(yv);
          if (want_lines) {
            double minx = R_PosInf, maxx = R_NegInf, miny = R_PosInf, maxy = R_NegInf;
            for (int k = 0; k < m; ++k) {
              const double xi = xv[k], yi = yv[k];
              if (not_nan(xi)) { if (xi < minx) minx = xi; if (xi > maxx) maxx = xi; }
              if (not_nan(yi)) { if (yi < miny) miny = yi; if (yi > maxy) maxy = yi; }
            }
            if (not_nan(cy) && not_nan(minx) && not_nan(maxx)) {
              H_x0.push_back(minx); H_x1.push_back(maxx); H_y.push_back(cy);
            }
            if (not_nan(cx) && not_nan(miny) && not_nan(maxy)) {
              V_x.push_back(cx); V_y0.push_back(miny); V_y1.push_back(maxy);
            }
          }
        }
        for (int k = lo; k < hi; ++k) {
          const int i = order[k];
          int q;
          if (!xonly) {
            const int lox = (not_nan(x[i]) && not_nan(cx)) ? (x[i] <= cx) : 0;
            const int loy = (not_nan(y[i]) && not_nan(cy)) ? (y[i] <= cy) : 0;
            q = 1 + lox + 2 * loy;
          } else {
            const int lox = (not_nan(x[i]) && not_nan(cx)) ? (x[i] > cx) : 0;
            q = 1 + lox;
          }
          digit[i] = (unsigned char)q;
          ++count[q];
        }
      }

      int start[5];
      start[1] = 0;
      for (int d = 2; d <= 4; ++d) start[d] = start[d - 1] + count[d - 1];
      if (!sorted_x) {
        stable_by_digit(order, lo, hi, digit, start, tmp);
        if (sorted_xy) stable_by_digit(order_y, lo, hi, digit, start, tmp);
      }
      for (int d = 1; d <= 4; ++d) {
        if (!count[d]) continue;
        next.push_back((int)nodes.size());
        nodes.push_back(Node{lo + start[d], lo + start[d] + count[d], id, (char)('0' + d)});
      }
    }
    cells.swap(next);
    ++depth;

    if (min_obs_stop) {
      int minc = n;
      for (int id : cells) minc = std::min(minc, nodes[id].hi - nodes[id].lo);
      if (minc <= obs_req) break;
    }
  }

  // ---- labels: one string per final cell ----
  std::string buf;
  auto label = [&](int id) {
    buf.clear();
    for (int c = id; c > 0; c = nodes[c].parent) buf.push_back(nodes[c].digit);
    buf.push_back('q');
    std::reverse(buf.begin(), buf.end());
  };

  CharacterVector q_cur(n);
  for (int id : cells) {
    if (nodes[id].hi == nodes[id].lo) continue;
    label(id);
    Shield<SEXP> s(Rf_mkCharLen(buf.data(), (int)buf.size()));
    for (int k = nodes[id].lo; k < nodes[id].hi; ++k) SET_STRING_ELT(q_cur, order[k], s);
  }
  if (quadrants_only) return List::create(_["quadrant"] = q_cur);

  // A cell's prior label is its parent's; sibling final cells share one
  // regression point. "pq" marks data that was never split.
  CharacterVector q_prior(n);
  std::vector<int> group_of_node(nodes.size(), -1);
  std::vector<int> group(n, 0);
  std::vector<std::string> rp_q;
  int root_group = -1;
  for (int id : cells) {
    if (nodes[id].hi == nodes[id].lo) continue;
    const int parent = nodes[id].parent;
    int g;
    if (parent < 0) {
      if (root_group < 0) { root_group = (int)rp_q.size(); rp_q.push_back("pq"); }
      g = root_group;
    } else {
      if (group_of_node[parent] < 0) {
        group_of_node[parent] = (int)rp_q.size();
        label(parent);
        rp_q.push_back(buf);
      }
      g = group_of_node[parent];
    }
    Shield<SEXP> s(Rf_mkCharLen(rp_q[g].data(), (int)rp_q[g].size()));
    for (int k = nodes[id].lo; k < nodes[id].hi; ++k) {
      SET_STRING_ELT(q_prior, order[k], s);
      group[order[k]] = g;
    }
  }

  List part = plain_data_frame(List::create(_["x"] = x, _["y"] = y, _["quadrant"] = q_cur,
                                            _["prior.quadrant"] = q_prior));

  // Regression points: each group's values in row order, as the reducers
  // that depend on it require.
  const int G = (int)rp_q.size();
  std::vector<int> g_start(G + 1, 0);
  for (int i = 0; i < n; ++i) ++g_start[group[i] + 1];
  for (int g = 0; g < G; ++g) g_start[g + 1] += g_start[g];
  std::vector<int> g_at(g_start.begin(), g_start.end() - 1);
  std::vector<double> gx(n), gy(n);
  for (int i = 0; i < n; ++i) {
    const int at = g_at[group[i]]++;
    gx[at] = x[i];
    gy[at] = y[i];
  }
  std::vector<double> rp_x(G), rp_y(G);
  for (int g = 0; g < G; ++g) {
    xv.assign(gx.begin() + g_start[g], gx.begin() + g_start[g + 1]);
    yv.assign(gy.begin() + g_start[g], gy.begin() + g_start[g + 1]);
    rp_x[g] = agg.for_x(xv);
    rp_y[g] = agg.for_y(yv);
  }

  List rp = plain_data_frame(List::create(_["quadrant"] = wrap(rp_q),
                                          _["x"] = wrap(rp_x),
                                          _["y"] = wrap(rp_y)));

  List seg_h = plain_data_frame(List::create(_["x0"] = wrap(H_x0),
                                             _["x1"] = wrap(H_x1),
                                             _["y"] = wrap(H_y)));
  List seg_v = plain_data_frame(List::create(_["x"] = wrap(V_x),
                                             _["y0"] = wrap(V_y0),
                                             _["y1"] = wrap(V_y1)));

  return List::create(_["order"] = depth,
                      _["dt"] = part,
                      _["regression.points"] = rp,
                      _["segments_h"] = seg_h,
                      _["segments_v"] = seg_v,
                      _["vlines"] = wrap(V_lines));
}
