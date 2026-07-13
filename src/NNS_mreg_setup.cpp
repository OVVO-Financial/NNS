#include <Rcpp.h>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <unordered_map>
using namespace Rcpp;

double gravity_value(std::vector<double> x, bool discrete);

namespace {
double reduce_numeric(std::vector<double>& v, int reducer, bool is_class) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x){ return !R_finite(x); }), v.end());
  if (v.empty()) return NA_REAL;
  std::sort(v.begin(), v.end());
  if (is_class || reducer == 2) {
    double best = v[0]; int bestn = 0;
    for (std::size_t i = 0; i < v.size();) {
      const double val = v[i]; int n = 0;
      while (i < v.size() && v[i] == val) { ++i; ++n; }
      if (n > bestn) { bestn = n; best = val; }
    }
    return best;
  }
  if (reducer == 1) {
    const std::size_t n = v.size();
    return (n % 2) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
  }
  if (reducer == 3) return gravity_value(v, false);
  double s = 0.0; for (double x : v) s += x;
  return s / static_cast<double>(v.size());
}

int find_interval_one(double x, const NumericVector& b) {
  const int n = b.size();
  if (n == 0) return 0;
  // Match R's one-based findInterval(..., left.open = FALSE,
  // rightmost.closed = TRUE) IDs exactly:
  // below first -> 0; exact final boundary -> n - 1; above final -> n.
  // The string IDs are sorted lexicographically by split(), so parity here is
  // required to preserve RPM row ordering and stable tie-break behavior.
  if (x < b[0]) return 0;
  if (x == b[n - 1]) return std::max(0, n - 1);
  return static_cast<int>(std::upper_bound(b.begin(), b.end(), x) - b.begin());
}

std::string id_for_row(const NumericMatrix& X, const List& boundaries, int r) {
  std::ostringstream os;
  for (int j = 0; j < X.ncol(); ++j) {
    if (j) os << ".";
    os << find_interval_one(X(r, j), boundaries[j]);
  }
  return os.str();
}
}

// [[Rcpp::export]]
IntegerVector NNS_duplicate_column_map_cpp(const NumericMatrix& X) {
  const int n = X.nrow(), p = X.ncol();
  IntegerVector rep(p);
  std::unordered_map<std::string, std::vector<int> > seen;
  for (int j = 0; j < p; ++j) {
    std::string key; key.reserve(n * sizeof(double));
    for (int i = 0; i < n; ++i) {
      const double v = X(i, j);
      key.append(reinterpret_cast<const char*>(&v), sizeof(double));
    }
    int r = j;
    auto it = seen.find(key);
    if (it != seen.end()) {
      for (int cand : it->second) {
        bool same = true;
        for (int i = 0; i < n && same; ++i) same = (X(i, j) == X(i, cand));
        if (same) { r = cand; break; }
      }
      it->second.push_back(j);
    } else {
      seen[key] = std::vector<int>(1, j);
    }
    rep[j] = r + 1;
  }
  return rep;
}

// [[Rcpp::export]]
List NNS_mreg_setup_cpp(const NumericMatrix& X, const NumericVector& y,
                        const List& boundaries, int reducer_code, bool is_class) {
  const int nr = X.nrow(), p = X.ncol();
  if (y.size() != nr) stop("y length must equal nrow(X)");
  CharacterVector ids(nr);
  std::unordered_map<std::string, std::vector<int> > groups;
  groups.reserve(static_cast<std::size_t>(nr) * 2);
  for (int i = 0; i < nr; ++i) {
    const std::string id = id_for_row(X, boundaries, i);
    ids[i] = id;
    groups[id].push_back(i);
  }
  std::vector<std::string> ordered_ids;
  ordered_ids.reserve(groups.size());
  for (const auto& kv : groups) ordered_ids.push_back(kv.first);
  std::sort(ordered_ids.begin(), ordered_ids.end());

  NumericMatrix rpm(ordered_ids.size(), p + 1);
  std::vector<double> vals;
  for (std::size_t g = 0; g < ordered_ids.size(); ++g) {
    const std::vector<int>& idx = groups[ordered_ids[g]];
    for (int j = 0; j < p; ++j) {
      vals.clear(); vals.reserve(idx.size());
      for (int row : idx) vals.push_back(X(row, j));
      rpm(g, j) = reduce_numeric(vals, reducer_code, false);
    }
    vals.clear(); vals.reserve(idx.size());
    for (int row : idx) vals.push_back(y[row]);
    rpm(g, p) = reduce_numeric(vals, reducer_code, is_class);
  }
  CharacterVector row_ids(ordered_ids.size());
  for (std::size_t i = 0; i < ordered_ids.size(); ++i) row_ids[i] = ordered_ids[i];
  return List::create(_["RPM"] = rpm, _["ids"] = ids, _["row_ids"] = row_ids);
}
