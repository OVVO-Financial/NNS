// src/internal_functions.cpp
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>

using namespace Rcpp;

// ---------- small utilities ----------

inline bool isFactor(SEXP x) {
  return TYPEOF(x) == INTSXP && Rf_isFactor(x);
}

inline int presentLevelCount(IntegerVector f) {
  int L = Rf_length(Rf_getAttrib(f, R_LevelsSymbol));
  std::vector<char> seen(L + 1, 0);
  for (int i = 0; i < f.size(); ++i) {
    int k = f[i];
    if (k != NA_INTEGER) seen[k] = 1;
  }
  int cnt = 0;
  for (int k = 1; k <= L; ++k) cnt += seen[k];
  return cnt;
}

NumericVector flattenNumericNoNA(SEXP x) {
  std::vector<double> out;
  if (Rf_isNull(x)) return NumericVector(0);
  if (Rf_isVectorAtomic(x) && TYPEOF(x) != STRSXP) {
    NumericVector v = as<NumericVector>(x);
    out.reserve(v.size());
    for (double z : v) if (R_finite(z)) out.push_back(z);
  } else if (TYPEOF(x) == VECSXP) {
    List L(x);
    for (int i = 0; i < L.size(); ++i) {
      NumericVector v = flattenNumericNoNA(L[i]);
      for (double z : v) if (R_finite(z)) out.push_back(z);
    }
  }
  return wrap(out);
}

inline void setColNames(NumericMatrix &m, const CharacterVector &nm) {
  colnames(m) = nm;
}

IntegerVector sampleWithoutReplacement(const IntegerVector& idx, int k) {
  IntegerVector pool = clone(idx);
  if (k >= pool.size()) return pool;
  for (int i = 0; i < k; ++i) {
    int j = i + (int)std::floor(R::runif(0.0, 1.0) * (pool.size() - i));
    std::swap(pool[i], pool[j]);
  }
  IntegerVector out(no_init(k));
  std::copy(pool.begin(), pool.begin() + k, out.begin());
  return out;
}

NumericVector colSd(const NumericMatrix &M) {
  int n = M.nrow(), p = M.ncol();
  NumericVector sds(no_init(p));
  if (n <= 1) {
    std::fill(sds.begin(), sds.end(), NA_REAL);
    return sds;
  }
  for (int j = 0; j < p; ++j) {
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += M(i, j);
    mu /= (double)n;
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      double d = M(i, j) - mu;
      ss += d * d;
    }
    sds[j] = std::sqrt(ss / (double)(n - 1));
  }
  return sds;
}

double vecSd(const NumericVector &x) {
  int n = x.size();
  if (n <= 1) return NA_REAL;
  double mu = 0.0;
  for (int i = 0; i < n; ++i) mu += x[i];
  mu /= (double)n;
  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = x[i] - mu;
    ss += d * d;
  }
  return std::sqrt(ss / (double)(n - 1));
}

// ---------- 1) is.fcl ----------

// [[Rcpp::export(name = "is.fcl")]]
bool is_fcl(SEXP x) {
  return isFactor(x) || TYPEOF(x) == STRSXP || TYPEOF(x) == LGLSXP;
}

// ---------- 2) is.discrete ----------
// [[Rcpp::export(name = "is.discrete")]]
bool is_discrete(SEXP x) {
  if (TYPEOF(x) == INTSXP || TYPEOF(x) == LGLSXP) return true;
  if (TYPEOF(x) == REALSXP) {
    NumericVector v(x);
    for (int i = 0; i < v.size(); ++i) {
      if (NumericVector::is_na(v[i])) continue;
      if (v[i] != std::trunc(v[i])) return false;
    }
    return true;
  }
  return false;
}

// ---------- 3) factor_2_dummy & factor_2_dummy_FR ----------
// [[Rcpp::export]]
SEXP factor_2_dummy(SEXP x) {
  while (TYPEOF(x) == VECSXP && !Rf_isFactor(x)) {
    List L(x);
    if (L.size() == 1) x = L[0];
    else break;
  }
  if (isFactor(x)) {
    IntegerVector f(x);
    int L = Rf_length(Rf_getAttrib(x, R_LevelsSymbol));
    int present = presentLevelCount(f);
    if (present <= 1) return as<NumericVector>(x);
    
    int n = f.size();
    int cols = std::max(0, L - 1);
    NumericMatrix out(n, cols);
    for (int i = 0; i < n; ++i) {
      int k = f[i];
      if (k != NA_INTEGER && k > 1) out(i, k - 2) = 1.0;
    }
    CharacterVector lev = Rf_getAttrib(x, R_LevelsSymbol);
    if (cols > 0) {
      CharacterVector cn(no_init(cols));
      for (int j = 0; j < cols; ++j) cn[j] = lev[j + 1];
      setColNames(out, cn);
    }
    return out;
  }
  return as<NumericVector>(x);
}

// [[Rcpp::export]]
SEXP factor_2_dummy_FR(SEXP x) {
  while (TYPEOF(x) == VECSXP && !Rf_isFactor(x)) {
    List L(x);
    if (L.size() == 1) x = L[0];
    else break;
  }
  if (isFactor(x)) {
    IntegerVector f(x);
    int L = Rf_length(Rf_getAttrib(x, R_LevelsSymbol));
    int present = presentLevelCount(f);
    if (present <= 1) return as<NumericVector>(x);
    
    int n = f.size();
    NumericMatrix out(n, L);
    for (int i = 0; i < n; ++i) {
      int k = f[i];
      if (k != NA_INTEGER) out(i, k - 1) = 1.0;
    }
    CharacterVector lev = Rf_getAttrib(x, R_LevelsSymbol);
    setColNames(out, lev);
    return out;
  }
  return as<NumericVector>(x);
}

// ---------- 4) generate.vectors ----------
// [[Rcpp::export(name = "generate.vectors")]]
List generate_vectors(NumericVector x, IntegerVector l) {
  int n = x.size();
  List comp_series(l.size()), comp_index(l.size());
  for (int t = 0; t < l.size(); ++t) {
    int lag = l[t];
    if (lag <= 0) {
      comp_series[t] = NumericVector(0);
      comp_index[t] = IntegerVector(0);
      continue;
    }
    int start = (n % lag) + 1;
    int m = ((n - start) / lag) + 1;
    
    NumericVector s(no_init(m));
    IntegerVector idx(no_init(m));
    
    int pos = start;
    for (int i = 0; i < m; ++i, pos += lag) {
      s[i] = x[pos - 1];
      idx[i] = i + 1;
    }
    comp_series[t] = s;
    comp_index[t] = idx;
  }
  return List::create(_["Component.index"] = comp_index, _["Component.series"] = comp_series);
}

// ---------- 5) generate.lin.vectors ----------
static List create_recycled_list_int(const IntegerVector& values, int list_length) {
  List res(list_length);
  std::vector< std::vector<int> > buckets(list_length);
  for (int i = 0; i < values.size(); ++i) buckets[(i % list_length)].push_back(values[i]);
  for (int j = 0; j < list_length; ++j) res[j] = buckets[j].empty() ? R_NilValue : wrap(buckets[j]);
  return res;
}

static List create_recycled_list_num(const NumericVector& values, int list_length) {
  List res(list_length);
  std::vector< std::vector<double> > buckets(list_length);
  for (int i = 0; i < values.size(); ++i) buckets[(i % list_length)].push_back(values[i]);
  for (int j = 0; j < list_length; ++j) res[j] = buckets[j].empty() ? R_NilValue : wrap(buckets[j]);
  return res;
}

// [[Rcpp::export(name = "generate.lin.vectors")]]
List generate_lin_vectors(NumericVector x, int l, int h = 1) {
  int n = x.size();
  int max_fcast = std::min(h, l);
  
  List comp_series(max_fcast), comp_index(max_fcast);
  for (int i = 1; i <= max_fcast; ++i) {
    int start = ((n + i - 1) % l) + 1;
    int m = ((n - start) / l) + 1;
    NumericVector s(no_init(m));
    IntegerVector idx(no_init(m));
    int pos = start;
    for (int k = 0; k < m; ++k, pos += l) {
      s[k] = x[pos - 1];
      idx[k] = k + 1;
    }
    comp_series[i - 1] = s;
    comp_index[i - 1] = idx;
  }
  
  IntegerVector one_to_h(no_init(h));
  for (int i = 0; i < h; ++i) one_to_h[i] = i + 1;
  List forecast_index = create_recycled_list_int(one_to_h, max_fcast);
  
  NumericVector fvals(no_init(h));
  for (int i = 1; i <= h; ++i) {
    int ci = ((((i - 1) % l)) % std::max(1, max_fcast)) + 1;
    int last_val = as<IntegerVector>(comp_index[ci - 1]).size();
    fvals[i - 1] = (double)last_val + std::ceil((double)i / (double)l);
  }
  
  return List::create(_["Component.index"] = comp_index,
                      _["Component.series"] = comp_series,
                      _["forecast.values"] = create_recycled_list_num(fvals, l),
                      _["forecast.index"] = forecast_index);
}

// ---------- 6) ARMA.seas.weighting ----------
// [[Rcpp::export(name = "ARMA.seas.weighting")]]
List ARMA_seas_weighting(bool sf, SEXP mat) {
  if (!Rf_isMatrix(mat) && !(Rf_inherits(mat, "data.frame")) && TYPEOF(mat) != VECSXP) {
    NumericVector M = as<NumericVector>(mat);
    double lag = (M.size() > 0 && !NumericVector::is_na(M[0])) ? M[0] : NA_REAL;
    return List::create(_["lag"] = lag, _["Weights"] = 1.0);
  }
  
  int n = NA_INTEGER;
  if (Rf_isMatrix(mat)) {
    IntegerVector dims = Rf_getAttrib(mat, R_DimSymbol);
    if (dims.size() == 2) n = dims[1];
  } else if (Rf_inherits(mat, "data.frame") || TYPEOF(mat) == VECSXP) {
    n = List(mat).size();
  }
  
  if (n == 1) return List::create(_["lag"] = 1, _["Weights"] = 1.0);
  List M(mat);
  
  if (sf) {
    if (M.containsElementNamed("all.periods")) {
      SEXP ap = M["all.periods"];
      if (Rf_inherits(ap, "data.frame") || TYPEOF(ap) == VECSXP) {
        List AP(ap);
        if (AP.containsElementNamed("Period")) {
          NumericVector Period = as<NumericVector>(AP["Period"]);
          double lag_scalar = (Period.size() > 0 && !NumericVector::is_na(Period[0])) ? Period[0] : NA_REAL;
          return List::create(_["lag"] = lag_scalar, _["Weights"] = 1.0);
        }
      }
    }
    return List::create(_["lag"] = 1.0, _["Weights"] = 1.0);
  }
  
  NumericVector lag = flattenNumericNoNA(M.containsElementNamed("Period") ? M["Period"] : R_NilValue);
  NumericVector observation_weighting(no_init(lag.size()));
  for (int i = 0; i < lag.size(); ++i) observation_weighting[i] = 1.0 / std::sqrt(lag[i]);
  
  NumericVector covar = flattenNumericNoNA(M.containsElementNamed("Coefficient.of.Variation") ? M["Coefficient.of.Variation"] : R_NilValue);
  NumericVector varcovar = flattenNumericNoNA(M.containsElementNamed("Variable.Coefficient.of.Variation") ? M["Variable.Coefficient.of.Variation"] : R_NilValue);
  
  NumericVector lag_weighting;
  if (covar.size() == 1 && NumericVector::is_na(covar[0])) {
    lag_weighting = NumericVector(varcovar.size(), 1.0);
  } else {
    int m = std::min(varcovar.size(), covar.size());
    lag_weighting = NumericVector(no_init(m));
    for (int i = 0; i < m; ++i) lag_weighting[i] = varcovar[i] - covar[i];
    observation_weighting = observation_weighting[Rcpp::Range(0, m - 1)];
  }
  
  NumericVector wprod(no_init(lag_weighting.size()));
  double denom = 0.0;
  for (int i = 0; i < wprod.size() && i < observation_weighting.size(); ++i) {
    wprod[i] = lag_weighting[i] * observation_weighting[i];
    denom += wprod[i];
  }
  return List::create(_["lag"] = lag, _["Weights"] = (denom == 0.0) ? NumericVector(wprod.size(), 0.0) : wprod / denom);
}

// ---------- 8) NNS.meboot.part ----------
// [[Rcpp::export(name = "NNS.meboot.part")]]
NumericVector NNS_meboot_part(NumericVector xx, int n, NumericVector z,
                              double xmin, double xmax, NumericVector desintxb, bool reachbnd) {
  NumericVector p = runif(n);
  int m = xx.size();
  NumericVector q(no_init(n));
  
  if (m == 0) {
    std::fill(q.begin(), q.end(), NA_REAL);
  } else if (m == 1) {
    std::fill(q.begin(), q.end(), xx[0]);
  } else {
    for (int i = 0; i < n; ++i) {
      double pi = p[i];
      if (pi <= 0.0) { q[i] = xx[0]; continue; }
      if (pi >= 1.0) { q[i] = xx[m - 1]; continue; }
      double h = 1.0 + (m - 1.0) * pi;
      int j = static_cast<int>(std::floor(h));
      if (j < 1) j = 1; else if (j > m - 1) j = m - 1;
      q[i] = (1.0 - (h - j)) * xx[j - 1] + (h - j) * xx[j];
    }
  }
  
  double invn = 1.0 / (double)n;
  for (int i = 0; i < p.size(); ++i) {
    if (p[i] <= invn) {
      double val = xmin + (p[i] - 0.0) * (z[0] - xmin) / (invn - 0.0);
      if (!reachbnd) val = val + desintxb[0] - 0.5 * (z[0] + xmin);
      q[i] = val;
    }
  }
  
  double edge = (double)(n - 1) / (double)n;
  for (int i = 0; i < p.size(); ++i) {
    if (p[i] >= edge) {
      double val = z[n - 2] + (p[i] - edge) * (xmax - z[n - 2]) / (1.0 - edge);
      if (!reachbnd) val = val + desintxb[n - 1] - 0.5 * (z[n - 2] + xmax);
      q[i] = val;
    }
  }
  return q;
}

// ---------- 9) NNS.meboot.expand.sd ----------
// [[Rcpp::export(name = "NNS.meboot.expand.sd")]]
SEXP NNS_meboot_expand_sd(SEXP x, NumericMatrix ensemble, double fiv = 5.0) {
  NumericVector sdx;
  if (Rf_isMatrix(x) || Rf_inherits(x, "data.frame")) {
    NumericMatrix X;
    if (Rf_inherits(x, "data.frame")) {
      DataFrame DF = as<DataFrame>(x);
      int nr = (Rf_length(x) == 0) ? 0 : Rf_length(VECTOR_ELT(x, 0));
      int p  = Rf_length(x);
      X = NumericMatrix(nr, p);
      for (int j = 0; j < p; ++j) X(_, j) = as<NumericVector>(DF[j]);
    } else X = as<NumericMatrix>(x);
    sdx = colSd(X);
  } else sdx = NumericVector::create(vecSd(as<NumericVector>(x)));
  
  NumericVector ens_sd = colSd(ensemble);
  NumericVector sdf(no_init(sdx.size() + ens_sd.size()));
  int pos = 0;
  for (double v : sdx) sdf[pos++] = v;
  for (double v : ens_sd) sdf[pos++] = v;
  
  NumericVector sdfa(no_init(sdf.size()));
  NumericVector sdfd(no_init(sdf.size()));
  for (int i = 0; i < sdf.size(); ++i) {
    sdfa[i] = sdf[i] / sdf[0];
    sdfd[i] = sdf[0] / sdf[i];
  }
  
  double mx = 1.0 + (fiv / 100.0);
  for (int i = 0; i < sdfa.size(); ++i) if (sdfa[i] < 1.0) sdfa[i] = R::runif(1.0, mx);
  
  int J = ensemble.ncol();
  for (int j = 0; j < J; ++j) {
    double a = sdfd[j + 1] * sdfa[j + 1];
    if (std::floor(a) > 0.0) {
      for (int i = 0; i < ensemble.nrow(); ++i) ensemble(i, j) *= a;
    }
  }
  
  if (Rf_inherits(x, "ts")) {
    ensemble.attr("class") = "ts";
    ensemble.attr("tsp") = Rf_getAttrib(x, Rf_install("tsp"));
  }
  return ensemble;
}

// ---------- 10) force.clt ----------
// [[Rcpp::export(name = "force.clt")]]
SEXP force_clt(SEXP x, NumericMatrix ensemble) {
  int n = ensemble.nrow(), bigj = ensemble.ncol();
  double gm = NA_REAL;
  NumericVector s;
  
  if (Rf_isMatrix(x) || Rf_inherits(x, "data.frame")) {
    NumericMatrix X;
    if (Rf_inherits(x, "data.frame")) {
      DataFrame DF = as<DataFrame>(x);
      int nr = (Rf_length(x) == 0) ? 0 : Rf_length(VECTOR_ELT(x, 0));
      int p  = Rf_length(x);
      X = NumericMatrix(nr, p);
      for (int j = 0; j < p; ++j) X(_, j) = as<NumericVector>(DF[j]);
    } else X = as<NumericMatrix>(x);
    double sumAll = 0.0;
    for (int i = 0; i < X.nrow(); ++i) for (int j = 0; j < X.ncol(); ++j) sumAll += X(i, j);
    gm = sumAll / (double)(X.nrow() * X.ncol());
    s = colSd(X);
  } else {
    NumericVector xv = as<NumericVector>(x);
    double sumAll = 0.0;
    for (int i = 0; i < xv.size(); ++i) sumAll += xv[i];
    gm = sumAll / (double)xv.size();
    s = NumericVector::create(vecSd(xv));
  }
  
  NumericVector xbar(no_init(bigj));
  for (int j = 0; j < bigj; ++j) {
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += ensemble(i, j);
    xbar[j] = mu / (double)n;
  }
  
  IntegerVector oo(no_init(bigj));
  for (int j = 0; j < bigj; ++j) oo[j] = j;
  std::sort(oo.begin(), oo.end(), [&](int a, int b){ return xbar[a] < xbar[b]; });
  
  NumericVector sortxbar = clone(xbar);
  std::sort(sortxbar.begin(), sortxbar.end());
  
  NumericVector smean = clone(s);
  for (int i = 0; i < smean.size(); ++i) smean[i] = s[i] / std::sqrt((double)bigj);
  double smean_scalar = smean.size() ? smean[0] : 0.0;
  
  NumericVector newbar(no_init(bigj));
  for (int j = 0; j < bigj; ++j) {
    double sm = (smean.size() == 1) ? smean_scalar : smean[j % smean.size()];
    newbar[j] = gm + R::qnorm((double)(j + 1) / (double)(bigj + 1), 0.0, 1.0, 1, 0) * sm;
  }
  
  double mu_nb = 0.0, ss_nb = 0.0;
  for (int j = 0; j < bigj; ++j) mu_nb += newbar[j];
  mu_nb /= (double)bigj;
  for (int j = 0; j < bigj; ++j) ss_nb += (newbar[j] - mu_nb) * (newbar[j] - mu_nb);
  double sd_nb = std::sqrt(ss_nb / (double)(bigj - 1));
  
  NumericMatrix out = clone(ensemble);
  for (int i = 0; i < bigj; ++i) {
    int col = oo[i];
    double sm = (smean.size() == 1) ? smean_scalar : smean[i % smean.size()];
    double add = (((newbar[i] - mu_nb) / sd_nb) * sm + gm) - sortxbar[i];
    for (int r = 0; r < n; ++r) out(r, col) = ensemble(r, col) + add;
  }
  
  if (Rf_inherits(x, "ts")) {
    out.attr("class") = "ts";
    out.attr("tsp") = Rf_getAttrib(x, Rf_install("tsp"));
  }
  return out;
}

// ---------- 11) downSample / upSample ----------

static IntegerVector subset_factor_codes(const IntegerVector& codes, const IntegerVector& rows) {
  IntegerVector out(no_init(rows.size()));
  for (int i = 0; i < rows.size(); ++i) {
    const int idx = rows[i] - 1;
    out[i] = (idx >= 0 && idx < codes.size()) ? codes[idx] : NA_INTEGER;
  }
  return out;
}

template <int RTYPE>
static Rcpp::Vector<RTYPE> subset_vec_template(const Rcpp::Vector<RTYPE>& v, const IntegerVector& rows) {
  Rcpp::Vector<RTYPE> out(no_init(rows.size()));
  for (int i = 0; i < rows.size(); ++i) {
    const int idx = rows[i] - 1;
    out[i] = (idx >= 0 && idx < v.size()) ? v[idx] : Rcpp::Vector<RTYPE>::get_na();
  }
  return out;
}

template <>
inline Rcpp::CharacterVector subset_vec_template<STRSXP>(const Rcpp::CharacterVector& v, const IntegerVector& rows) {
  Rcpp::CharacterVector out(no_init(rows.size()));
  for (int i = 0; i < rows.size(); ++i) {
    const int idx = rows[i] - 1;
    if (idx >= 0 && idx < v.size()) {
      out[i] = v[idx];
    } else {
      out[i] = NA_STRING;
    }
  }
  return out;
}

static DataFrame subset_df_rows_with_y(const DataFrame& X, const IntegerVector& rows, SEXP y_factor, const std::string& yname, bool include_y) {
  const int p = X.size(), m = rows.size();
  List out(include_y ? (p + 1) : p);
  CharacterVector out_names(include_y ? (p + 1) : p);
  CharacterVector in_names = X.names();
  
  for (int j = 0; j < p; ++j) {
    SEXP col = X[j];
    out_names[j] = in_names[j];
    switch (TYPEOF(col)) {
    case INTSXP: {
      IntegerVector iv(col);
      RObject cls = iv.attr("class");
      if (!cls.isNULL() && as<CharacterVector>(cls).size() > 0 && as<CharacterVector>(cls)[0] == "factor") {
        IntegerVector sub = subset_factor_codes(iv, rows);
        sub.attr("class")  = iv.attr("class");
        sub.attr("levels") = iv.attr("levels");
        out[j] = sub;
      } else out[j] = subset_vec_template<INTSXP>(iv, rows);
      break;
    }
    case REALSXP: out[j] = subset_vec_template<REALSXP>(NumericVector(col), rows); break;
    case LGLSXP:  out[j] = subset_vec_template<LGLSXP>(LogicalVector(col), rows);  break;
    case STRSXP:  out[j] = subset_vec_template<STRSXP>(CharacterVector(col), rows);break;
    default: {
      CharacterVector cv = as<CharacterVector>(col);
      out[j] = subset_vec_template<STRSXP>(cv, rows);
      break;
    }}
  }
  if (include_y) {
    IntegerVector ycodes = as<IntegerVector>(y_factor);
    IntegerVector ysub   = subset_factor_codes(ycodes, rows);
    ysub.attr("class")   = CharacterVector::create("factor");
    ysub.attr("levels")  = Rf_getAttrib(y_factor, R_LevelsSymbol);
    out[p] = ysub;
    out_names[p] = yname;
  }
  out.attr("names") = out_names;
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -m);
  return DataFrame(out);
}

static IntegerVector sample_indices(int N, int k, bool replace) {
  IntegerVector res(no_init(k));
  if (N <= 0 || k <= 0) return res;
  if (!replace && k > N) k = N;
  
  if (replace) {
    for (int i = 0; i < k; ++i) {
      int draw = 1 + (int)floor(R::runif(0.0, 1.0) * N);
      res[i] = (draw > N) ? N : draw;
    }
  } else {
    std::vector<int> a(N);
    for (int i = 0; i < N; ++i) a[i] = i + 1;
    for (int i = 0; i < k; ++i) {
      int j = i + (int)floor(R::runif(0.0, 1.0) * (N - i));
      if (j >= N) j = N - 1;
      std::swap(a[i], a[j]);
      res[i] = a[i];
    }
  }
  return res;
}

// ---------- downSample ------------------------------------

// [[Rcpp::export]]
SEXP downSample(SEXP x, SEXP y, bool list = false, std::string yname = "Class") {
  RNGScope scope;
  if (!Rf_isFactor(y)) {
    Rcpp::warning("Down-sampling requires a factor variable as the response. The original data was returned.");
    return List::create(_["x"] = as<DataFrame>(x), _["y"] = y);
  }
  
  DataFrame X = as<DataFrame>(x);
  IntegerVector fy = as<IntegerVector>(y);
  CharacterVector lev = Rf_getAttrib(y, R_LevelsSymbol);
  const int n = X.nrows(), L = lev.size();
  if (fy.size() != n) stop("downSample: nrow(x) != length(y)");
  
  std::vector< std::vector<int> > perClass(L);
  for (int i = 0; i < n; ++i) if (fy[i] != NA_INTEGER) perClass[fy[i] - 1].push_back(i + 1);
  
  int minClass = n;
  bool any_ok = false;
  for (int k = 0; k < L; ++k) {
    int sz = (int)perClass[k].size();
    if (sz > 0) { any_ok = true; if (sz < minClass) minClass = sz; }
  }
  if (!any_ok || minClass <= 0) stop("downSample: no non-empty classes.");
  
  std::vector<int> rows_out;
  rows_out.reserve(minClass * L);
  for (int k = 0; k < L; ++k) {
    if (perClass[k].empty()) continue;
    IntegerVector s = sample_indices(perClass[k].size(), minClass, false);
    for (int j = 0; j < s.size(); ++j) rows_out.push_back(perClass[k][ s[j] - 1 ]);
  }
  IntegerVector rows = wrap(rows_out);
  
  if (list) {
    DataFrame Xout = subset_df_rows_with_y(X, rows, R_NilValue, yname, false);
    IntegerVector ysub = subset_factor_codes(fy, rows);
    ysub.attr("class")  = CharacterVector::create("factor");
    ysub.attr("levels") = lev;
    return List::create(_["x"] = Xout, _["y"] = ysub);
  }
  return subset_df_rows_with_y(X, rows, y, yname, true);
}

// ---------- upSample --------------------------------------

// [[Rcpp::export]]
SEXP upSample(SEXP x, SEXP y, bool list = false, std::string yname = "Class") {
  RNGScope scope;
  if (!Rf_isFactor(y)) {
    Rcpp::warning("Up-sampling requires a factor variable as the response. The original data was returned.");
    return List::create(_["x"] = as<DataFrame>(x), _["y"] = y);
  }
  
  DataFrame X = as<DataFrame>(x);
  IntegerVector fy = as<IntegerVector>(y);
  CharacterVector lev = Rf_getAttrib(y, R_LevelsSymbol);
  const int n = X.nrows(), L = lev.size();
  if (fy.size() != n) stop("upSample: nrow(x) != length(y)");
  
  std::vector< std::vector<int> > perClass(L);
  for (int i = 0; i < n; ++i) if (fy[i] != NA_INTEGER) perClass[fy[i] - 1].push_back(i + 1);
  
  int maxClass = 0;
  bool any_ok = false;
  for (int k = 0; k < L; ++k) {
    int sz = (int)perClass[k].size();
    if (sz > 0) { any_ok = true; if (sz > maxClass) maxClass = sz; }
  }
  if (!any_ok || maxClass <= 0) stop("upSample: no non-empty classes.");
  
  std::vector<int> rows_out;
  rows_out.reserve(maxClass * L);
  for (int k = 0; k < L; ++k) {
    if (perClass[k].empty()) continue;
    IntegerVector s = sample_indices(perClass[k].size(), maxClass, true);
    for (int j = 0; j < s.size(); ++j) rows_out.push_back(perClass[k][ s[j] - 1 ]);
  }
  IntegerVector rows = wrap(rows_out);
  
  if (list) {
    DataFrame Xout = subset_df_rows_with_y(X, rows, R_NilValue, yname, false);
    IntegerVector ysub = subset_factor_codes(fy, rows);
    ysub.attr("class")  = CharacterVector::create("factor");
    ysub.attr("levels") = lev;
    return List::create(_["x"] = Xout, _["y"] = ysub);
  }
  return subset_df_rows_with_y(X, rows, y, yname, true);
}
