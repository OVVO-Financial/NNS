// File: src/internal_functions.cpp
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ---------- small utilities ----------

// test if SEXP is a factor
inline bool isFactor(SEXP x) {
  return TYPEOF(x) == INTSXP && Rf_isFactor(x);
}

// get present (observed) level count for a factor
inline int presentLevelCount(IntegerVector f) {
  int L = Rf_length(Rf_getAttrib(f, R_LevelsSymbol));
  std::vector<char> seen(L + 1, 0); // codes are 1..L, NA is NA_INTEGER
  for (int i = 0; i < f.size(); ++i) {
    int k = f[i];
    if (k != NA_INTEGER) seen[k] = 1;
  }
  int cnt = 0;
  for (int k = 1; k <= L; ++k) cnt += seen[k];
  return cnt;
}

// flatten numeric content from possibly nested list / vector (na.omit)
NumericVector flattenNumericNoNA(SEXP x) {
  std::vector<double> out;
  if (Rf_isNull(x)) return NumericVector(0);
  if (Rf_isVectorAtomic(x) && TYPEOF(x) != STRSXP) {
    NumericVector v = as<NumericVector>(x);
    for (double z : v) if (!Rcpp::NumericVector::is_na(z)) out.push_back(z);
  } else if (TYPEOF(x) == VECSXP) {
    List L(x);
    for (int i = 0; i < L.size(); ++i) {
      NumericVector v = flattenNumericNoNA(L[i]);
      for (double z : v) if (!Rcpp::NumericVector::is_na(z)) out.push_back(z);
    }
  }
  return wrap(out);
}

// set colnames helper
inline void setColNames(NumericMatrix &m, const CharacterVector &nm) {
  Rf_setAttrib(m, R_DimNamesSymbol, List::create(R_NilValue, nm));
}

// sample k distinct indices from IntegerVector idx (without replacement)
IntegerVector sampleWithoutReplacement(const IntegerVector& idx, int k) {
  IntegerVector pool = clone(idx);
  if (k >= pool.size()) return pool;
  // Fisher-Yates partial shuffle using R RNG
  for (int i = 0; i < k; ++i) {
    int j = i + (int)std::floor(R::runif(0.0, 1.0) * (pool.size() - i));
    std::swap(pool[i], pool[j]);
  }
  IntegerVector out(k);
  std::copy(pool.begin(), pool.begin() + k, out.begin());
  return out;
}

// compute column-wise sd for a numeric matrix (unbiased, denom n-1)
NumericVector colSd(const NumericMatrix &M) {
  int n = M.nrow(), p = M.ncol();
  NumericVector sds(p);
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

// compute sd for a numeric vector (unbiased, denom n-1)
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

// Replicates: sum(as.numeric(x) %% 1) == 0
// [[Rcpp::export(name = "is.discrete")]]
bool is_discrete(SEXP x) {
  NumericVector v = as<NumericVector>(x); // as.numeric
  long double acc = 0.0L;
  for (int i = 0; i < v.size(); ++i) {
    if (NumericVector::is_na(v[i])) continue; // mimic na.rm=TRUE effect
    long double rem = (long double)v[i] - std::floor((long double)v[i]); // x %% 1
    acc += rem;
    if (acc > 0) return false; // early exit
  }
  return acc == 0.0L;
}

// ---------- 3) factor_2_dummy & factor_2_dummy_FR ----------

// [[Rcpp::export]]
SEXP factor_2_dummy(SEXP x) {
  // unlist(x)
  while (TYPEOF(x) == VECSXP && !Rf_isFactor(x)) {
    List L(x);
    if (L.size() == 1) x = L[0];
    else break;
  }
  if (isFactor(x)) {
    IntegerVector f(x);
    int L = Rf_length(Rf_getAttrib(x, R_LevelsSymbol));
    int present = presentLevelCount(f);
    if (present <= 1) {
      return as<NumericVector>(x);
    }
    int n = f.size();
    int cols = std::max(0, L - 1);
    NumericMatrix out(n, cols);
    std::fill(out.begin(), out.end(), 0.0);
    for (int i = 0; i < n; ++i) {
      int k = f[i]; // 1..L or NA
      if (k != NA_INTEGER && k > 1) out(i, k - 2) = 1.0;
    }
    CharacterVector lev = Rf_getAttrib(x, R_LevelsSymbol);
    if (cols > 0) {
      CharacterVector cn(cols);
      for (int j = 0; j < cols; ++j) cn[j] = lev[j + 1];
      setColNames(out, cn);
    }
    return out;
  } else {
    return as<NumericVector>(x);
  }
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
    if (present <= 1) {
      return as<NumericVector>(x);
    }
    int n = f.size();
    NumericMatrix out(n, L);
    std::fill(out.begin(), out.end(), 0.0);
    for (int i = 0; i < n; ++i) {
      int k = f[i]; // 1..L or NA
      if (k != NA_INTEGER) out(i, k - 1) = 1.0;
    }
    CharacterVector lev = Rf_getAttrib(x, R_LevelsSymbol);
    setColNames(out, lev);
    return out;
  } else {
    return as<NumericVector>(x);
  }
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
    int start = (n % lag) + 1; // 1-based
    std::vector<double> ser;
    for (int pos = start; pos <= n; pos += lag) {
      ser.push_back(x[pos - 1]);
    }
    int m = (int)ser.size();
    NumericVector s(m);
    IntegerVector idx(m);
    for (int i = 0; i < m; ++i) { s[i] = ser[i]; idx[i] = i + 1; }
    comp_series[t] = s;
    comp_index[t] = idx;
  }
  return List::create(_["Component.index"] = comp_index,
                      _["Component.series"] = comp_series);
}

// ---------- 5) generate.lin.vectors (+ recycled-list helpers) ----------

static List create_recycled_list_int(const IntegerVector& values, int list_length) {
  List res(list_length);
  std::vector< std::vector<int> > buckets(list_length);
  for (int i = 0; i < values.size(); ++i) {
    int index = (i % list_length); // 0-based
    buckets[index].push_back(values[i]);
  }
  for (int j = 0; j < list_length; ++j) {
    IntegerVector v(buckets[j].begin(), buckets[j].end());
    if (v.size() == 0) res[j] = R_NilValue; else res[j] = v;
  }
  return res;
}

static List create_recycled_list_num(const NumericVector& values, int list_length) {
  List res(list_length);
  std::vector< std::vector<double> > buckets(list_length);
  for (int i = 0; i < values.size(); ++i) {
    int index = (i % list_length); // 0-based
    buckets[index].push_back(values[i]);
  }
  for (int j = 0; j < list_length; ++j) {
    NumericVector v(buckets[j].begin(), buckets[j].end());
    if (v.size() == 0) res[j] = R_NilValue; else res[j] = v;
  }
  return res;
}

// [[Rcpp::export(name = "generate.lin.vectors")]]
List generate_lin_vectors(NumericVector x, int l, int h = 1) {
  int n = x.size();
  int max_fcast = std::min(h, l);
  
  // Build Component.series for i = 1..max_fcast
  List comp_series(max_fcast), comp_index(max_fcast);
  for (int i = 1; i <= max_fcast; ++i) {
    int start = ((n + i - 1) % l) + 1; // 1-based
    std::vector<double> ser;
    for (int pos = start; pos <= n; pos += l) ser.push_back(x[pos - 1]);
    int m = (int)ser.size();
    NumericVector s(m);
    IntegerVector idx(m);
    for (int k = 0; k < m; ++k) { s[k] = ser[k]; idx[k] = k + 1; }
    comp_series[i - 1] = s;
    comp_index[i - 1] = idx;
  }
  
  // forecast.index: recycled buckets of 1..h over max_fcast lists
  IntegerVector one_to_h(h);
  for (int i = 0; i < h; ++i) one_to_h[i] = i + 1;
  List forecast_index = create_recycled_list_int(one_to_h, max_fcast);
  
  // forecast.values.raw and forecast.values (see R logic)
  std::vector<double> fvals_raw(h);
  for (int i = 1; i <= h; ++i) {
    int recycled_index = ((i - 1) % l) + 1;
    int ci = ((recycled_index - 1) % std::max(1, max_fcast)) + 1;
    IntegerVector last_idx = comp_index[ci - 1];
    int last_val = last_idx.size();
    int forecast_increment = (int)std::ceil((double)i / (double)l);
    fvals_raw[i - 1] = (double)last_val + (double)forecast_increment;
  }
  NumericVector fvals = wrap(fvals_raw);
  List forecast_values = create_recycled_list_num(fvals, l);
  
  return List::create(_["Component.index"] = comp_index,
                      _["Component.series"] = comp_series,
                      _["forecast.values"] = forecast_values,
                      _["forecast.index"] = forecast_index);
}

// ---------- 6) ARMA.seas.weighting ----------

// [[Rcpp::export(name = "ARMA.seas.weighting")]]
List ARMA_seas_weighting(bool sf, SEXP mat) {
  // Scalar/vector => lag=M[1], Weights=1
  if (!Rf_isMatrix(mat) && !(Rf_inherits(mat, "data.frame")) && TYPEOF(mat) != VECSXP) {
    NumericVector M = as<NumericVector>(mat);
    double lag = (M.size() > 0 && !NumericVector::is_na(M[0])) ? M[0] : NA_REAL;
    return List::create(_["lag"] = lag, _["Weights"] = 1.0);
  }
  
  // compute ncol(mat)
  int n = NA_INTEGER;
  if (Rf_isMatrix(mat)) {
    IntegerVector dims = Rf_getAttrib(mat, R_DimSymbol);
    if (dims.size() == 2) n = dims[1];
  } else if (Rf_inherits(mat, "data.frame")) {
    List DF(mat);
    n = DF.size();
  } else if (TYPEOF(mat) == VECSXP) {
    List L(mat);
    n = L.size();
  }
  
  if (n == 1) {
    return List::create(_["lag"] = 1, _["Weights"] = 1.0);
  }
  
  // general case: n > 1, mat is a list/data.frame with named components
  List M(mat);
  
  if (sf) {
    // lag <- M$all.periods$Period[1]
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
  
  // Determine lag from seasonality test
  NumericVector lag = flattenNumericNoNA(M.containsElementNamed("Period") ? M["Period"] : R_NilValue);
  NumericVector observation_weighting(lag.size());
  for (int i = 0; i < lag.size(); ++i) observation_weighting[i] = 1.0 / std::sqrt(lag[i]);
  
  NumericVector covar = flattenNumericNoNA(M.containsElementNamed("Coefficient.of.Variation") ? M["Coefficient.of.Variation"] : R_NilValue);
  NumericVector varcovar = flattenNumericNoNA(M.containsElementNamed("Variable.Coefficient.of.Variation") ? M["Variable.Coefficient.of.Variation"] : R_NilValue);
  
  NumericVector lag_weighting;
  if (covar.size() == 1 && NumericVector::is_na(covar[0])) {
    lag_weighting = NumericVector(varcovar.size(), 1.0);
  } else {
    int m = std::min(varcovar.size(), covar.size());
    lag_weighting = NumericVector(m);
    for (int i = 0; i < m; ++i) lag_weighting[i] = varcovar[i] - covar[i];
    observation_weighting = observation_weighting[Rcpp::Range(0, m - 1)];
  }
  
  NumericVector wprod(lag_weighting.size());
  double denom = 0.0;
  for (int i = 0; i < wprod.size() && i < observation_weighting.size(); ++i) {
    wprod[i] = lag_weighting[i] * observation_weighting[i];
    denom += wprod[i];
  }
  NumericVector Weights = (denom == 0.0) ? NumericVector(wprod.size(), 0.0) : wprod / denom;
  
  return List::create(_["lag"] = lag, _["Weights"] = Weights);
}

// ---------- 7) lag.mtx ----------

// [[Rcpp::export(name = "lag.mtx")]]
SEXP lag_mtx(SEXP x, SEXP tau) {
  NumericMatrix X;
  CharacterVector colheads;
  bool x_was_vector = false;
  
  if (Rf_isMatrix(x)) {
    X = as<NumericMatrix>(x);
  } else if (Rf_inherits(x, "data.frame")) {
    DataFrame DF = as<DataFrame>(x);
    int n = DF.nrows(), p = DF.size();
    X = NumericMatrix(n, p);
    for (int j = 0; j < p; ++j) X(_, j) = as<NumericVector>(DF[j]);
    colheads = DF.names();
  } else {
    // vector -> single-column matrix
    NumericVector v = as<NumericVector>(x);
    X = NumericMatrix(v.size(), 1);
    X(_, 0) = v;
    x_was_vector = true;
  }
  
  int n = X.nrow();
  int p = X.ncol();
  
  // compute max_tau and normalize tau into a List-of-numeric
  List tau_list;
  NumericVector tau_all;
  
  if (TYPEOF(tau) == VECSXP) {
    tau_list = List(tau);
    for (int i = 0; i < tau_list.size(); ++i) {
      NumericVector v = as<NumericVector>(tau_list[i]);
      for (double z : v) tau_all.push_back(z);
    }
  } else {
    NumericVector v = as<NumericVector>(tau);
    tau_list = List(p);
    for (int j = 0; j < p; ++j) tau_list[j] = v; // same set for each column
    tau_all = v;
  }
  
  int max_tau = 0;
  for (double z : tau_all) if (!NumericVector::is_na(z)) max_tau = std::max(max_tau, (int)std::round(z));
  
  // Build per-column "0..max_tau" lag blocks
  std::vector< NumericMatrix > blocks;
  blocks.reserve(p);
  std::vector< CharacterVector > block_names;
  
  for (int j = 0; j < p; ++j) {
    int rows = n - max_tau;
    if (rows < 0) rows = 0;
    NumericMatrix B(rows, max_tau + 1);
    CharacterVector names(max_tau + 1);
    std::string base;
    if (x_was_vector) {
      base = "x";
    } else if (colheads.size() > 0) {
      base = as<std::string>(colheads[j]);
    } else {
      base = "x";
    }
    for (int k = 0; k <= max_tau; ++k) {
      names[k] = base + "_tau_" + std::to_string(k);
      int start = max_tau - k;         // 0-based row start
      for (int i = 0; i < rows; ++i) { // rows = n - max_tau
        B(i, k) = X(i + start, j);
      }
    }
    setColNames(B, names);
    blocks.push_back(B);
    block_names.push_back(names);
  }
  
  // Combine blocks column-wise
  int rows = (n - max_tau < 0) ? 0 : n - max_tau;
  int total_cols = p * (max_tau + 1);
  NumericMatrix M(rows, total_cols);
  CharacterVector Mnames(total_cols);
  int pos = 0;
  for (int j = 0; j < p; ++j) {
    NumericMatrix& B = blocks[j];
    CharacterVector nm = block_names[j];
    for (int k = 0; k <= max_tau; ++k) {
      M(_, pos) = B(_, k);
      Mnames[pos] = nm[k];
      ++pos;
    }
  }
  setColNames(M, Mnames);
  
  // If length(unlist(tau)) > 1: select relevant lags
  if (tau_all.size() > 1) {
    std::vector<int> keep;
    for (int i = 0; i < p; ++i) {
      NumericVector tvec = as<NumericVector>(tau_list[i]);
      int base = i * (max_tau + 1);
      keep.push_back(base + 0); // always keep tau_0
      for (int u = 0; u < tvec.size(); ++u) {
        int tk = (int)std::round(tvec[u]);
        if (tk >= 0 && tk <= max_tau) keep.push_back(base + tk);
      }
    }
    std::sort(keep.begin(), keep.end());
    // build reduced matrix
    NumericMatrix Rm(rows, keep.size());
    CharacterVector Rn(keep.size());
    for (int k = 0; k < (int)keep.size(); ++k) {
      Rm(_, k) = M(_, keep[k]);
      Rn[k] = Mnames[keep[k]];
    }
    setColNames(Rm, Rn);
    M = Rm;
    Mnames = Rn;
  }
  
  // Reorder so that "tau_0" columns come first
  std::vector<int> vars, others;
  for (int j = 0; j < M.ncol(); ++j) {
    std::string nm = as<std::string>(Mnames[j]);
    if (nm.find("tau_0") != std::string::npos) vars.push_back(j);
    else others.push_back(j);
  }
  std::vector<int> order; order.reserve(M.ncol());
  order.insert(order.end(), vars.begin(), vars.end());
  order.insert(order.end(), others.begin(), others.end());
  
  NumericMatrix Out(rows, M.ncol());
  CharacterVector Onames(M.ncol());
  for (int jj = 0; jj < (int)order.size(); ++jj) {
    Out(_, jj) = M(_, order[jj]);
    Onames[jj] = Mnames[order[jj]];
  }
  setColNames(Out, Onames);
  
  // Return as data.frame, like the R function
  Out.attr("class") = CharacterVector::create("data.frame");
  Out.attr("row.names") = IntegerVector::create(NA_INTEGER, -(rows));
  return Out;
}

// ---------- 8) NNS.meboot.part ----------

// [[Rcpp::export(name = "NNS.meboot.part")]]
NumericVector NNS_meboot_part(NumericVector x, int n, NumericVector z,
                              double xmin, double xmax,
                              NumericVector desintxb, bool reachbnd) {
  NumericVector p = runif(n);
  Function quantile("quantile");
  NumericVector q = as<NumericVector>(quantile(x, _["probs"] = p));
  
  std::vector<int> ref1;
  double invn = 1.0 / (double)n;
  for (int i = 0; i < p.size(); ++i) if (p[i] <= invn) ref1.push_back(i);
  if (!ref1.empty()) {
    NumericVector px(ref1.size());
    for (int i = 0; i < (int)ref1.size(); ++i) px[i] = p[ref1[i]];
    Function approx("approx");
    NumericVector x0 = NumericVector::create(0.0, invn);
    NumericVector y0 = NumericVector::create(xmin, z[0]);
    List aout = approx(_["x"] = x0, _["y"] = y0, _["xout"] = px);
    NumericVector qq = aout["y"];
    for (int i = 0; i < (int)ref1.size(); ++i) {
      double val = qq[i];
      if (!reachbnd) val = qq[i] + desintxb[0] - 0.5 * (z[0] + xmin);
      q[ref1[i]] = val;
    }
  }
  
  std::vector<int> ref4;
  double edge = (double)(n - 1) / (double)n;
  for (int i = 0; i < p.size(); ++i) if (p[i] == edge) ref4.push_back(i);
  for (int idx : ref4) q[idx] = z[n - 2];
  
  std::vector<int> ref5;
  for (int i = 0; i < p.size(); ++i) if (p[i] > edge) ref5.push_back(i);
  if (!ref5.empty()) {
    NumericVector px(ref5.size());
    for (int i = 0; i < (int)ref5.size(); ++i) px[i] = p[ref5[i]];
    Function approx("approx");
    NumericVector x1 = NumericVector::create(edge, 1.0);
    NumericVector y1 = NumericVector::create(z[n - 2], xmax);
    List aout = approx(_["x"] = x1, _["y"] = y1, _["xout"] = px);
    NumericVector qq = aout["y"];
    for (int i = 0; i < (int)ref5.size(); ++i) {
      double val = qq[i];
      if (!reachbnd) val = qq[i] + desintxb[n - 1] - 0.5 * (z[n - 2] + xmax);
      q[ref5[i]] = val;
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
      int nr = DF.nrows(), p = DF.size();
      X = NumericMatrix(nr, p);
      for (int j = 0; j < p; ++j) X(_, j) = as<NumericVector>(DF[j]);
    } else {
      X = as<NumericMatrix>(x);
    }
    sdx = colSd(X);
  } else {
    sdx = NumericVector::create(vecSd(as<NumericVector>(x)));
  }
  
  NumericVector ens_sd = colSd(ensemble);
  NumericVector sdf(sdx.size() + ens_sd.size());
  int pos = 0;
  for (double v : sdx) sdf[pos++] = v;
  for (double v : ens_sd) sdf[pos++] = v;
  
  NumericVector sdfa = clone(sdf); // actual / original
  for (int i = 0; i < sdfa.size(); ++i) sdfa[i] = sdf[i] / sdf[0];
  NumericVector sdfd = clone(sdf); // desired / actual
  for (int i = 0; i < sdfd.size(); ++i) sdfd[i] = sdf[0] / sdf[i];
  
  double mx = 1.0 + (fiv / 100.0);
  for (int i = 0; i < sdfa.size(); ++i) {
    if (sdfa[i] < 1.0) {
      sdfa[i] = R::runif(1.0, mx);
    }
  }
  
  int J = ensemble.ncol();
  NumericVector sdfdXsdfa(J);
  for (int j = 0; j < J; ++j) sdfdXsdfa[j] = sdfd[j + 1] * sdfa[j + 1];
  
  for (int j = 0; j < J; ++j) {
    if (std::floor(sdfdXsdfa[j]) > 0.0) {
      double a = sdfdXsdfa[j];
      for (int i = 0; i < ensemble.nrow(); ++i) ensemble(i, j) *= a;
    }
  }
  
  Function is_ts("is.ts"), frequency("frequency"), startF("start"), tsF("ts");
  LogicalVector its = is_ts(x);
  if (its.size() && its[0]) {
    SEXP freq = frequency(x);
    SEXP start0 = startF(x);
    return tsF(ensemble, _["frequency"] = freq, _["start"] = start0);
  }
  return ensemble;
}

// ---------- 10) force.clt ----------

// [[Rcpp::export(name = "force.clt")]]
SEXP force_clt(SEXP x, NumericMatrix ensemble) {
  int n = ensemble.nrow();
  int bigj = ensemble.ncol();
  
  // desired grand mean gm and s (sd of x)
  double gm = NA_REAL;
  NumericVector s;
  
  if (Rf_isMatrix(x) || Rf_inherits(x, "data.frame")) {
    NumericMatrix X;
    if (Rf_inherits(x, "data.frame")) {
      DataFrame DF = as<DataFrame>(x);
      int nr = DF.nrows(), p = DF.size();
      X = NumericMatrix(nr, p);
      for (int j = 0; j < p; ++j) X(_, j) = as<NumericVector>(DF[j]);
    } else {
      X = as<NumericMatrix>(x);
    }
    // gm: mean of all entries in X
    double sumAll = 0.0;
    for (int i = 0; i < X.nrow(); ++i)
      for (int j = 0; j < X.ncol(); ++j)
        sumAll += X(i, j);
    gm = sumAll / (double)(X.nrow() * X.ncol());
    s = colSd(X); // vector
  } else {
    NumericVector xv = as<NumericVector>(x);
    double sumAll = 0.0;
    for (int i = 0; i < xv.size(); ++i) sumAll += xv[i];
    gm = sumAll / (double)xv.size();
    s = NumericVector::create(vecSd(xv));
  }
  
  // xbar: column means of ensemble
  NumericVector xbar(bigj);
  for (int j = 0; j < bigj; ++j) {
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += ensemble(i, j);
    xbar[j] = mu / (double)n;
  }
  
  // order indices by xbar
  IntegerVector oo(bigj);
  for (int j = 0; j < bigj; ++j) oo[j] = j;
  std::sort(oo.begin(), oo.end(), [&](int a, int b){ return xbar[a] < xbar[b]; });
  
  // sortxbar for difference calc
  NumericVector sortxbar = clone(xbar);
  std::sort(sortxbar.begin(), sortxbar.end());
  
  // smean = s / sqrt(bigj)  (recycle if s is vector)
  NumericVector smean = clone(s);
  for (int i = 0; i < smean.size(); ++i) smean[i] = s[i] / std::sqrt((double)bigj);
  double smean_scalar = smean.size() ? smean[0] : 0.0;
  
  // newbar = gm + qnorm(1:bigj/(bigj+1)) * smean
  NumericVector qn(bigj);
  for (int j = 0; j < bigj; ++j) qn[j] = R::qnorm((double)(j + 1) / (double)(bigj + 1), 0.0, 1.0, 1, 0);
  NumericVector newbar(bigj);
  for (int j = 0; j < bigj; ++j) {
    double sm = (smean.size() == 1) ? smean_scalar : smean[j % smean.size()];
    newbar[j] = gm + qn[j] * sm;
  }
  
  // scale newbar to zero-mean / unit-sd
  double mu_nb = 0.0;
  for (int j = 0; j < bigj; ++j) mu_nb += newbar[j];
  mu_nb /= (double)bigj;
  double ss_nb = 0.0;
  for (int j = 0; j < bigj; ++j) { double d = newbar[j] - mu_nb; ss_nb += d * d; }
  double sd_nb = std::sqrt(ss_nb / (double)(bigj - 1));
  NumericVector scn(bigj);
  for (int j = 0; j < bigj; ++j) scn[j] = (newbar[j] - mu_nb) / sd_nb;
  
  // newm = scn * smean + gm
  NumericVector newm(bigj);
  for (int j = 0; j < bigj; ++j) {
    double sm = (smean.size() == 1) ? smean_scalar : smean[j % smean.size()];
    newm[j] = scn[j] * sm + gm;
  }
  
  // meanfix = newm - sortxbar
  NumericVector meanfix(bigj);
  for (int j = 0; j < bigj; ++j) meanfix[j] = newm[j] - sortxbar[j];
  
  // apply fixes in the sorted order
  NumericMatrix out = clone(ensemble);
  for (int i = 0; i < bigj; ++i) {
    int col = oo[i];
    double add = meanfix[i];
    for (int r = 0; r < n; ++r) out(r, col) = ensemble(r, col) + add;
  }
  
  // preserve ts attributes if ensemble is ts
  Function is_ts("is.ts"), frequency("frequency"), startF("start"), tsF("ts");
  LogicalVector its = is_ts(ensemble);
  if (its.size() && its[0]) {
    SEXP freq = frequency(ensemble);
    SEXP start0 = startF(ensemble);
    return tsF(out, _["frequency"] = freq, _["start"] = start0);
  }
  return out;
}

// ---------- 11) downSample / upSample ----------

// [[Rcpp::export]]
SEXP downSample(SEXP x, SEXP y, bool list = false, std::string yname = "Class") {
  // x -> data.frame(stringsAsFactors = TRUE)
  Function as_df("as.data.frame");
  DataFrame X = as<DataFrame>(as_df(x, _["stringsAsFactors"] = true));
  
  // y must be factor
  Function is_factor("is.factor");
  LogicalVector ok = is_factor(y);
  if (!(ok.size() && ok[0])) {
    Rcpp::warning("Down-sampling requires a factor variable as the response. The original data was returned.");
    return List::create(_["x"] = X, _["y"] = y);
  }
  
  IntegerVector fy = as<IntegerVector>(y); // factor codes 1..L
  CharacterVector lev = Rf_getAttrib(y, R_LevelsSymbol);
  int L = lev.size();
  
  // class counts (ignore NA)
  std::vector< IntegerVector > perClass(L);
  for (int i = 0; i < fy.size(); ++i) {
    int k = fy[i];
    if (k == NA_INTEGER) continue;
    perClass[k - 1].push_back(i + 1); // 1-based indices for R subsetting
  }
  int minClass = INT_MAX;
  for (int k = 0; k < L; ++k) if (perClass[k].size() > 0) minClass = std::min(minClass, (int)perClass[k].size());
  if (minClass == INT_MAX) minClass = 0;
  
  // sample each class
  Function rbind_("rbind"), cbind_("cbind");
  List sampled_frames;
  CharacterVector presentLevels;
  for (int k = 0; k < L; ++k) {
    if (perClass[k].size() == 0) continue;
    IntegerVector idx = sampleWithoutReplacement(perClass[k], minClass);
    Language subset_call("[", X, idx, R_MissingArg); // rows only; drop not needed
    SEXP Xi = subset_call.eval();
    sampled_frames.push_back(Xi);
    presentLevels.push_back(lev[k]);
  }
  
  // combine
  SEXP Xc;
  if (sampled_frames.size() == 0) {
    Xc = X;
  } else if (sampled_frames.size() == 1) {
    Xc = sampled_frames[0];
  } else {
    Xc = rbind_(sampled_frames);
  }
  
  // new y vector
  int rows = Rf_length(Rf_getAttrib(Xc, R_RowNamesSymbol));
  if (rows < 0) rows = -rows; // compact row.names negative encoding
  IntegerVector ynew(rows);
  CharacterVector ylev = presentLevels;
  int block = (minClass > 0) ? minClass : 0;
  int pos = 0;
  for (int k = 0; k < ylev.size(); ++k) {
    for (int i = 0; i < block; ++i) ynew[pos++] = k + 1;
  }
  ynew.attr("class") = "factor";
  ynew.attr("levels") = ylev;
  
  if (list) {
    return List::create(_["x"] = Xc, _["y"] = ynew);
  } else {
    SEXP out = cbind_(Xc, ynew);
    List DF(out);
    CharacterVector names = DF.names();
    names[names.size() - 1] = yname;
    DF.names() = names;
    return DF;
  }
}

// [[Rcpp::export]]
SEXP upSample(SEXP x, SEXP y, bool list = false, std::string yname = "Class") {
  // x -> data.frame(stringsAsFactors = TRUE)
  Function as_df("as.data.frame");
  DataFrame X = as<DataFrame>(as_df(x, _["stringsAsFactors"] = true));
  
  // y must be factor
  Function is_factor("is.factor");
  LogicalVector ok = is_factor(y);
  if (!(ok.size() && ok[0])) {
    Rcpp::warning("Up-sampling requires a factor variable as the response. The original data was returned.");
    return List::create(_["x"] = X, _["y"] = y);
  }
  
  IntegerVector fy = as<IntegerVector>(y); // factor codes 1..L
  CharacterVector lev = Rf_getAttrib(y, R_LevelsSymbol);
  int L = lev.size();
  
  // per-class indices
  std::vector< IntegerVector > perClass(L);
  for (int i = 0; i < fy.size(); ++i) {
    int k = fy[i];
    if (k == NA_INTEGER) continue;
    perClass[k - 1].push_back(i + 1);
  }
  int maxClass = 0;
  for (int k = 0; k < L; ++k) maxClass = std::max(maxClass, (int)perClass[k].size());
  
  // assemble upsampled frames
  Function rbind_("rbind"), cbind_("cbind");
  List sampled_frames;
  for (int k = 0; k < L; ++k) {
    IntegerVector idx = perClass[k];
    if (idx.size() == 0) continue;
    if (idx.size() < maxClass) {
      IntegerVector extra(maxClass - idx.size());
      for (int i = 0; i < extra.size(); ++i) {
        int pick = (int)std::floor(R::runif(0.0, 1.0) * idx.size()); // 0..size-1
        extra[i] = idx[pick];
      }
      IntegerVector idx2(idx.size() + extra.size());
      std::copy(idx.begin(), idx.end(), idx2.begin());
      std::copy(extra.begin(), extra.end(), idx2.begin() + idx.size());
      idx = idx2;
    }
    Language subset_call("[", X, idx, R_MissingArg); // rows only; drop not needed
    SEXP Xi = subset_call.eval();
    sampled_frames.push_back(Xi);
  }
  
  // combine
  SEXP Xc;
  if (sampled_frames.size() == 0) {
    Xc = X;
  } else if (sampled_frames.size() == 1) {
    Xc = sampled_frames[0];
  } else {
    Xc = rbind_(sampled_frames);
  }
  
  // new y vector: repeat each level maxClass times
  int rows = Rf_length(Rf_getAttrib(Xc, R_RowNamesSymbol));
  if (rows < 0) rows = -rows;
  IntegerVector ynew(rows);
  int pos = 0;
  for (int k = 0; k < L; ++k) {
    if (perClass[k].size() == 0) continue;
    for (int i = 0; i < maxClass; ++i) ynew[pos++] = k + 1;
  }
  ynew.attr("class") = "factor";
  ynew.attr("levels") = lev;
  
  if (list) {
    return List::create(_["x"] = Xc, _["y"] = ynew);
  } else {
    SEXP out = cbind_(Xc, ynew);
    List DF(out);
    CharacterVector names = DF.names();
    names[names.size() - 1] = yname;
    DF.names() = names;
    return DF;
  }
}
