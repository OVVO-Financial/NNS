// src/internal_functions.cpp
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>

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
  //Rf_setAttrib(m, R_DimNamesSymbol, List::create(R_NilValue, nm));
  colnames(m) = nm;
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



// ---------- 8) NNS.meboot.part ----------

// [[Rcpp::export(name = "NNS.meboot.part")]]
NumericVector NNS_meboot_part(NumericVector xx, int n, NumericVector z,
                              double xmin, double xmax,
                              NumericVector desintxb, bool reachbnd) {
  NumericVector p = runif(n);
  int m = xx.size();
  NumericVector q(n);

  if (m == 0) {
    std::fill(q.begin(), q.end(), NA_REAL);
  } else if (m == 1) {
    std::fill(q.begin(), q.end(), xx[0]);
  } else {
    for (int i = 0; i < n; ++i) {
      double pi = p[i];
      if (pi <= 0.0) {
        q[i] = xx[0];
        continue;
      }
      if (pi >= 1.0) {
        q[i] = xx[m - 1];
        continue;
      }

      double h = 1.0 + (m - 1.0) * pi;
      int j = static_cast<int>(std::floor(h));
      double g = h - static_cast<double>(j);

      if (j < 1) j = 1;
      if (j > m - 1) j = m - 1;

      q[i] = (1.0 - g) * xx[j - 1] + g * xx[j];
    }
  }
  
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
  
  double edge = (double)(n - 1) / (double)n;
  std::vector<int> ref5;
  for (int i = 0; i < p.size(); ++i) if (p[i] >= edge) ref5.push_back(i);
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
      //int nr = DF.nrows(), p = DF.size();
      int nr = (Rf_length(x) == 0) ? 0 : Rf_length(VECTOR_ELT(x, 0));
      int p  = Rf_length(x);
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
      //int nr = DF.nrows(), p = DF.size();
      int nr = (Rf_length(x) == 0) ? 0 : Rf_length(VECTOR_ELT(x, 0));
      int p  = Rf_length(x);
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

// ---- helpers ---------------------------------------------------------------

// copy rows for a factor (IntegerVector with class "factor" and "levels")
static IntegerVector subset_factor_codes(const IntegerVector& codes,
                                         const IntegerVector& rows) {
  const int m = rows.size();
  IntegerVector out(m);
  for (int i = 0; i < m; ++i) {
    const int idx = rows[i] - 1; // rows are 1-based
    out[i] = (idx >= 0 && idx < codes.size()) ? codes[idx] : NA_INTEGER;
  }
  return out;
}

// generic: copy rows for a simple vector while preserving NA and type
template <int RTYPE>
static Rcpp::Vector<RTYPE> subset_vec_template(const Rcpp::Vector<RTYPE>& v,
                                               const IntegerVector& rows) {
  const int m = rows.size();
  Rcpp::Vector<RTYPE> out(m);
  for (int i = 0; i < m; ++i) {
    const int idx = rows[i] - 1;
    if (idx >= 0 && idx < v.size()) {
      out[i] = v[idx];
    } else {
      out[i] = Rcpp::Vector<RTYPE>::get_na();   // works for INT/REAL/LGL
    }
  }
  return out;
}

// specialization for character vectors (STRSXP) to avoid proxy/SEXP mismatch
template <>
inline Rcpp::CharacterVector
subset_vec_template<STRSXP>(const Rcpp::CharacterVector& v,
                            const IntegerVector& rows) {
  const int m = rows.size();
  Rcpp::CharacterVector out(m);
  for (int i = 0; i < m; ++i) {
    const int idx = rows[i] - 1;
    if (idx >= 0 && idx < v.size()) {
      out[i] = v[idx];
    } else {
      out[i] = NA_STRING;                      // explicit NA for character
    }
  }
  return out;
}

// build a data.frame from X for row indices `rows`,
// optionally appending the factor `y_factor` as last column named `yname`.
static DataFrame subset_df_rows_with_y(const DataFrame& X,
                                       const IntegerVector& rows,
                                       SEXP y_factor,                 // factor (or R_NilValue)
                                       const std::string& yname,
                                       bool include_y) {
  const int p = X.size();
  const int m = rows.size();
  
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
      // if factor, carry class + levels over
      if (!cls.isNULL() && as<CharacterVector>(cls).size() > 0 &&
          as<CharacterVector>(cls)[0] == "factor") {
        IntegerVector sub = subset_factor_codes(iv, rows);
        sub.attr("class")  = iv.attr("class");
        sub.attr("levels") = iv.attr("levels");
        out[j] = sub;
      } else {
        out[j] = subset_vec_template<INTSXP>(iv, rows);
      }
      break;
    }
    case REALSXP: out[j] = subset_vec_template<REALSXP>(NumericVector(col), rows); break;
    case LGLSXP:  out[j] = subset_vec_template<LGLSXP>(LogicalVector(col), rows);  break;
    case STRSXP:  out[j] = subset_vec_template<STRSXP>(CharacterVector(col), rows);break;
    default: {
      // fallback: coerce exotic types to character
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
    out[p]       = ysub;
    out_names[p] = yname;
  }
  
  out.attr("names")     = out_names;
  out.attr("class")     = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -m);
  return DataFrame(out);
}

// sample k integers in 1..N (with/without replacement) using R RNG
static IntegerVector sample_indices(int N, int k, bool replace) {
  IntegerVector res(k);
  if (N <= 0 || k <= 0) return res;
  if (!replace && k > N) k = N;
  
  if (replace) {
    for (int i = 0; i < k; ++i) {
      int draw = 1 + (int)floor(R::runif(0.0, 1.0) * N);
      if (draw > N) draw = N;
      res[i] = draw;
    }
  } else {
    // partial Fisher–Yates for first k positions
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
  
  // Coerce x -> data.frame (stringsAsFactors = TRUE equivalent)
  Function as_df("as.data.frame");
  DataFrame X = as<DataFrame>(as_df(x));
  
  // Require factor y (caret warns & returns original when not a factor)
  Function is_factor("is.factor");
  if (!as<bool>(is_factor(y))) {
    Rcpp::warning("Down-sampling requires a factor variable as the response. The original data was returned.");
    return List::create(_["x"] = X, _["y"] = y);
  }
  
  // y as factor codes (1..L) and levels
  IntegerVector fy = as<IntegerVector>(y);
  CharacterVector lev = Rf_getAttrib(y, R_LevelsSymbol);
  //const int n = X.nrows();
  const int n = Rf_nrows(X);
  if (fy.size() != n) stop("downSample: nrow(x) != length(y)");
  const int L = lev.size();
  
  // collect row indices per class (1-based)
  std::vector< std::vector<int> > perClass(L);
  for (int i = 0; i < n; ++i) {
    int k = fy[i];
    if (k != NA_INTEGER) perClass[k - 1].push_back(i + 1);
  }
  
  // target size: min class count (caret::min(table(y)))
  int minClass = n;
  bool any_ok = false;
  for (int k = 0; k < L; ++k) {
    int sz = (int)perClass[k].size();
    if (sz > 0) { any_ok = true; if (sz < minClass) minClass = sz; }
  }
  if (!any_ok || minClass <= 0) stop("downSample: no non-empty classes.");
  
  // sample minClass rows within each class without replacement
  std::vector<int> rows_out;
  rows_out.reserve(minClass * L);
  for (int k = 0; k < L; ++k) {
    const int gsz = (int)perClass[k].size();
    if (gsz == 0) continue;
    IntegerVector s = sample_indices(gsz, minClass, /*replace*/false); // 1..gsz
    for (int j = 0; j < s.size(); ++j) rows_out.push_back(perClass[k][ s[j] - 1 ]);
  }
  IntegerVector rows = wrap(rows_out);
  
  if (list) {
    DataFrame Xout = subset_df_rows_with_y(X, rows, R_NilValue, yname, /*include_y=*/false);
    IntegerVector ysub = subset_factor_codes(fy, rows);
    ysub.attr("class")  = CharacterVector::create("factor");
    ysub.attr("levels") = lev;
    return List::create(_["x"] = Xout, _["y"] = ysub);
  } else {
    return subset_df_rows_with_y(X, rows, y, yname, /*include_y=*/true);
  }
}

// ---------- upSample --------------------------------------

// [[Rcpp::export]]
SEXP upSample(SEXP x, SEXP y, bool list = false, std::string yname = "Class") {
  RNGScope scope;
  
  // Coerce x -> data.frame
  Function as_df("as.data.frame");
  DataFrame X = as<DataFrame>(as_df(x));
  
  // Require factor y
  Function is_factor("is.factor");
  if (!as<bool>(is_factor(y))) {
    Rcpp::warning("Up-sampling requires a factor variable as the response. The original data was returned.");
    return List::create(_["x"] = X, _["y"] = y);
  }
  
  IntegerVector fy = as<IntegerVector>(y);
  CharacterVector lev = Rf_getAttrib(y, R_LevelsSymbol);
  //const int n = X.nrows();
  const int n = Rf_nrows(X);
  if (fy.size() != n) stop("upSample: nrow(x) != length(y)");
  const int L = lev.size();
  
  // collect row indices per class (1-based)
  std::vector< std::vector<int> > perClass(L);
  for (int i = 0; i < n; ++i) {
    int k = fy[i];
    if (k != NA_INTEGER) perClass[k - 1].push_back(i + 1);
  }
  
  // target size: max class count (caret::max(table(y)))
  int maxClass = 0;
  bool any_ok = false;
  for (int k = 0; k < L; ++k) {
    int sz = (int)perClass[k].size();
    if (sz > 0) { any_ok = true; if (sz > maxClass) maxClass = sz; }
  }
  if (!any_ok || maxClass <= 0) stop("upSample: no non-empty classes.");
  
  // sample up to maxClass within each class with replacement
  std::vector<int> rows_out;
  rows_out.reserve(maxClass * L);
  for (int k = 0; k < L; ++k) {
    const int gsz = (int)perClass[k].size();
    if (gsz == 0) continue;
    IntegerVector s = sample_indices(gsz, maxClass, /*replace*/true); // 1..gsz
    for (int j = 0; j < s.size(); ++j) rows_out.push_back(perClass[k][ s[j] - 1 ]);
  }
  IntegerVector rows = wrap(rows_out);
  
  if (list) {
    DataFrame Xout = subset_df_rows_with_y(X, rows, R_NilValue, yname, /*include_y=*/false);
    IntegerVector ysub = subset_factor_codes(fy, rows);
    ysub.attr("class")  = CharacterVector::create("factor");
    ysub.attr("levels") = lev;
    return List::create(_["x"] = Xout, _["y"] = ysub);
  } else {
    return subset_df_rows_with_y(X, rows, y, yname, /*include_y=*/true);
  }
}
