#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Small numerical tolerance for positive-definite checks
static inline bool is_pos(double x) { return x > 0.0 && std::isfinite(x); }

// [[Rcpp::export]]
List fast_lm(NumericVector x, NumericVector y) {
  int nx = x.size();
  int ny = y.size();
  
  if (nx != ny) {
    stop("fast_lm: length(x) != length(y) (got %i vs %i).", nx, ny);
  }
  
  // Means
  double mean_x = mean(x);
  double mean_y = mean(y);
  
  // Variance of x and covariance of (x,y)
  double var_x = 0.0, cov_xy = 0.0;
  for (int i = 0; i < nx; i++) {
    double dx = x[i] - mean_x;
    double dy = y[i] - mean_y;
    var_x += dx * dx;
    cov_xy += dx * dy;
  }
  
  NumericVector coef(2);
  NumericVector fitted(ny);
  NumericVector residuals(ny);
  
  if (var_x == 0.0) {
    // All x are identical -> slope = 0, intercept = mean(y)
    coef[0] = mean_y;   // intercept
    coef[1] = 0.0;      // slope
    
    for (int i = 0; i < ny; i++) {
      fitted[i] = mean_y;
      residuals[i] = y[i] - mean_y;
    }
    
  } else {
    // Standard OLS slope + intercept
    double slope = cov_xy / var_x;
    double intercept = mean_y - slope * mean_x;
    
    coef[0] = intercept;
    coef[1] = slope;
    
    for (int i = 0; i < ny; i++) {
      fitted[i] = intercept + slope * x[i];
      residuals[i] = y[i] - fitted[i];
    }
  }
  
  int df_resid = ny - 2;
  
  return List::create(
    Named("coef") = coef,
    Named("residuals") = residuals,
    Named("fitted.values") = fitted,
    Named("df.residual") = df_resid
  );
}

// --- Linear algebra helpers for multiple regression ---

// Cholesky decomposition of a symmetric positive-definite matrix A
// Returns lower triangular matrix L such that A = L * L^T
static NumericMatrix cholesky_decomposition(const NumericMatrix& A) {
  const R_xlen_t n = A.nrow();
  if (n != A.ncol()) stop("cholesky_decomposition: matrix must be square.");
  NumericMatrix L(n, n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    // Compute L(i, i)
    double sum = A(i, i);
    for (R_xlen_t k = 0; k < i; ++k) sum -= L(i, k) * L(i, k);
    if (!is_pos(sum)) stop("cholesky_decomposition: matrix not positive-definite (nonpositive pivot at %lld).", static_cast<long long>(i+1));
    L(i, i) = std::sqrt(sum);
    
    // Compute L(j, i) for j > i
    const double Lii = L(i, i);
    for (R_xlen_t j = i + 1; j < n; ++j) {
      double s = A(j, i);
      for (R_xlen_t k = 0; k < i; ++k) s -= L(j, k) * L(i, k);
      L(j, i) = s / Lii;
    }
  }
  return L;
}

// Solve L * z = b (forward substitution, L is lower triangular)
static NumericVector forward_substitution(const NumericMatrix& L, const NumericVector& b) {
  const R_xlen_t n = L.nrow();
  if (b.size() != n) stop("forward_substitution: incompatible dimensions.");
  NumericVector z(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    double sum = b[i];
    for (R_xlen_t j = 0; j < i; ++j) sum -= L(i, j) * z[j];
    const double Lii = L(i, i);
    if (Lii == 0.0 || !std::isfinite(Lii)) stop("forward_substitution: singular pivot.");
    z[i] = sum / Lii;
  }
  return z;
}

// Solve L^T * x = z (back substitution, L^T is upper triangular)
static NumericVector back_substitution(const NumericMatrix& L, const NumericVector& z) {
  const R_xlen_t n = L.nrow();
  if (z.size() != n) stop("back_substitution: incompatible dimensions.");
  NumericVector x(n);
  
  for (R_xlen_t i = n; i-- > 0; ) { // i = n-1 ... 0
    double sum = z[i];
    for (R_xlen_t j = i + 1; j < n; ++j) sum -= L(j, i) * x[j]; // L^T(i, j) = L(j, i)
    const double Lii = L(i, i);
    if (Lii == 0.0 || !std::isfinite(Lii)) stop("back_substitution: singular pivot.");
    x[i] = sum / Lii;
  }
  return x;
}

// [[Rcpp::export]]
List fast_lm_mult(NumericMatrix x, NumericVector y) {
  const R_xlen_t n = x.nrow();
  const R_xlen_t p = x.ncol();
  if (n == 0) stop("fast_lm_mult: 'x' has zero rows.");
  if (p == 0) stop("fast_lm_mult: 'x' has zero columns.");
  if (y.size() != n) stop("fast_lm_mult: length(y) != nrow(x) (got %lld vs %lld).",
      static_cast<long long>(y.size()), static_cast<long long>(n));
  
  // Design matrix with intercept
  NumericMatrix X(n, p + 1);
  for (R_xlen_t i = 0; i < n; ++i) {
    X(i, 0) = 1.0; // Intercept
    for (R_xlen_t j = 0; j < p; ++j) X(i, j + 1) = x(i, j);
  }
  
  // Compute XtX and Xty
  const R_xlen_t q = p + 1;
  NumericMatrix XtX(q, q);
  NumericVector Xty(q);
  
  for (R_xlen_t i = 0; i < q; ++i) {
    for (R_xlen_t j = 0; j <= i; ++j) { // fill lower triangle, then mirror
      double s = 0.0;
      for (R_xlen_t k = 0; k < n; ++k) s += X(k, i) * X(k, j);
      XtX(i, j) = s;
      if (i != j) XtX(j, i) = s;
    }
    double sy = 0.0;
    for (R_xlen_t k = 0; k < n; ++k) sy += X(k, i) * y[k];
    Xty[i] = sy;
  }
  
  // Solve normal equations via Cholesky
  NumericMatrix L = cholesky_decomposition(XtX);
  NumericVector z = forward_substitution(L, Xty);
  NumericVector coef = back_substitution(L, z);
  
  // Fitted values and residuals
  NumericVector fitted_values(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (R_xlen_t j = 0; j < q; ++j) s += coef[j] * X(i, j);
    fitted_values[i] = s;
  }
  NumericVector residuals = y - fitted_values;
  
  // R-squared
  const double y_mean = mean(y);
  double TSS = 0.0, RSS = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double dy = y[i] - y_mean;
    TSS += dy * dy;
    const double re = residuals[i];
    RSS += re * re;
  }
  const double R2 = (TSS == 0.0) ? NA_REAL : (1.0 - RSS / TSS);
  
  return List::create(
    _["coefficients"] = coef,
    _["fitted.values"] = fitted_values,
    _["residuals"]     = residuals,
    _["r.squared"]     = R2
  );
}
