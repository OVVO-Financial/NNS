#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List fast_lm(NumericVector x, NumericVector y) {
  int n = x.size();
  
  // Calculate means
  double mean_x = 0;
  double mean_y = 0;
  for (int i = 0; i < n; ++i) {
    mean_x += x[i];
    mean_y += y[i];
  }
  mean_x /= n;
  mean_y /= n;
  
  // Calculate coefficients
  double b1 = 0;
  double b0 = 0;
  double numerator = 0;
  double denominator = 0;
  for (int i = 0; i < n; ++i) {
    numerator += (x[i] - mean_x) * (y[i] - mean_y);
    denominator += (x[i] - mean_x) * (x[i] - mean_x);
  }
  b1 = numerator / denominator;
  b0 = mean_y - b1 * mean_x;
  
  // Calculate fitted values and residuals
  NumericVector fitted_values(n);
  for (int i = 0; i < n; ++i) {
    fitted_values[i] = b0 + b1 * x[i];
  }
  
  
  // Return coefficients, fitted values, residuals, and R-squared
  return List::create(
    _["coef"] = NumericVector::create(b0, b1),
    _["fitted.values"] = fitted_values
  );
}

// Cholesky decomposition of a symmetric positive-definite matrix A
// Returns lower triangular matrix L such that A = L * L^T
NumericMatrix cholesky_decomposition(NumericMatrix A) {
  int n = A.nrow();
  NumericMatrix L(n, n); // Initialize L with zeros
  
  for (int i = 0; i < n; ++i) {
    // Compute L(i, i)
    double sum = A(i, i);
    for (int k = 0; k < i; ++k) {
      sum -= L(i, k) * L(i, k);
    }
    L(i, i) = sqrt(sum > 0 ? sum : 0); // Ensure non-negative for stability
    
    // Compute L(j, i) for j > i
    for (int j = i + 1; j < n; ++j) {
      sum = A(j, i);
      for (int k = 0; k < i; ++k) {
        sum -= L(j, k) * L(i, k);
      }
      L(j, i) = sum / L(i, i);
    }
  }
  return L;
}

// Solve L * z = b (forward substitution, L is lower triangular)
NumericVector forward_substitution(NumericMatrix L, NumericVector b) {
  int n = L.nrow();
  NumericVector z(n);
  
  for (int i = 0; i < n; ++i) {
    double sum = b[i];
    for (int j = 0; j < i; ++j) {
      sum -= L(i, j) * z[j];
    }
    z[i] = sum / L(i, i);
  }
  return z;
}

// Solve L^T * x = z (back substitution, L^T is upper triangular)
NumericVector back_substitution(NumericMatrix L, NumericVector z) {
  int n = L.nrow();
  NumericVector x(n);
  
  for (int i = n - 1; i >= 0; --i) {
    double sum = z[i];
    for (int j = i + 1; j < n; ++j) {
      sum -= L(j, i) * x[j]; // L(j, i) is L^T(i, j)
    }
    x[i] = sum / L(i, i);
  }
  return x;
}

// [[Rcpp::export]]
List fast_lm_mult(NumericMatrix x, NumericVector y) {
  int n = x.nrow(); // Number of observations
  int p = x.ncol(); // Number of predictors
  
  // Add intercept term to the design matrix
  NumericMatrix X(n, p + 1);
  for (int i = 0; i < n; ++i) {
    X(i, 0) = 1; // Intercept column
    for (int j = 0; j < p; ++j) {
      X(i, j + 1) = x(i, j);
    }
  }
  
  // Compute X'X and X'y
  NumericMatrix XtX(p + 1, p + 1);
  NumericVector Xty(p + 1);
  for (int i = 0; i < p + 1; ++i) {
    for (int j = 0; j < p + 1; ++j) {
      double sum = 0;
      for (int k = 0; k < n; ++k) {
        sum += X(k, i) * X(k, j);
      }
      XtX(i, j) = sum;
    }
    double sum = 0;
    for (int k = 0; k < n; ++k) {
      sum += X(k, i) * y[k];
    }
    Xty[i] = sum;
  }
  
  // Cholesky decomposition: X'X = L * L^T
  NumericMatrix L = cholesky_decomposition(XtX);
  
  // Solve L * z = X'y (forward substitution)
  NumericVector z = forward_substitution(L, Xty);
  
  // Solve L^T * beta = z (back substitution)
  NumericVector coef = back_substitution(L, z);
  
  // Compute fitted values
  NumericVector fitted_values(n);
  for (int i = 0; i < n; ++i) {
    double sum = 0;
    for (int j = 0; j < p + 1; ++j) {
      sum += coef[j] * X(i, j);
    }
    fitted_values[i] = sum;
  }
  
  // Compute residuals
  NumericVector residuals = fitted_values - y;
  
  // Compute R-squared
  double y_mean = mean(y);
  double TSS = sum(pow(y - y_mean, 2));
  double RSS = sum(pow(residuals, 2));
  double R2 = 1 - RSS / TSS;
  
  return List::create(
    _["coefficients"] = coef,
    _["fitted.values"] = fitted_values,
    _["residuals"] = residuals,
    _["r.squared"] = R2
  );
}