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
  
  // Solve the normal equations (X'X)beta = X'y
  NumericVector coef(p + 1);
  for (int i = 0; i < p + 1; ++i) {
    double sum = Xty[i];
    for (int j = 0; j < p + 1; ++j) {
      sum -= XtX(i, j) * coef[j];
    }
    coef[i] = sum / XtX(i, i);
  }
  
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
  NumericVector residuals = y - fitted_values;
  
  // Compute R-squared
  double TSS = sum(pow(y - mean(y), 2)); // Total sum of squares
  double RSS = sum(pow(residuals, 2)); // Residual sum of squares
  double R2 = 1 - RSS / TSS;
  
  // Return coefficients, fitted values, residuals, and R-squared value
  return List::create(
    _["coefficients"] = coef,
    _["fitted.values"] = fitted_values,
    _["residuals"] = residuals,
    _["r.squared"] = R2
  );
}
