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



