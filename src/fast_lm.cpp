#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


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
