#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List stoch_superiority_cpp(NumericVector x, NumericVector y) {
  NumericVector xs = clone(x).sort();
  NumericVector ys = clone(y).sort();
  
  const int n_x = xs.size();
  const int n_y = ys.size();
  
  if (n_x == 0 || n_y == 0) {
    stop("x and y must both have positive length.");
  }
  
  long double less_count = 0.0L;
  long double tie_count  = 0.0L;
  
  int left  = 0;   // number of y strictly less than x[i]
  int right = 0;   // number of y less than or equal to x[i]
  
  for (int i = 0; i < n_x; ++i) {
    const double xi = xs[i];
    
    while (left < n_y && ys[left] < xi) {
      ++left;
    }
    while (right < n_y && ys[right] <= xi) {
      ++right;
    }
    
    less_count += left;
    tie_count  += (right - left);
  }
  
  const long double denom = static_cast<long double>(n_x) *
    static_cast<long double>(n_y);
  
  const double p_gt   = static_cast<double>(less_count / denom);
  const double p_tie  = static_cast<double>(tie_count  / denom);
  const double p_star = p_gt + 0.5 * p_tie;
  
  return List::create(
    Named("p_gt")   = p_gt,
    Named("p_tie")  = p_tie,
    Named("p_star") = p_star
  );
}
