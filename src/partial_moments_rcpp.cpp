// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::export(rng = false)]]
double LPM_C(const double degree, const double target, const NumericVector variable) {
  size_t n = variable.size();
  double out=0;
  for (size_t i = 0; i < n; i++) {
    if (variable[i] <= target) 
      out += std::pow(target-variable[i], degree);
  }
  out /= n;
  return(out);
}
// [[Rcpp::export(rng = false)]]
double UPM_C(const double degree, const double target, const NumericVector variable) {
  size_t n = variable.size();
  double out=0;
  for (size_t i = 0; i < n; i++) {
    if (variable[i] > target) 
      out += std::pow(variable[i]-target, degree);
  }
  out /= n;
  return(out);
}

struct LPM_Worker : public Worker
{
  const double degree;
  const NumericVector target;
  const NumericVector variable;
  NumericVector output;
  LPM_Worker(const double degree, const NumericVector target, const NumericVector variable, NumericVector output): degree(degree), target(target), variable(variable), output(output) {}
  void operator()(std::size_t begin, std::size_t end) { for (size_t i = begin; i < end; i++) output[i] = LPM_C(degree, target[i], variable); }
};

struct UPM_Worker : public Worker
{
  const double degree;
  const NumericVector target;
  const NumericVector variable;
  NumericVector output;
  UPM_Worker(const double degree, const NumericVector target, const NumericVector variable, NumericVector output) : degree(degree), target(target), variable(variable), output(output) {}
  void operator()(std::size_t begin, std::size_t end) { for (size_t i = begin; i < end; i++) output[i] = UPM_C(degree, target[i], variable); }
};


//' Lower Partial Moment
//'
//' This function generates a univariate lower partial moment for any degree or target.
//'
//' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
//' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
//' @param variable a numeric vector.
//' @return LPM of variable
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100)
//' LPM(0, mean(x), x)
//' @export
// [[Rcpp::export("LPM", rng = false)]]
NumericVector LPM_CPv(const double degree, const NumericVector target, const NumericVector variable) {
  size_t target_size=target.size();
  NumericVector output = NumericVector(target_size);
  LPM_Worker tmp_func(degree,target,variable,output);
  parallelFor(0, target_size, tmp_func);
  return(output);
}


//' Upper Partial Moment
//'
//' This function generates a univariate upper partial moment for any degree or target.
//' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
//' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
//' @param variable a numeric vector.
//' @return UPM of variable
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100)
//' UPM(0, mean(x), x)
//' @export
// [[Rcpp::export("UPM", rng = false)]]
NumericVector UPM_CPv(const double degree, const NumericVector target, const NumericVector variable) {
  size_t target_size=target.size();
  NumericVector output = NumericVector(target_size);
  UPM_Worker tmp_func(degree,target,variable,output);
  parallelFor(0, target_size, tmp_func);
  return(output);
}




////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export(rng = false)]]
double CoUPM_C(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const double target_x, const double target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i] - target_x);
    if (x1<=0) continue;
    
    double y1=(y[i] - target_y);
    if (y1<=0) continue;
    
    out += std::pow(x1, degree_x) * std::pow(y1, degree_y);
  }
  return out/max_size;
}

// [[Rcpp::export(rng = false)]]
double CoLPM_C(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const double target_x, const double target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    if (x1<=0) continue;
    
    double y1=(target_y-y[i]);
    if (y1<=0) continue;
    
    out += std::pow(x1, degree_x) * std::pow(y1, degree_y);
  }
  return out/max_size;
}


// [[Rcpp::export(rng = false)]]
double DLPM_C(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const double target_x, const double target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    if (x1<=0) continue;
    
    double y1=(target_y-y[i]);
    if (y1<=0) continue;
    
    out += std::pow(x1, degree_x) * std::pow(y1, degree_y);
  }
  return out/max_size;
}

// [[Rcpp::export(rng = false)]]
double DUPM_C(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const double target_x, const double target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    if (x1<=0) continue;
    
    double y1=(y[i]-target_y);
    if (y1<=0) continue;
    
    out += std::pow(x1, degree_x) * std::pow(y1, degree_y);
  }
  return out/max_size;
}





struct CoLPM_Worker : public Worker
{
  const double degree_x;
  const double degree_y;
  const NumericVector x;
  const NumericVector y;
  const NumericVector target_x;
  const NumericVector target_y;
  const size_t n_tx;
  const size_t n_ty;
  NumericVector output;
  CoLPM_Worker(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y,
    NumericVector output
  ): degree_x(degree_x), degree_y(degree_y), x(x), y(y), target_x(target_x), target_y(target_y), 
  n_tx(target_x.size()), n_ty(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++) {
      output[i] = CoLPM_C(degree_x, degree_y, x, y, target_x[i%n_tx], target_y[i%n_ty]); 
    }
  }
};

struct CoUPM_Worker : public Worker
{
  const double degree_x;
  const double degree_y;
  const NumericVector x;
  const NumericVector y;
  const NumericVector target_x;
  const NumericVector target_y;
  const size_t n_tx;
  const size_t n_ty;
  NumericVector output;
  CoUPM_Worker(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y,
    NumericVector output
  ): degree_x(degree_x), degree_y(degree_y), x(x), y(y), target_x(target_x), target_y(target_y), 
  n_tx(target_x.size()), n_ty(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++) {
      output[i] = CoUPM_C(degree_x, degree_y, x, y, target_x[i%n_tx], target_y[i%n_ty]); 
    }
  }
};
struct DLPM_Worker : public Worker
{
  const double degree_x;
  const double degree_y;
  const NumericVector x;
  const NumericVector y;
  const NumericVector target_x;
  const NumericVector target_y;
  const size_t n_tx;
  const size_t n_ty;
  NumericVector output;
  DLPM_Worker(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y,
    NumericVector output
  ): degree_x(degree_x), degree_y(degree_y), x(x), y(y), target_x(target_x), target_y(target_y), 
  n_tx(target_x.size()), n_ty(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++) {
      output[i] = DLPM_C(degree_x, degree_y, x, y, target_x[i%n_tx], target_y[i%n_ty]); 
    }
  }
};
struct DUPM_Worker : public Worker
{
  const double degree_x;
  const double degree_y;
  const NumericVector x;
  const NumericVector y;
  const NumericVector target_x;
  const NumericVector target_y;
  const size_t n_tx;
  const size_t n_ty;
  NumericVector output;
  DUPM_Worker(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y,
    NumericVector output
  ): degree_x(degree_x), degree_y(degree_y), x(x), y(y), target_x(target_x), target_y(target_y), 
  n_tx(target_x.size()), n_ty(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++){
      output[i] = DUPM_C(degree_x, degree_y, x, y, target_x[i%n_tx], target_y[i%n_ty]); 
	}
  }
};


// [[Rcpp::export("Co.LPM_rcpp", rng = false)]]
NumericVector CoLPM_CPv(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  CoLPM_Worker tmp_func(degree_x, degree_y, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}


// [[Rcpp::export("Co.UPM_rcpp", rng = false)]]
NumericVector CoUPM_CPv(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  CoUPM_Worker tmp_func(degree_x, degree_y, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}


// [[Rcpp::export("D.LPM_rcpp", rng = false)]]
NumericVector DLPM_CPv(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  DLPM_Worker tmp_func(degree_x, degree_y, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}


// [[Rcpp::export("D.UPM_rcpp", rng = false)]]
NumericVector DUPM_CPv(
    const double degree_x, const double degree_y, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  DUPM_Worker tmp_func(degree_x, degree_y, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}
