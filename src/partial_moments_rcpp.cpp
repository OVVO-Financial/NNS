// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;


double LPM_C(const double degree, const double target, const RVector<double> variable) {
  size_t n = variable.size();
  double out=0;
  for (size_t i = 0; i < n; i++) {
    if (variable[i] <= target) 
      out += std::pow(target-variable[i], degree);
  }
  out /= n;
  return(out);
}

double UPM_C(const double degree, const double target, const RVector<double> variable) {
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
  const RVector<double> target;
  const RVector<double> variable;
  RVector<double> output;
  LPM_Worker(const double degree, const NumericVector target, const NumericVector variable, NumericVector output): degree(degree), target(target), variable(variable), output(output) {}
  void operator()(std::size_t begin, std::size_t end) { for (size_t i = begin; i < end; i++) output[i] = LPM_C(degree, target[i], variable); }
};

struct UPM_Worker : public Worker
{
  const double degree;
  const RVector<double> target;
  const RVector<double> variable;
  RVector<double> output;
  UPM_Worker(const double degree, const NumericVector target, const NumericVector variable, NumericVector output) : degree(degree), target(target), variable(variable), output(output) {}
  void operator()(std::size_t begin, std::size_t end) { for (size_t i = begin; i < end; i++) output[i] = UPM_C(degree, target[i], variable); }
};


//' Lower Partial Moment
//'
//' This function generates a univariate lower partial moment for any degree or target.
//'
//' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
//' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
//' @param variable a numeric vector.  \link{data.frame} or \link{list} type objects are not permissible.
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
//' @param variable a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
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
double CoUPM_C(
    const double degree_x, const double degree_y, 
    const RVector<double> x, const RVector<double> y, 
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

double CoLPM_C(
    const double degree_x, const double degree_y, 
    const RVector<double> x, const RVector<double> y, 
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


double DLPM_C(
    const double degree_x, const double degree_y, 
    const RVector<double> x, const RVector<double> y, 
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


double DUPM_C(
    const double degree_x, const double degree_y, 
    const RVector<double> x, const RVector<double> y, 
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
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_tx;
  const size_t n_ty;
  RVector<double> output;
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
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_tx;
  const size_t n_ty;
  RVector<double> output;
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
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_tx;
  const size_t n_ty;
  RVector<double> output;
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
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_tx;
  const size_t n_ty;
  RVector<double> output;
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


//' Co-Lower Partial Moment
//' (Lower Left Quadrant 4)
//'
//' This function generates a co-lower partial moment for between two equal length variables for any degree or target.
//' @param degree_x integer; Degree for variable X.  \code{(degree_x = 0)} is frequency, \code{(degree_x = 1)} is area.
//' @param degree_y integer; Degree for variable Y.  \code{(degree_y = 0)} is frequency, \code{(degree_y = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
//' @return Co-LPM of two variables
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' Co.LPM(0, 0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("Co.LPM", rng = false)]]
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


//' Co-Upper Partial Moment
//' (Upper Right Quadrant 1)
//'
//' This function generates a co-upper partial moment between two equal length variables for any degree or target.
//' @param degree_x integer; Degree for variable X.  \code{(degree_x = 0)} is frequency, \code{(degree_x = 1)} is area.
//' @param degree_y integer; Degree for variable Y.  \code{(degree_y = 0)} is frequency, \code{(degree_y = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
//' @return Co-UPM of two variables
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' Co.UPM(0, 0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("Co.UPM", rng = false)]]
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


//' Divergent-Lower Partial Moment
//' (Lower Right Quadrant 3)
//'
//' This function generates a divergent lower partial moment between two equal length variables for any degree or target.
//' @param degree_x integer; Degree for variable X.  \code{(degree_x = 0)} is frequency, \code{(degree_x = 1)} is area.
//' @param degree_y integer; Degree for variable Y.  \code{(degree_y = 0)} is frequency, \code{(degree_y = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
//' @return Divergent LPM of two variables
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' D.LPM(0, 0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("D.LPM", rng = false)]]
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


//' Divergent-Upper Partial Moment
//' (Upper Left Quadrant 2)
//'
//' This function generates a divergent upper partial moment between two equal length variables for any degree or target.
//' @param degree_x integer; Degree for variable X.  \code{(degree_x = 0)} is frequency, \code{(degree_x = 1)} is area.
//' @param degree_y integer; Degree for variable Y.  \code{(degree_y = 0)} is frequency, \code{(degree_y = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
//' @return Divergent UPM of two variables
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' D.UPM(0, 0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("D.UPM", rng = false)]]
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
