// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;


double LPM_C(const double &degree, const double &target, const RVector<double> &variable) {
  size_t n = variable.size();
  double out=0;
  if (degree==0){
    for (size_t i = 0; i < n; i++)
      if (variable[i] <= target)
        out += 1;
  }else if(degree==1){
    for (size_t i = 0; i < n; i++) {
      if (variable[i] <= target)
        out += target-variable[i];
    }
  }else{
    for (size_t i = 0; i < n; i++) {
      if (variable[i] <= target)
        out += std::pow(target-variable[i], degree);
    }
  }
  out /= n;
  return(out);
}

double UPM_C(const double &degree, const double &target, const RVector<double> &variable) {
  size_t n = variable.size();
  double out=0;
  if (degree==0){
    for (size_t i = 0; i < n; i++)
      if (variable[i] > target)
        out += 1;
  }else if(degree==1){
    for (size_t i = 0; i < n; i++) {
      if (variable[i] > target)
        out += variable[i]-target;
    }
  }else{
    for (size_t i = 0; i < n; i++) {
      if (variable[i] > target)
        out += std::pow(variable[i]-target, degree);
    }
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
  LPM_Worker(const double degree, const NumericVector &target, const NumericVector &variable, NumericVector &output): 
    degree(degree), target(target), variable(variable), output(output) {}
  void operator()(std::size_t begin, std::size_t end) { for (size_t i = begin; i < end; i++) output[i] = LPM_C(degree, target[i], variable); }
};

struct UPM_Worker : public Worker
{
  const double degree;
  const RVector<double> target;
  const RVector<double> variable;
  RVector<double> output;
  UPM_Worker(const double degree, const NumericVector &target, const NumericVector &variable, NumericVector &output): 
    degree(degree), target(target), variable(variable), output(output) {}
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
NumericVector LPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  size_t target_size=target.size();
  NumericVector output = NumericVector(target_size);
  LPM_Worker tmp_func(degree, target, variable, output);
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
NumericVector UPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  size_t target_size=target.size();
  NumericVector output = NumericVector(target_size);
  UPM_Worker tmp_func(degree, target, variable,output);
  parallelFor(0, target_size, tmp_func);
  return(output);
}




////////////////////////////////////////////////////////////////////////////////
double CoUPM_C(
    const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  
  double out=0;
  bool dont_use_pow=(degree_upm==1 || degree_upm==0),
    d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    if (x1<=0) continue;
    
    double y1=(y[i]-target_y);
    if (y1<=0) continue;
    
    if(d_upm_0){
        x1=(x1==0?0:1);
        y1=(y1==0?0:1);
    }
    
    if(dont_use_pow)
      out += x1 * y1;
    else
      out += std::pow(x1, degree_upm) * std::pow(y1, degree_upm);
  }
  return out/max_size;
}

double CoLPM_C(
    const double &degree_lpm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  bool dont_use_pow=(degree_lpm==1 || degree_lpm==0),
    d_lpm_0=(degree_lpm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    if (x1<=0) continue;
    
    double y1=(target_y-y[i]);
    if (y1<=0) continue;
    
    if(d_lpm_0){
        x1=(x1==0?0:1);
        y1=(y1==0?0:1);
    }
    
    if(dont_use_pow)
      out += x1 * y1;
    else
      out += std::pow(x1, degree_lpm) * std::pow(y1, degree_lpm);
  }
  return out/max_size;
}


double DLPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  bool dont_use_pow_lpm=(degree_lpm==1 || degree_lpm==0), dont_use_pow_upm=(degree_upm==1 || degree_upm==0),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(x[i]-target_x);
    if (x1<=0) continue;
    
    double y1=(target_y-y[i]);
    if (y1<=0) continue;
    
    if(d_lpm_0) x1=(x1==0?0:1);
    if(d_upm_0) y1=(y1==0?0:1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm)
      out += x1 * y1;
    else if(dont_use_pow_lpm)
      out += x1 * std::pow(y1, degree_upm);
    else if(dont_use_pow_upm)
      out += std::pow(x1, degree_lpm) * y1;
    else
      out += std::pow(x1, degree_lpm) * std::pow(y1, degree_upm);
  }
  return out/max_size;
}


double DUPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y
){
  size_t n_x=x.size(), n_y=y.size();
  size_t max_size=(n_x>n_y?n_x:n_y);
  size_t min_size=(n_x<n_y?n_x:n_y);
  if (n_x!=n_y)
    Rcpp::warning("x vector length != y vector length");
  if (min_size<=0)   // if len = 0, return 0
    return 0;
  double out=0;
  bool dont_use_pow_lpm=(degree_lpm==1 || degree_lpm==0), dont_use_pow_upm=(degree_upm==1 || degree_upm==0),
    d_lpm_0=(degree_lpm==0), d_upm_0=(degree_upm==0);
  for(size_t i=0; i<min_size; i++){
    double x1=(target_x-x[i]);
    if (x1<=0) continue;
    
    double y1=(y[i]-target_y);
    if (y1<=0) continue;
    
    if(d_lpm_0) x1=(x1==0?0:1);
    if(d_upm_0) y1=(y1==0?0:1);
    
    if(dont_use_pow_lpm && dont_use_pow_upm)
      out += x1 * y1;
    else if(dont_use_pow_lpm)
      out += x1 * std::pow(y1, degree_upm);
    else if(dont_use_pow_upm)
      out += std::pow(x1, degree_lpm) * y1;
    else
      out += std::pow(x1, degree_lpm) * std::pow(y1, degree_upm);
  }
  return out/max_size;
}



struct CoLPM_Worker : public Worker
{
  const double degree_lpm;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_t_x;
  const size_t n_t_y;
  RVector<double> output;
  CoLPM_Worker(
    const double degree_lpm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y,
    NumericVector output
  ): 
    degree_lpm(degree_lpm), x(x), y(y), target_x(target_x), target_y(target_y),
    n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++) {
      output[i] = CoLPM_C(degree_lpm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); 
    }
  }
};

struct CoUPM_Worker : public Worker
{
  const double degree_upm;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_t_x;
  const size_t n_t_y;
  RVector<double> output;
  CoUPM_Worker(
    const double &degree_upm,
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y,
    NumericVector &output
  ): 
    degree_upm(degree_upm), x(x), y(y), target_x(target_x), target_y(target_y),
    n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++) {
      output[i] = CoUPM_C(degree_upm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); 
    }
  }
};

struct DLPM_Worker : public Worker
{
  const double degree_lpm;
  const double degree_upm;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_t_x;
  const size_t n_t_y;
  RVector<double> output;
  DLPM_Worker(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y,
    NumericVector &output
  ): 
    degree_lpm(degree_lpm), degree_upm(degree_upm), x(x), y(y), target_x(target_x), target_y(target_y), 
    n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++) {
      output[i] = DLPM_C(degree_lpm, degree_upm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); 
    }
  }
};

struct DUPM_Worker : public Worker
{
  const double degree_lpm;
  const double degree_upm;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<double> target_x;
  const RVector<double> target_y;
  const size_t n_t_x;
  const size_t n_t_y;
  RVector<double> output;
  DUPM_Worker(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y,
    NumericVector &output
  ): 
    degree_lpm(degree_lpm), degree_upm(degree_upm), x(x), y(y), target_x(target_x), target_y(target_y), 
    n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)
  {}
  void operator()(std::size_t begin, std::size_t end) { 
    for (size_t i = begin; i < end; i++){
      output[i] = DUPM_C(degree_lpm, degree_upm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); 
    }
  }
};



//' Co-Lower Partial Moment
//' (Lower Left Quadrant 4)
//'
//' This function generates a co-lower partial moment for between two equal length variables for any degree or target.
//' @param degree_lpm integer; Degree for variable X and Y.  \code{(degree_lpm = 0)} is frequency, \code{(degree_lpm = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Target for lower deviations of variable X.  Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Target for lower deviations of variable Y.  Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
//' @return Co-LPM of two variables
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' Co.LPM(0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("Co.LPM", rng = false)]]
NumericVector CoLPM_CPv(
    const double degree_lpm, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  CoLPM_Worker tmp_func(degree_lpm, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}


//' Co-Upper Partial Moment
//' (Upper Right Quadrant 1)
//'
//' This function generates a co-upper partial moment between two equal length variables for any degree or target.
//' @param degree_upm integer; Degree for variable X and Y.  \code{(degree_upm = 0)} is frequency, \code{(degree_upm = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Target for upside deviations of variable X.  Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Target for upside deviations of variable Y.  Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
//' @return Co-UPM of two variables
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100)
//' Co.UPM(0, x, y, mean(x), mean(y))
//' @export
// [[Rcpp::export("Co.UPM", rng = false)]]
NumericVector CoUPM_CPv(
    const double degree_upm, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y 
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  CoUPM_Worker tmp_func(degree_upm, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}


//' Divergent-Lower Partial Moment
//' (Lower Right Quadrant 3)
//'
//' This function generates a divergent lower partial moment between two equal length variables for any degree or target.
//' @param degree_lpm integer; Degree for variable X.  \code{(degree_lpm = 0)} is frequency, \code{(degree_lpm = 1)} is area.
//' @param degree_upm integer; Degree for variable Y.  \code{(degree_upm = 0)} is frequency, \code{(degree_upm = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Target for upside deviations of variable X.  Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Target for lower deviations of variable Y.  Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
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
    const double degree_lpm, const double degree_upm, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  DLPM_Worker tmp_func(degree_lpm, degree_upm, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}


//' Divergent-Upper Partial Moment
//' (Upper Left Quadrant 2)
//'
//' This function generates a divergent upper partial moment between two equal length variables for any degree or target.
//' @param degree_lpm integer; Degree for variable X.  \code{(degree_lpm = 0)} is frequency, \code{(degree_lpm = 1)} is area.
//' @param degree_upm integer; Degree for variable Y.  \code{(degree_upm = 0)} is frequency, \code{(degree_upm = 1)} is area.
//' @param x a numeric vector.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param y a numeric vector of equal length to \code{x}.   \link{data.frame} or \link{list} type objects are not permissible.
//' @param target_x numeric; Target for lower deviations of variable X.  Typically the mean of Variable X for classical statistics equivalences, but does not have to be.
//' @param target_y numeric; Target for upper deviations of variable Y.  Typically the mean of Variable Y for classical statistics equivalences, but does not have to be.
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
    const double degree_lpm, const double degree_upm, 
    const NumericVector x, const NumericVector y, 
    const NumericVector target_x, const NumericVector target_y
) {
  size_t target_x_size=target_x.size();
  size_t target_y_size=target_y.size();
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size);
  NumericVector output = NumericVector(max_target_size);
  DUPM_Worker tmp_func(degree_lpm, degree_upm, x, y, target_x, target_y, output);
  parallelFor(0, output.size(), tmp_func);
  return(output);
}






void PMMatrix_CPv(
  const double &degree_lpm, 
  const double &degree_upm, 
  const RMatrix<double>::Column &x, 
  const RMatrix<double>::Column &y, 
  const double &target_x,
  const double &target_y, 
  const bool &pop_adj, 
  const size_t &rows, 
  double &coLpm,
  double &coUpm,   
  double &dLpm, 
  double &dUpm
){
  coLpm=0;
  coUpm=0;
  dLpm=0;
  dUpm=0; 
  if(rows == 0)
    return;
  bool dont_use_pow_lpm=(degree_lpm==1 || degree_lpm==0),
       dont_use_pow_upm=(degree_upm==1 || degree_upm==0),
       d_lpm_0=(degree_lpm==0),
       d_upm_0=(degree_upm==0);
  bool dont_use_pow = (dont_use_pow_lpm && dont_use_pow_upm);
  //double x_less_target_upm, target_less_x_upm, y_less_target_upm, target_less_y_upm;
  double x_less_target_upm, y_less_target_upm, target_less_y_upm;
  for(size_t i = 0; i < rows; i++){
    double x_less_target_lpm = (target_x - x[i]),
           target_less_x_lpm = (x[i] - target_x),
           y_less_target_lpm = (target_y - y[i]),
           target_less_y_lpm = (y[i] - target_y);
    bool   less_equal_zero_target_less_x = (target_less_x_lpm <= 0),
           less_equal_zero_target_less_y = (target_less_y_lpm <= 0),
           less_equal_zero_x_less_target = (x_less_target_lpm <= 0),
           less_equal_zero_y_less_target = (y_less_target_lpm <= 0);
    if(
		less_equal_zero_target_less_x && 
		less_equal_zero_target_less_y && 
		less_equal_zero_x_less_target && 
		less_equal_zero_y_less_target
	)	// nothing to do with this row
        continue;
    // UPM values
	if(d_upm_0){ // pow(0,0) == 0, pow(x,0) == 1
        x_less_target_upm=(x_less_target_lpm==0?0:1);
        //target_less_x_upm=(target_less_x_lpm==0?0:1);
        y_less_target_upm=(y_less_target_lpm==0?0:1);
        target_less_y_upm=(target_less_y_lpm==0?0:1);
    }else{
        x_less_target_upm=x_less_target_lpm;
        //target_less_x_upm=target_less_x_lpm;
        y_less_target_upm=y_less_target_lpm;
        target_less_y_upm=target_less_y_lpm;
    }
	// LPM values
    if(d_lpm_0){ // pow(0,0) == 0, pow(x,0) == 1
        x_less_target_lpm=(x_less_target_lpm==0?0:1);
        target_less_x_lpm=(target_less_x_lpm==0?0:1);
        y_less_target_lpm=(y_less_target_lpm==0?0:1);
        target_less_y_lpm=(target_less_y_lpm==0?0:1);
    }

	/*
	CoLPM_C
		double x1=(target_x-x[i]);
		double y1=(target_y-y[i]);
		out += std::pow(x1, degree_lpm) * std::pow(y1, degree_lpm);
	*/
    if(!less_equal_zero_target_less_x && !less_equal_zero_target_less_y){
        if(dont_use_pow_lpm)
            coLpm += target_less_x_lpm * target_less_y_lpm;
        else 
            coLpm += std::pow(target_less_x_lpm, degree_lpm) * std::pow(target_less_y_lpm, degree_lpm);
    }

	/*
	CoUPM_C
		double x1=(x[i]-target_x);
		double y1=(y[i]-target_y);
		out += std::pow(x1, degree_upm) * std::pow(y1, degree_upm);
	*/
    if(!less_equal_zero_x_less_target && !less_equal_zero_y_less_target){
        if(dont_use_pow_upm)
            coUpm += x_less_target_upm * y_less_target_upm;
        else 
            coUpm += std::pow(x_less_target_upm, degree_upm) * std::pow(y_less_target_upm, degree_upm);
    }

	/*
	DLPM_C(
		double x1=(x[i]-target_x);
		double y1=(target_y-y[i]);
		out += std::pow(x1, degree_lpm) * std::pow(y1, degree_upm);
	*/
    if(!less_equal_zero_x_less_target && !less_equal_zero_target_less_y){
        if(dont_use_pow)
            dLpm += x_less_target_lpm * target_less_y_upm;
        else if(dont_use_pow_lpm)
            dLpm += x_less_target_lpm * std::pow(target_less_y_upm, degree_upm);
        else if(dont_use_pow_upm)
            dLpm += std::pow(x_less_target_lpm, degree_lpm) * target_less_y_upm;
        else 
            dLpm += std::pow(x_less_target_lpm, degree_lpm) * std::pow(target_less_y_upm, degree_upm);
    }

	/*
	DUPM_C(
		double x1=(target_x-x[i]);
		double y1=(y[i]-target_y);
		out += std::pow(x1, degree_lpm) * std::pow(y1, degree_upm);

	*/
    if(!less_equal_zero_target_less_x && !less_equal_zero_y_less_target){
        if(dont_use_pow)
            dUpm += target_less_x_lpm * y_less_target_upm;
        else if(dont_use_pow_lpm)
            dUpm += target_less_x_lpm * std::pow(y_less_target_upm, degree_upm);
        else if(dont_use_pow_upm)
            dUpm += std::pow(target_less_x_lpm, degree_lpm) * y_less_target_upm;
        else 
            dUpm += std::pow(target_less_x_lpm, degree_lpm) * std::pow(y_less_target_upm, degree_upm);
    }
  }
  coLpm/=rows;
  coUpm/=rows;
  dLpm/=rows;
  dUpm/=rows;
  if(pop_adj && rows>1){
    double adjust=rows/(rows-1);
    coLpm*=adjust;
    coUpm*=adjust;
    dLpm*=adjust;
    dUpm*=adjust;
  }
}

struct PMMatrix_Worker : public Worker
{
  const double degree_lpm;
  const double degree_upm;
  const RMatrix<double> variable;
  const RVector<double> target;
  const size_t variable_cols;
  const size_t variable_rows;
  const size_t target_length;
  const bool pop_adj;
  RMatrix<double> coLpm;
  RMatrix<double> coUpm;
  RMatrix<double> dLpm;
  RMatrix<double> dUpm;
  PMMatrix_Worker(
    const double &degree_lpm, const double &degree_upm, 
    const NumericMatrix &variable, 
    const NumericVector &target,
    const bool &pop_adj,
    NumericMatrix &coLpm, NumericMatrix &coUpm,
    NumericMatrix &dLpm,  NumericMatrix &dUpm
  ): 
    degree_lpm(degree_lpm), degree_upm(degree_upm), 
    variable(variable), target(target),
    variable_cols(variable.cols()), variable_rows(variable.rows()), target_length(target.size()), 
	pop_adj(pop_adj),
    coLpm(coLpm), coUpm(coUpm), 
    dLpm(dLpm), dUpm(dUpm)
  {
    if(variable_cols != target_length)
      Rcpp::stop("varible matrix cols != target vector length");
  }
  void operator()(std::size_t begin, std::size_t end) {
    for (size_t i = begin; i < end; i++){
      for (size_t l = 0; l < variable_cols; l++){
        PMMatrix_CPv(
          degree_lpm,
          degree_upm,
          variable.column(i),
          variable.column(l),
          target[i],
          target[l],
		  pop_adj,
          variable_rows,
          coLpm(i,l),
          coUpm(i,l),
          dLpm(i,l),
          dUpm(i,l)
        );
      }
    }
  }
};


//' Partial Moment Matrix - Internal Use, without covariance and pop_adjust
//'
//'
//' This function generates a co-partial moment matrix for the specified co-partial moment.
//' @param LPM_degree integer; Degree for \code{variable} below \code{target} deviations.  \code{(LPM_degree = 0)} is frequency, \code{(LPM_degree = 1)} is area.
//' @param UPM_degree integer; Degree for \code{variable} above \code{target} deviations.  \code{(UPM_degree = 0)} is frequency, \code{(UPM_degree = 1)} is area.
//' @param target numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)  \code{(target = NULL)} (default) will set the target as the mean of every variable.
//' @param variable a numeric matrix or data.frame.
//' @param pop_adj logical; \code{FALSE} (default) Adjusts the sample co-partial moment matrices for population statistics.
//' @return Matrix of partial moment quadrant values (CUPM, DUPM, DLPM, CLPM)
//' @note For divergent asymmetical \code{"D.LPM" and "D.UPM"} matrices, matrix is \code{D.LPM(column,row,...)}.
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @references Viole, F. (2017) "Bayes' Theorem From Partial Moments"
//' \url{https://www.ssrn.com/abstract=3457377}
//' @examples
//' set.seed(123)
//' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
//' A <- cbind(x,y,z)
//' PMMatrix_CPv(LPM_degree = 1, UPM_degree = 1, variable = A)
//'
//' ## Use of vectorized numeric targets (target_x, target_y, target_z)
//' PMMatrix_CPv(LPM_degree = 1, UPM_degree = 1, target = c(0, 0.15, .25), variable = A)
//'
//' ## Calling Individual Partial Moment Quadrants
//' cov.mtx <- PMMatrix_CPv(LPM_degree = 1, UPM_degree = 1, variable = A)
//' cov.mtx$cupm
//'
//' @export
// [[Rcpp::export("PMMatrix_CPv", rng = false)]]
List PMMatrix_CPv(
    const double LPM_degree, 
    const double UPM_degree, 
    const NumericVector target, 
    const NumericMatrix variable,
	const bool pop_adj
) {
  size_t variable_cols=variable.cols();
  size_t target_length=target.size();
  if(variable_cols != target_length){
    Rcpp::stop("varible matrix cols != target vector length");
    return List::create();
  }
  NumericMatrix coLpm(variable_cols, variable_cols);
  NumericMatrix coUpm(variable_cols, variable_cols);
  NumericMatrix dLpm(variable_cols, variable_cols);
  NumericMatrix dUpm(variable_cols, variable_cols);
  PMMatrix_Worker tmp_func(LPM_degree, UPM_degree, variable, target, pop_adj, coLpm, coUpm, dLpm, dUpm);
  parallelFor(0, variable_cols, tmp_func);
  return(
    List::create(
      Named("cupm") = coUpm,
      Named("dupm") = dUpm,
      Named("dlpm") = dLpm,
      Named("clpm") = coLpm
    )
  );
}
