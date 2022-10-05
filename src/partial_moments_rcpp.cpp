// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "partial_moments_rcpp.h"
#include "partial_moments.h"
using namespace Rcpp;

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
NumericVector LPM_RCPP(const double &degree, const RObject &target, const RObject &variable) {
  NumericVector target_vec, variable_vec;
  if (is<NumericVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else if (is<IntegerVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else
    variable_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(variable, "unlist"), "as.vector");
  if (is<NumericVector>(target) && !target.isNULL()){
	target_vec = as<NumericVector>(target);
  }else{
	target_vec = NumericVector(1);
	target_vec[0] = mean(variable_vec);
  }
  return LPM_CPv(degree, target_vec, variable_vec);
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
NumericVector UPM_RCPP(const double &degree, const RObject &target, const RObject &variable) {
  NumericVector target_vec, variable_vec;
  if (is<NumericVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else if (is<IntegerVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else
    variable_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(variable, "unlist"), "as.vector");
  if (is<NumericVector>(target) && !target.isNULL()){
	target_vec = as<NumericVector>(target);
  }else{
	target_vec = NumericVector(1);
	target_vec[0] = mean(variable_vec);
  }
  return UPM_CPv(degree, target_vec, variable_vec);
}


//' Lower Partial Moment RATIO
//'
//' This function generates a standardized univariate lower partial moment for any degree or target.
//' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
//' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
//' @param variable a numeric vector.
//' @return Standardized LPM of variable
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @references Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"
//' \url{https://www.ssrn.com/abstract=3007373}
//' @examples
//' set.seed(123)
//' x <- rnorm(100)
//' LPM.ratio(0, mean(x), x)
//'
//' \dontrun{
//' ## Empirical CDF (degree = 0)
//' lpm_cdf <- LPM.ratio(0, sort(x), x)
//' plot(sort(x), lpm_cdf)
//'
//' ## Continuous CDF (degree = 1)
//' lpm_cdf_1 <- LPM.ratio(1, sort(x), x)
//' plot(sort(x), lpm_cdf_1)
//'
//' ## Joint CDF
//' x <- rnorm(5000) ; y <- rnorm(5000)
//' plot3d(x, y, Co.LPM(0, sort(x), sort(y), x, y), col = "blue", xlab = "X", ylab = "Y",
//' zlab = "Probability", box = FALSE)
//' }
//' @export
// [[Rcpp::export("LPM.ratio", rng = false)]]
NumericVector LPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable) {
  NumericVector target_vec, variable_vec;
  if (is<NumericVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else if (is<IntegerVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else if (is<DataFrame>(variable))
    variable_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(variable, "unlist"), "as.vector");
  else
	Rcpp::stop("variable should be numeric vector, or data table");
  if (is<NumericVector>(target) && !target.isNULL()){
	target_vec = as<NumericVector>(target);
  }else{
	target_vec = NumericVector(1);
	target_vec[0] = mean(variable_vec);
  }
  return LPM_ratio_CPv(degree, target_vec, variable_vec);
}


//' Upper Partial Moment RATIO
//'
//' This function generates a standardized univariate upper partial moment for any degree or target.
//' @param degree integer; \code{(degree = 0)} is frequency, \code{(degree = 1)} is area.
//' @param target numeric; Typically set to mean, but does not have to be. (Vectorized)
//' @param variable a numeric vector.
//' @return Standardized UPM of variable
//' @author Fred Viole, OVVO Financial Systems
//' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
//' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
//' @examples
//' set.seed(123)
//' x <- rnorm(100)
//' UPM.ratio(0, mean(x), x)
//'
//' ## Joint Upper CDF
//' \dontrun{
//' x <- rnorm(5000) ; y <- rnorm(5000)
//' plot3d(x, y, Co.UPM(0, sort(x), sort(y), x, y), col = "blue", xlab = "X", ylab = "Y",
//' zlab = "Probability", box = FALSE)
//' }
//' @export
// [[Rcpp::export("UPM.ratio", rng = false)]]
NumericVector UPM_ratio_RCPP(const double &degree, const RObject &target, const RObject &variable) {
  NumericVector target_vec, variable_vec;
  if (is<NumericVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else if (is<IntegerVector>(variable))
    variable_vec=as<NumericVector>(variable);
  else if (is<DataFrame>(variable))
    variable_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(variable, "unlist"), "as.vector");
  else
	Rcpp::stop("variable should be numeric vector, or data table");
  if (is<NumericVector>(target) && !target.isNULL()){
	target_vec = as<NumericVector>(target);
  }else{
	target_vec = NumericVector(1);
	target_vec[0] = mean(variable_vec);
  }
  return UPM_ratio_CPv(degree, target_vec, variable_vec);
}


//' Co-Lower Partial Moment
//' (Lower Left Quadrant 4)
//'
//' This function generates a co-lower partial moment for between two equal length variables for any degree or target.
//' @param degree_lpm integer; Degree for lower deviations of both variable X and Y.  \code{(degree_lpm = 0)} is frequency, \code{(degree_lpm = 1)} is area.
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
NumericVector CoLPM_RCPP(
    const double &degree_lpm, 
    const RObject &x, const RObject &y, 
    const RObject &target_x, const RObject &target_y
) {
  NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
  if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
  else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
  else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
  else                         Rcpp::stop("x should be numeric vector, or data table");

  if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
  else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
  else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
  else                         Rcpp::stop("y should be numeric vector, or data table");

  if (is<NumericVector>(target_x) && !target_x.isNULL()){
	target_x_vec = as<NumericVector>(target_x);
  }else{
	target_x_vec = NumericVector(1);
	target_x_vec[0] = mean(x_vec);
  }
  if (is<NumericVector>(target_y) && !target_y.isNULL()){
	target_y_vec = as<NumericVector>(target_y);
  }else{
	target_y_vec = NumericVector(1);
	target_y_vec[0] = mean(y_vec);
  }
  return CoLPM_CPv(degree_lpm, x_vec, y_vec, target_x_vec, target_y_vec);
}


//' Co-Upper Partial Moment
//' (Upper Right Quadrant 1)
//'
//' This function generates a co-upper partial moment between two equal length variables for any degree or target.
//' @param degree_upm integer; Degree for upper variations of both variable X and Y.  \code{(degree_upm = 0)} is frequency, \code{(degree_upm = 1)} is area.
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
NumericVector CoUPM_RCPP(
    const double &degree_upm, 
    const RObject &x, const RObject &y, 
    const RObject &target_x, const RObject &target_y
) {
  NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
  if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
  else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
  else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
  else                         Rcpp::stop("x should be numeric vector, or data table");

  if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
  else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
  else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
  else                         Rcpp::stop("y should be numeric vector, or data table");

  if (is<NumericVector>(target_x) && !target_x.isNULL()){
	target_x_vec = as<NumericVector>(target_x);
  }else{
	target_x_vec = NumericVector(1);
	target_x_vec[0] = mean(x_vec);
  }
  if (is<NumericVector>(target_y) && !target_y.isNULL()){
	target_y_vec = as<NumericVector>(target_y);
  }else{
	target_y_vec = NumericVector(1);
	target_y_vec[0] = mean(y_vec);
  }
  return CoUPM_CPv(degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
}


//' Divergent-Lower Partial Moment
//' (Lower Right Quadrant 3)
//'
//' This function generates a divergent lower partial moment between two equal length variables for any degree or target.
//' @param degree_lpm integer; Degree for lower deviations of variable Y.  \code{(degree_lpm = 0)} is frequency, \code{(degree_lpm = 1)} is area.
//' @param degree_upm integer; Degree for upper deviations of variable X.  \code{(degree_upm = 0)} is frequency, \code{(degree_upm = 1)} is area.
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
NumericVector DLPM_RCPP(
    const double &degree_lpm, const double &degree_upm, 
    const RObject &x, const RObject &y, 
    const RObject &target_x, const RObject &target_y
) {
  NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
  if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
  else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
  else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
  else                         Rcpp::stop("x should be numeric vector, or data table");

  if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
  else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
  else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
  else                         Rcpp::stop("y should be numeric vector, or data table");

  if (is<NumericVector>(target_x) && !target_x.isNULL()){
	target_x_vec = as<NumericVector>(target_x);
  }else{
	target_x_vec = NumericVector(1);
	target_x_vec[0] = mean(x_vec);
  }
  if (is<NumericVector>(target_y) && !target_y.isNULL()){
	target_y_vec = as<NumericVector>(target_y);
  }else{
	target_y_vec = NumericVector(1);
	target_y_vec[0] = mean(y_vec);
  }
  return DLPM_CPv(degree_lpm, degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
}


//' Divergent-Upper Partial Moment
//' (Upper Left Quadrant 2)
//'
//' This function generates a divergent upper partial moment between two equal length variables for any degree or target.
//' @param degree_lpm integer; Degree for lower deviations of variable X.  \code{(degree_lpm = 0)} is frequency, \code{(degree_lpm = 1)} is area.
//' @param degree_upm integer; Degree for upper deviations of variable Y.  \code{(degree_upm = 0)} is frequency, \code{(degree_upm = 1)} is area.
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
NumericVector DUPM_RCPP(
    const double &degree_lpm, const double &degree_upm, 
    const RObject &x, const RObject &y, 
    const RObject &target_x, const RObject &target_y
) {
  NumericVector target_x_vec, target_y_vec, x_vec, y_vec;
  if (is<NumericVector>(x))    x_vec=as<NumericVector>(x);
  else if (is<IntegerVector>(x))	x_vec=as<NumericVector>(x);
  else if (is<DataFrame>(x))   x_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(x, "unlist"), "as.vector");
  else                         Rcpp::stop("x should be numeric vector, or data table");

  if (is<NumericVector>(y))    y_vec=as<NumericVector>(y);
  else if (is<IntegerVector>(y))	y_vec=as<NumericVector>(y);
  else if (is<DataFrame>(y))   y_vec=Rcpp::internal::convert_using_rfunction(Rcpp::internal::convert_using_rfunction(y, "unlist"), "as.vector");
  else                         Rcpp::stop("y should be numeric vector, or data table");

  if (is<NumericVector>(target_x) && !target_x.isNULL()){
	target_x_vec = as<NumericVector>(target_x);
  }else{
	target_x_vec = NumericVector(1);
	target_x_vec[0] = mean(x_vec);
  }
  if (is<NumericVector>(target_y) && !target_y.isNULL()){
	target_y_vec = as<NumericVector>(target_y);
  }else{
	target_y_vec = NumericVector(1);
	target_y_vec[0] = mean(y_vec);
  }
  return DUPM_CPv(degree_lpm, degree_upm, x_vec, y_vec, target_x_vec, target_y_vec);
}

//' Partial Moment Matrix
//'
//'
//' This function generates a co-partial moment matrix for the specified co-partial moment.
//' @param LPM_degree integer; Degree for \code{variable} below \code{target} deviations.  \code{(LPM_degree = 0)} is frequency, \code{(LPM_degree = 1)} is area.
//' @param UPM_degree integer; Degree for \code{variable} above \code{target} deviations.  \code{(UPM_degree = 0)} is frequency, \code{(UPM_degree = 1)} is area.
//' @param target numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)  \code{(target = NULL)} (default) will set the target as the mean of every variable.
//' @param variable a numeric matrix or data.frame.
//' @param pop_adj logical; \code{FALSE} (default) Adjusts the sample co-partial moment matrices for population statistics.
//' @return Matrix of partial moment quadrant values (CUPM, DUPM, DLPM, CLPM), and overall covariance matrix.  Uncalled quadrants will return a matrix of zeros.
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
//' PM.matrix(LPM_degree = 1, UPM_degree = 1, variable = A, target = colMeans(A), pop_adj = TRUE)
//'
//' ## Use of vectorized numeric targets (target_x, target_y, target_z)
//' PM.matrix(LPM_degree = 1, UPM_degree = 1, target = c(0, 0.15, .25), variable = A, pop_adj = TRUE)
//'
//' ## Calling Individual Partial Moment Quadrants
//' cov.mtx <- PM.matrix(LPM_degree = 1, UPM_degree = 1, variable = A, target = colMeans(A), 
//'                      pop_adj = TRUE)
//' cov.mtx$cupm
//'
//' ## Full covariance matrix
//' cov.mtx$cov.matrix
//' @export
// [[Rcpp::export("PM.matrix", rng = false)]]
List PMMatrix_RCPP(
    const double &LPM_degree,
    const double &UPM_degree,
    const RObject &target,
    const RObject &variable,
    const bool pop_adj
) {
  if(variable.isNULL()){
    Rcpp::stop("varible can't be null");
    return List::create();
  }
  NumericMatrix variable_matrix;
  if (is<NumericMatrix>(variable))
    variable_matrix = as<NumericMatrix>(variable);
  else if (is<IntegerMatrix>(variable))
    variable_matrix = as<NumericMatrix>(variable);
  else
    variable_matrix = Rcpp::internal::convert_using_rfunction(variable, "as.matrix");

  size_t variable_cols=variable_matrix.cols();
  NumericVector tgt;
  if((is<NumericVector>(target) || is<DataFrame>(target)) && !target.isNULL()){
      tgt=as<NumericVector>(target);
  }else{
      tgt=colMeans(variable_matrix);
  }
  
  size_t target_length=tgt.size();
  if(variable_cols != target_length){
    Rcpp::stop("varible matrix cols != target vector length");
    return List::create();
  }
  
  return PMMatrix_CPv(LPM_degree, UPM_degree, tgt, variable_matrix, pop_adj);
}
