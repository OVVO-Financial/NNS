// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "partial_moments.h"
/////////////////
// UPM / LPM
// single thread
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

// parallelFor
#define NNS_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS) \
  size_t target_size=target.size(); \
  NumericVector output = NumericVector(target_size); \
  WORKER_CLASS tmp_func(degree, target, variable, output); \
  parallelFor(0, target_size, tmp_func); \
  return(output);
NumericVector LPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(LPM_Worker);
}
NumericVector UPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  NNS_LPM_UPM_PARALLEL_FOR_FUNC(UPM_Worker);
}

NumericVector LPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  if (degree>0) {
    NumericVector lpm_output = LPM_CPv(degree, target, variable);
    NumericVector upm_output = UPM_CPv(degree, target, variable);
    NumericVector area = lpm_output+upm_output;
    return(lpm_output / area);
  } else {
    return LPM_CPv(degree, target, variable);
  }
}
NumericVector UPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable) {
  if (degree>0) {
    NumericVector lpm_output = LPM_CPv(degree, target, variable);
    NumericVector upm_output = UPM_CPv(degree, target, variable);
    NumericVector area = lpm_output+upm_output;
    return(upm_output / area);
  } else {
    return UPM_CPv(degree, target, variable);
  }
}

/////////////////
// CoUPM / CoLPM / DUPM / DLPM
// single thread
double CoUPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
){
  size_t n_x = x.size(), n_y = y.size();
  size_t max_size = (n_x>n_y ? n_x : n_y);
  size_t min_size = (n_x<n_y ? n_x : n_y);
  if (n_x != n_y)
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
  bool dont_use_pow_lpm=(degree_lpm==1 || degree_lpm==0), 
       dont_use_pow_upm=(degree_upm==1 || degree_upm==0),
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
  bool dont_use_pow_lpm=(degree_lpm==1 || degree_lpm==0), 
       dont_use_pow_upm=(degree_upm==1 || degree_upm==0),
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
// parallelFor
#define NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(WORKER_CLASS, LPM_DEGREE_VARIABLE, UPM_DEGREE_VARIABLE) \
  size_t target_x_size=target_x.size(); \
  size_t target_y_size=target_y.size(); \
  size_t max_target_size=(target_x_size>target_y_size?target_x_size:target_y_size); \
  NumericVector output = NumericVector(max_target_size); \
  WORKER_CLASS tmp_func(LPM_DEGREE_VARIABLE, UPM_DEGREE_VARIABLE, x, y, target_x, target_y, output); \
  parallelFor(0, output.size(), tmp_func); \
  return(output);
NumericVector CoLPM_CPv(
    const double &degree_lpm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoLPM_Worker, degree_lpm, degree_lpm);
}
NumericVector CoUPM_CPv(
    const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y 
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(CoUPM_Worker, degree_upm, degree_upm);
}
NumericVector DLPM_CPv(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(DLPM_Worker, degree_lpm, degree_upm);
}
NumericVector DUPM_CPv(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
) {
  NNS_CO_DE_LPM_UPM_PARALLEL_FOR_FUNC(DUPM_Worker, degree_lpm, degree_upm);
}

/////////////////
// PM MATRIX
// single thread
void PMMatrix_Cv(
  const double &degree_lpm, 
  const double &degree_upm, 
  const RMatrix<double>::Column &x, 
  const RMatrix<double>::Column &y, 
  const double &target_x,
  const double &target_y, 
  const bool &pop_adj, 
  const double &adjust,
  const size_t &rows, 
  double &coLpm,
  double &coUpm,   
  double &dLpm, 
  double &dUpm,
  double &covMat
){
  coLpm=0;
  coUpm=0;
  dLpm=0;
  dUpm=0; 
  covMat=0;
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
    double x_less_target_lpm = (x[i] - target_x),
           target_less_x_lpm = (target_x - x[i]),
           y_less_target_lpm = (y[i] - target_y),
           target_less_y_lpm = (target_y - y[i]);
    bool   less_equal_zero_target_less_x = (target_less_x_lpm <= 0),
           less_equal_zero_target_less_y = (target_less_y_lpm <= 0),
           less_equal_zero_x_less_target = (x_less_target_lpm <= 0),
           less_equal_zero_y_less_target = (y_less_target_lpm <= 0);
    if(
        less_equal_zero_target_less_x && 
        less_equal_zero_target_less_y && 
        less_equal_zero_x_less_target && 
        less_equal_zero_y_less_target
    )   // nothing to do with this row
        continue;
    // UPM values
    if(d_upm_0){ // pow(0,0) == 0, pow(!=0,0) == 1
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
    if(d_lpm_0){ // pow(0,0) == 0, pow(!=0,0) == 1
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
  dLpm /=rows;
  dUpm /=rows;
  if(pop_adj && rows>1){
    coLpm*=adjust;
    coUpm*=adjust;
    dLpm *=adjust;
    dUpm *=adjust;
  }
  covMat = coUpm + coLpm - dUpm - dLpm;
}

// parallelFor
List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const NumericVector &target,
    const NumericMatrix &variable,
    const bool &pop_adj
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
  NumericMatrix covMat(variable_cols, variable_cols);
  
  PMMatrix_Worker tmp_func(LPM_degree, UPM_degree, variable, target, pop_adj, coLpm, coUpm, dLpm, dUpm, covMat);
  parallelFor(0, variable_cols, tmp_func);
  
  rownames(coLpm) = colnames(variable);
  colnames(coLpm) = colnames(variable);
  
  rownames(coUpm) = colnames(variable);
  colnames(coUpm) = colnames(variable);
  
  rownames(dLpm) = colnames(variable);
  colnames(dLpm) = colnames(variable);
  
  rownames(dUpm) = colnames(variable);
  colnames(dUpm) = colnames(variable);
  
  rownames(covMat) = colnames(variable);
  colnames(covMat) = colnames(variable);
  
  return(
    List::create(
      Named("cupm") = coUpm,
      Named("dupm") = dUpm,
      Named("dlpm") = dLpm,
      Named("clpm") = coLpm,
      Named("cov.matrix") = covMat
    )
  );
}
