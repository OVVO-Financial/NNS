#ifndef NNS_partial_moments_H
#define NNS_partial_moments_H

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <functional>
using namespace Rcpp;
using namespace RcppParallel;

/////////////////
// UPM / LPM
// single thread
double LPM_C(const double &degree, const double &target, const RVector<double> &variable);
double UPM_C(const double &degree, const double &target, const RVector<double> &variable);
// parallelFor
#define NNS_PM_SINGLE_VARIABLE_WORKER(NAME, FUNC) \
  struct NAME : public Worker\
  {\
    const double degree;\
    const RVector<double> target;\
    const RVector<double> variable;\
    RVector<double> output;\
    NAME (\
      const double degree,\
      const NumericVector &target,\
      const NumericVector &variable,\
      NumericVector &output\
    ): degree(degree), target(target), variable(variable), output(output) {}\
    void operator()(std::size_t begin, std::size_t end) {\
      for (size_t i = begin; i < end; i++)\
        output[i] = FUNC(degree, target[i], variable); \
    }\
  }
NNS_PM_SINGLE_VARIABLE_WORKER(LPM_Worker, LPM_C);
NNS_PM_SINGLE_VARIABLE_WORKER(UPM_Worker, UPM_C);
NumericVector LPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable);
NumericVector UPM_CPv(const double &degree, const NumericVector &target, const NumericVector &variable);
NumericVector LPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable);
NumericVector UPM_ratio_CPv(const double &degree, const NumericVector &target, const NumericVector &variable);


/////////////////
// CoUPM / CoLPM / DUPM / DLPM
// single thread
double CoUPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
);
double CoLPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y 
);
double DLPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y
);
double DUPM_C(
    const double &degree_lpm, const double &degree_upm, 
    const RVector<double> &x, const RVector<double> &y, 
    const double &target_x, const double &target_y
);

// parallelFor
#define NNS_PM_TWO_VARIABLES_WORKER(NAME, FUNC) \
  struct NAME : public Worker\
  {\
    const double degree_lpm;\
    const double degree_upm;\
    const RVector<double> x;\
    const RVector<double> y;\
    const RVector<double> target_x;\
    const RVector<double> target_y;\
    const size_t n_t_x;\
    const size_t n_t_y;\
    RVector<double> output;\
    NAME (\
      const double degree_lpm,\
      const double degree_upm,\
      const NumericVector &x, const NumericVector &y,\
      const NumericVector &target_x, const NumericVector &target_y,\
      NumericVector output\
    ):\
      degree_lpm(degree_lpm), degree_upm(degree_upm),  \
	  x(x), y(y), target_x(target_x), target_y(target_y),\
      n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)\
    {}\
    void operator()(std::size_t begin, std::size_t end) { \
      for (size_t i = begin; i < end; i++) {\
       output[i] = FUNC(degree_lpm, degree_upm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); \
     } \
    }\
  }
  
NNS_PM_TWO_VARIABLES_WORKER(CoLPM_Worker, CoLPM_C);
NNS_PM_TWO_VARIABLES_WORKER(CoUPM_Worker, CoUPM_C);
NNS_PM_TWO_VARIABLES_WORKER(DLPM_Worker, DLPM_C);
NNS_PM_TWO_VARIABLES_WORKER(DUPM_Worker, DUPM_C);
NumericVector CoLPM_CPv(
    const double &degree_lpm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
);
NumericVector CoUPM_CPv(
    const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y 
);
NumericVector DLPM_CPv(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
);
NumericVector DUPM_CPv(
    const double &degree_lpm, const double &degree_upm, 
    const NumericVector &x, const NumericVector &y, 
    const NumericVector &target_x, const NumericVector &target_y
);

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
);
// parallelFor
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
  double adjust;
  RMatrix<double> coLpm;
  RMatrix<double> coUpm;
  RMatrix<double> dLpm;
  RMatrix<double> dUpm;
  RMatrix<double> covMat;
  PMMatrix_Worker(
    const double &degree_lpm, const double &degree_upm, 
    const NumericMatrix &variable, 
    const NumericVector &target,
    const bool &pop_adj,
    NumericMatrix &coLpm, NumericMatrix &coUpm,
    NumericMatrix &dLpm,  NumericMatrix &dUpm,
    NumericMatrix &covMat
  ): 
    degree_lpm(degree_lpm), degree_upm(degree_upm), 
    variable(variable), target(target),
    variable_cols(variable.cols()), variable_rows(variable.rows()), target_length(target.size()), 
    pop_adj(pop_adj),
    coLpm(coLpm), coUpm(coUpm), 
    dLpm(dLpm), dUpm(dUpm),
    covMat(covMat)
  {
    if(variable_cols != target_length)
      Rcpp::stop("varible matrix cols != target vector length");
    adjust = 1;
    if (variable_rows > 1)
	  adjust=((double)variable_rows)/((double)variable_rows-1);
  }
  void operator()(std::size_t begin, std::size_t end) {
	for (size_t i = begin; i < end; i++){
      for (size_t l = 0; l < variable_cols; l++){
        PMMatrix_Cv(
          degree_lpm,
          degree_upm,
          variable.column(i),
          variable.column(l),
          target[i],
          target[l],
          pop_adj,
          adjust,
          variable_rows,
          coLpm(i,l),
          coUpm(i,l),
          dLpm(i,l),
          dUpm(i,l),
          covMat(i,l)
        );
      }
    }
  }
};
List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const NumericVector &target,
    const NumericMatrix &variable,
    const bool &pop_adj
);
#endif  //NNS_partial_moments_H
