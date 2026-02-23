// partial_moments.h
#ifndef NNS_partial_moments_H
#define NNS_partial_moments_H

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

// Backend API for the partial moment computations. These routines operate on
// RcppParallel vector and matrix proxies so they can be reused from serial and
// parallel workers. Higher-level wrappers that accept generic R objects live in
// partial_moments_rcpp.h/cpp.

/////////////////
// UPM / LPM
// single thread
double LPM_C(const double &degree,
             const double &target,
             const RcppParallel::RVector<double> &variable);
double UPM_C(const double &degree,
             const double &target,
             const RcppParallel::RVector<double> &variable);
// parallelFor
#define NNS_PM_SINGLE_VARIABLE_WORKER(NAME, FUNC)                         \
struct NAME : public RcppParallel::Worker                                 \
{                                                                         \
  const double degree;                                                    \
  const RcppParallel::RVector<double> target;                             \
  const RcppParallel::RVector<double> variable;                           \
  RcppParallel::RVector<double> output;                                   \
  NAME (                                                                  \
      const double degree,                                                \
      const Rcpp::NumericVector &target,                                  \
      const Rcpp::NumericVector &variable,                                \
      Rcpp::NumericVector &output                                         \
  ): degree(degree), target(target), variable(variable), output(output) {}\
  void operator()(std::size_t begin, std::size_t end) {                   \
    for (size_t i = begin; i < end; i++)                                  \
      output[i] = FUNC(degree, target[i], variable);                      \
  }                                                                       \
}
NNS_PM_SINGLE_VARIABLE_WORKER(LPM_Worker, LPM_C);
NNS_PM_SINGLE_VARIABLE_WORKER(UPM_Worker, UPM_C);
Rcpp::NumericVector LPM_CPv(const double &degree,
                            const Rcpp::NumericVector &target,
                            const Rcpp::NumericVector &variable);
Rcpp::NumericVector UPM_CPv(const double &degree,
                            const Rcpp::NumericVector &target,
                            const Rcpp::NumericVector &variable);
Rcpp::NumericVector LPM_ratio_CPv(const double &degree,
                                  const Rcpp::NumericVector &target,
                                  const Rcpp::NumericVector &variable);
Rcpp::NumericVector UPM_ratio_CPv(const double &degree,
                                  const Rcpp::NumericVector &target,
                                  const Rcpp::NumericVector &variable);

/////////////////
// CoUPM / CoLPM / DUPM / DLPM
// single thread
double CoUPM_C(
    const double &degree_x, const double &degree_y,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);
double CoLPM_C(
    const double &degree_x, const double &degree_y,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);
double DLPM_C(
    const double &degree_lpm, const double &degree_upm,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);
double DUPM_C(
    const double &degree_lpm, const double &degree_upm,
    const RcppParallel::RVector<double> &x, const RcppParallel::RVector<double> &y,
    const double &target_x, const double &target_y
);

// parallelFor
#define NNS_PM_TWO_VARIABLES_WORKER(NAME, FUNC)                                             \
struct NAME : public RcppParallel::Worker                                                   \
{                                                                                           \
  const double degree_lpm;                                                                  \
  const double degree_upm;                                                                  \
  const RcppParallel::RVector<double> x;                                                    \
  const RcppParallel::RVector<double> y;                                                    \
  const RcppParallel::RVector<double> target_x;                                             \
  const RcppParallel::RVector<double> target_y;                                             \
  const size_t n_t_x;                                                                       \
  const size_t n_t_y;                                                                       \
  RcppParallel::RVector<double> output;                                                     \
  NAME (                                                                                    \
      const double degree_lpm,                                                              \
      const double degree_upm,                                                              \
      const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,                           \
      const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y,             \
      Rcpp::NumericVector &output                                                           \
  ):                                                                                        \
    degree_lpm(degree_lpm), degree_upm(degree_upm),                                         \
    x(x), y(y), target_x(target_x), target_y(target_y),                                     \
    n_t_x(target_x.size()), n_t_y(target_y.size()), output(output)                          \
  {}                                                                                        \
  void operator()(std::size_t begin, std::size_t end) {                                     \
    for (size_t i = begin; i < end; i++) {                                                  \
      output[i] = FUNC(degree_lpm, degree_upm, x, y, target_x[i%n_t_x], target_y[i%n_t_y]); \
    }                                                                                       \
  }                                                                                         \
}

NNS_PM_TWO_VARIABLES_WORKER(CoLPM_Worker, CoLPM_C);
NNS_PM_TWO_VARIABLES_WORKER(CoUPM_Worker, CoUPM_C);
NNS_PM_TWO_VARIABLES_WORKER(DLPM_Worker, DLPM_C);
NNS_PM_TWO_VARIABLES_WORKER(DUPM_Worker, DUPM_C);
Rcpp::NumericVector CoLPM_CPv(
    const double &degree_x, const double &degree_y,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);
Rcpp::NumericVector CoUPM_CPv(
    const double &degree_x, const double &degree_y,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);
Rcpp::NumericVector DLPM_CPv(
    const double &degree_lpm, const double &degree_upm,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);
Rcpp::NumericVector DUPM_CPv(
    const double &degree_lpm, const double &degree_upm,
    const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
    const Rcpp::NumericVector &target_x, const Rcpp::NumericVector &target_y
);

/////////////////
// PM MATRIX
// single thread
void PMMatrix_Cv(
    const double &degree_lpm,
    const double &degree_upm,
    const RcppParallel::RMatrix<double>::Column &x,
    const RcppParallel::RMatrix<double>::Column &y,
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
struct PMMatrix_Worker : public RcppParallel::Worker
{
  const double degree_lpm;
  const double degree_upm;
  const RcppParallel::RMatrix<double> variable;
  const RcppParallel::RVector<double> target;
  const size_t variable_cols;
  const size_t variable_rows;
  const size_t target_length;
  const bool pop_adj;
  double adjust;
  RcppParallel::RMatrix<double> coLpm;
  RcppParallel::RMatrix<double> coUpm;
  RcppParallel::RMatrix<double> dLpm;
  RcppParallel::RMatrix<double> dUpm;
  RcppParallel::RMatrix<double> covMat;
  PMMatrix_Worker(
    const double &degree_lpm, const double &degree_upm,
    const Rcpp::NumericMatrix &variable,
    const Rcpp::NumericVector &target,
    const bool &pop_adj,
    Rcpp::NumericMatrix &coLpm, Rcpp::NumericMatrix &coUpm,
    Rcpp::NumericMatrix &dLpm,  Rcpp::NumericMatrix &dUpm,
    Rcpp::NumericMatrix &covMat
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
      Rcpp::stop("variable matrix cols != target vector length");
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
Rcpp::List PMMatrix_CPv(
    const double &LPM_degree,
    const double &UPM_degree,
    const Rcpp::NumericVector &target,
    const Rcpp::NumericMatrix &variable,
    const bool &pop_adj,
    const bool &norm
);

// n‐D co‐partial‐moments prototypes (parallel back‐ends)
double clpm_nD_cpp(const Rcpp::NumericMatrix &data,
                   const Rcpp::NumericVector &target,
                   double degree,
                   bool norm);

double cupm_nD_cpp(const Rcpp::NumericMatrix &data,
                   const Rcpp::NumericVector &target,
                   double degree,
                   bool norm);

double dpm_nD_cpp(const Rcpp::NumericMatrix &data,
                  const Rcpp::NumericVector &target,
                  double degree,
                  bool norm);

#endif  //NNS_partial_moments_H
