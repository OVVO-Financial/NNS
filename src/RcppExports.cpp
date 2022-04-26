// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LPM_C
double LPM_C(const double degree, const double target, const NumericVector variable);
RcppExport SEXP _NNS_LPM_C(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(LPM_C(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// UPM_C
double UPM_C(const double degree, const double target, const NumericVector variable);
RcppExport SEXP _NNS_UPM_C(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(UPM_C(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// LPM_CPv
NumericVector LPM_CPv(const double degree, const NumericVector target, const NumericVector variable);
RcppExport SEXP _NNS_LPM_CPv(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(LPM_CPv(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// UPM_CPv
NumericVector UPM_CPv(const double degree, const NumericVector target, const NumericVector variable);
RcppExport SEXP _NNS_UPM_CPv(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(UPM_CPv(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// CoUPM_C
double CoUPM_C(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const double target_x, const double target_y);
RcppExport SEXP _NNS_CoUPM_C(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const double >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(CoUPM_C(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// CoLPM_C
double CoLPM_C(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const double target_x, const double target_y);
RcppExport SEXP _NNS_CoLPM_C(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const double >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(CoLPM_C(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// DLPM_C
double DLPM_C(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const double target_x, const double target_y);
RcppExport SEXP _NNS_DLPM_C(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const double >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(DLPM_C(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// DUPM_C
double DUPM_C(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const double target_x, const double target_y);
RcppExport SEXP _NNS_DUPM_C(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const double >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(DUPM_C(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// CoLPM_CPv
NumericVector CoLPM_CPv(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const NumericVector target_x, const NumericVector target_y);
RcppExport SEXP _NNS_CoLPM_CPv(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(CoLPM_CPv(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// CoUPM_CPv
NumericVector CoUPM_CPv(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const NumericVector target_x, const NumericVector target_y);
RcppExport SEXP _NNS_CoUPM_CPv(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(CoUPM_CPv(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// DLPM_CPv
NumericVector DLPM_CPv(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const NumericVector target_x, const NumericVector target_y);
RcppExport SEXP _NNS_DLPM_CPv(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(DLPM_CPv(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// DUPM_CPv
NumericVector DUPM_CPv(const double degree_x, const double degree_y, const NumericVector x, const NumericVector y, const NumericVector target_x, const NumericVector target_y);
RcppExport SEXP _NNS_DUPM_CPv(SEXP degree_xSEXP, SEXP degree_ySEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type degree_x(degree_xSEXP);
    Rcpp::traits::input_parameter< const double >::type degree_y(degree_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(DUPM_CPv(degree_x, degree_y, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NNS_LPM_C", (DL_FUNC) &_NNS_LPM_C, 3},
    {"_NNS_UPM_C", (DL_FUNC) &_NNS_UPM_C, 3},
    {"_NNS_LPM_CPv", (DL_FUNC) &_NNS_LPM_CPv, 3},
    {"_NNS_UPM_CPv", (DL_FUNC) &_NNS_UPM_CPv, 3},
    {"_NNS_CoUPM_C", (DL_FUNC) &_NNS_CoUPM_C, 6},
    {"_NNS_CoLPM_C", (DL_FUNC) &_NNS_CoLPM_C, 6},
    {"_NNS_DLPM_C", (DL_FUNC) &_NNS_DLPM_C, 6},
    {"_NNS_DUPM_C", (DL_FUNC) &_NNS_DUPM_C, 6},
    {"_NNS_CoLPM_CPv", (DL_FUNC) &_NNS_CoLPM_CPv, 6},
    {"_NNS_CoUPM_CPv", (DL_FUNC) &_NNS_CoUPM_CPv, 6},
    {"_NNS_DLPM_CPv", (DL_FUNC) &_NNS_DLPM_CPv, 6},
    {"_NNS_DUPM_CPv", (DL_FUNC) &_NNS_DUPM_CPv, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_NNS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
