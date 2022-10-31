// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LPM_RCPP
NumericVector LPM_RCPP(const double& degree, const RObject& target, const RObject& variable);
RcppExport SEXP _NNS_LPM_RCPP(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(LPM_RCPP(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// UPM_RCPP
NumericVector UPM_RCPP(const double& degree, const RObject& target, const RObject& variable);
RcppExport SEXP _NNS_UPM_RCPP(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(UPM_RCPP(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// LPM_ratio_RCPP
NumericVector LPM_ratio_RCPP(const double& degree, const RObject& target, const RObject& variable);
RcppExport SEXP _NNS_LPM_ratio_RCPP(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(LPM_ratio_RCPP(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// UPM_ratio_RCPP
NumericVector UPM_ratio_RCPP(const double& degree, const RObject& target, const RObject& variable);
RcppExport SEXP _NNS_UPM_ratio_RCPP(SEXP degreeSEXP, SEXP targetSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type variable(variableSEXP);
    rcpp_result_gen = Rcpp::wrap(UPM_ratio_RCPP(degree, target, variable));
    return rcpp_result_gen;
END_RCPP
}
// CoLPM_RCPP
NumericVector CoLPM_RCPP(const double& degree_lpm, const RObject& x, const RObject& y, const RObject& target_x, const RObject& target_y);
RcppExport SEXP _NNS_CoLPM_RCPP(SEXP degree_lpmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree_lpm(degree_lpmSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(CoLPM_RCPP(degree_lpm, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// CoUPM_RCPP
NumericVector CoUPM_RCPP(const double& degree_upm, const RObject& x, const RObject& y, const RObject& target_x, const RObject& target_y);
RcppExport SEXP _NNS_CoUPM_RCPP(SEXP degree_upmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree_upm(degree_upmSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(CoUPM_RCPP(degree_upm, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// DLPM_RCPP
NumericVector DLPM_RCPP(const double& degree_lpm, const double& degree_upm, const RObject& x, const RObject& y, const RObject& target_x, const RObject& target_y);
RcppExport SEXP _NNS_DLPM_RCPP(SEXP degree_lpmSEXP, SEXP degree_upmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree_lpm(degree_lpmSEXP);
    Rcpp::traits::input_parameter< const double& >::type degree_upm(degree_upmSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(DLPM_RCPP(degree_lpm, degree_upm, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// DUPM_RCPP
NumericVector DUPM_RCPP(const double& degree_lpm, const double& degree_upm, const RObject& x, const RObject& y, const RObject& target_x, const RObject& target_y);
RcppExport SEXP _NNS_DUPM_RCPP(SEXP degree_lpmSEXP, SEXP degree_upmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP target_xSEXP, SEXP target_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type degree_lpm(degree_lpmSEXP);
    Rcpp::traits::input_parameter< const double& >::type degree_upm(degree_upmSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_x(target_xSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target_y(target_ySEXP);
    rcpp_result_gen = Rcpp::wrap(DUPM_RCPP(degree_lpm, degree_upm, x, y, target_x, target_y));
    return rcpp_result_gen;
END_RCPP
}
// PMMatrix_RCPP
List PMMatrix_RCPP(const double& LPM_degree, const double& UPM_degree, const RObject& target, const RObject& variable, const bool pop_adj);
RcppExport SEXP _NNS_PMMatrix_RCPP(SEXP LPM_degreeSEXP, SEXP UPM_degreeSEXP, SEXP targetSEXP, SEXP variableSEXP, SEXP pop_adjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const double& >::type LPM_degree(LPM_degreeSEXP);
    Rcpp::traits::input_parameter< const double& >::type UPM_degree(UPM_degreeSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type variable(variableSEXP);
    Rcpp::traits::input_parameter< const bool >::type pop_adj(pop_adjSEXP);
    rcpp_result_gen = Rcpp::wrap(PMMatrix_RCPP(LPM_degree, UPM_degree, target, variable, pop_adj));
    return rcpp_result_gen;
END_RCPP
}
// NNS_bin
List NNS_bin(NumericVector x, double width, double origin, bool missinglast);
RcppExport SEXP _NNS_NNS_bin(SEXP xSEXP, SEXP widthSEXP, SEXP originSEXP, SEXP missinglastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type width(widthSEXP);
    Rcpp::traits::input_parameter< double >::type origin(originSEXP);
    Rcpp::traits::input_parameter< bool >::type missinglast(missinglastSEXP);
    rcpp_result_gen = Rcpp::wrap(NNS_bin(x, width, origin, missinglast));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NNS_LPM_RCPP", (DL_FUNC) &_NNS_LPM_RCPP, 3},
    {"_NNS_UPM_RCPP", (DL_FUNC) &_NNS_UPM_RCPP, 3},
    {"_NNS_LPM_ratio_RCPP", (DL_FUNC) &_NNS_LPM_ratio_RCPP, 3},
    {"_NNS_UPM_ratio_RCPP", (DL_FUNC) &_NNS_UPM_ratio_RCPP, 3},
    {"_NNS_CoLPM_RCPP", (DL_FUNC) &_NNS_CoLPM_RCPP, 5},
    {"_NNS_CoUPM_RCPP", (DL_FUNC) &_NNS_CoUPM_RCPP, 5},
    {"_NNS_DLPM_RCPP", (DL_FUNC) &_NNS_DLPM_RCPP, 6},
    {"_NNS_DUPM_RCPP", (DL_FUNC) &_NNS_DUPM_RCPP, 6},
    {"_NNS_PMMatrix_RCPP", (DL_FUNC) &_NNS_PMMatrix_RCPP, 5},
    {"_NNS_NNS_bin", (DL_FUNC) &_NNS_NNS_bin, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_NNS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
