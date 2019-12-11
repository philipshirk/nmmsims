// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pbinsumRow
NumericVector pbinsumRow(NumericVector y, double N, NumericVector p);
RcppExport SEXP _nmmsims_pbinsumRow(SEXP ySEXP, SEXP NSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pbinsumRow(y, N, p));
    return rcpp_result_gen;
END_RCPP
}
// pbinsum
NumericMatrix pbinsum(NumericMatrix y, NumericVector N, NumericMatrix p);
RcppExport SEXP _nmmsims_pbinsum(SEXP ySEXP, SEXP NSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pbinsum(y, N, p));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _nmmsims_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nmmsims_pbinsumRow", (DL_FUNC) &_nmmsims_pbinsumRow, 3},
    {"_nmmsims_pbinsum", (DL_FUNC) &_nmmsims_pbinsum, 3},
    {"_nmmsims_rcpp_hello_world", (DL_FUNC) &_nmmsims_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_nmmsims(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
