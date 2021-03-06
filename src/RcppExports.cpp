// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// declusterCpp
int declusterCpp(NumericVector complete_series, NumericVector y, double series_mean);
RcppExport SEXP _potMax_declusterCpp(SEXP complete_seriesSEXP, SEXP ySEXP, SEXP series_meanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type complete_series(complete_seriesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type series_mean(series_meanSEXP);
    rcpp_result_gen = Rcpp::wrap(declusterCpp(complete_series, y, series_mean));
    return rcpp_result_gen;
END_RCPP
}
// declusterWithTimeCpp
int declusterWithTimeCpp(NumericVector complete_series, NumericVector obs_times, NumericVector y, NumericVector t, double series_mean);
RcppExport SEXP _potMax_declusterWithTimeCpp(SEXP complete_seriesSEXP, SEXP obs_timesSEXP, SEXP ySEXP, SEXP tSEXP, SEXP series_meanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type complete_series(complete_seriesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type obs_times(obs_timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type series_mean(series_meanSEXP);
    rcpp_result_gen = Rcpp::wrap(declusterWithTimeCpp(complete_series, obs_times, y, t, series_mean));
    return rcpp_result_gen;
END_RCPP
}
// gumbelMaxDistCpp
NumericVector gumbelMaxDistCpp(double mu, double sigma, double Lambda, double integration_constant, int n_mc, bool progress_tf);
RcppExport SEXP _potMax_gumbelMaxDistCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP LambdaSEXP, SEXP integration_constantSEXP, SEXP n_mcSEXP, SEXP progress_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< double >::type integration_constant(integration_constantSEXP);
    Rcpp::traits::input_parameter< int >::type n_mc(n_mcSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_tf(progress_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(gumbelMaxDistCpp(mu, sigma, Lambda, integration_constant, n_mc, progress_tf));
    return rcpp_result_gen;
END_RCPP
}
// gumbelMaxDistUncertCpp
NumericMatrix gumbelMaxDistUncertCpp(NumericVector mu, NumericVector sigma, NumericVector Lambda, NumericVector integration_constant, int n_mc, int n_boot, bool progress_tf);
RcppExport SEXP _potMax_gumbelMaxDistUncertCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP LambdaSEXP, SEXP integration_constantSEXP, SEXP n_mcSEXP, SEXP n_bootSEXP, SEXP progress_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type integration_constant(integration_constantSEXP);
    Rcpp::traits::input_parameter< int >::type n_mc(n_mcSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_tf(progress_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(gumbelMaxDistUncertCpp(mu, sigma, Lambda, integration_constant, n_mc, n_boot, progress_tf));
    return rcpp_result_gen;
END_RCPP
}
// fullMaxDistCpp
NumericVector fullMaxDistCpp(double mu, double sigma, double k, double Lambda, double integration_constant, int n_mc, bool progress_tf);
RcppExport SEXP _potMax_fullMaxDistCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP kSEXP, SEXP LambdaSEXP, SEXP integration_constantSEXP, SEXP n_mcSEXP, SEXP progress_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< double >::type integration_constant(integration_constantSEXP);
    Rcpp::traits::input_parameter< int >::type n_mc(n_mcSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_tf(progress_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(fullMaxDistCpp(mu, sigma, k, Lambda, integration_constant, n_mc, progress_tf));
    return rcpp_result_gen;
END_RCPP
}
// fullMaxDistUncertCpp
NumericMatrix fullMaxDistUncertCpp(NumericVector mu, NumericVector sigma, NumericVector k, NumericVector Lambda, NumericVector integration_constant, int n_mc, int n_boot, bool progress_tf);
RcppExport SEXP _potMax_fullMaxDistUncertCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP kSEXP, SEXP LambdaSEXP, SEXP integration_constantSEXP, SEXP n_mcSEXP, SEXP n_bootSEXP, SEXP progress_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type integration_constant(integration_constantSEXP);
    Rcpp::traits::input_parameter< int >::type n_mc(n_mcSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_tf(progress_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(fullMaxDistUncertCpp(mu, sigma, k, Lambda, integration_constant, n_mc, n_boot, progress_tf));
    return rcpp_result_gen;
END_RCPP
}
// gumbelMaxDistMultiCpp
NumericVector gumbelMaxDistMultiCpp(NumericVector mu, NumericVector sigma, NumericVector Lambda, NumericVector integration_constant, int n_mc, bool progress_tf);
RcppExport SEXP _potMax_gumbelMaxDistMultiCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP LambdaSEXP, SEXP integration_constantSEXP, SEXP n_mcSEXP, SEXP progress_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type integration_constant(integration_constantSEXP);
    Rcpp::traits::input_parameter< int >::type n_mc(n_mcSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_tf(progress_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(gumbelMaxDistMultiCpp(mu, sigma, Lambda, integration_constant, n_mc, progress_tf));
    return rcpp_result_gen;
END_RCPP
}
// fullMaxDistMultiCpp
NumericVector fullMaxDistMultiCpp(NumericVector mu, NumericVector sigma, NumericVector k, NumericVector Lambda, NumericVector integration_constant, int n_mc, bool progress_tf);
RcppExport SEXP _potMax_fullMaxDistMultiCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP kSEXP, SEXP LambdaSEXP, SEXP integration_constantSEXP, SEXP n_mcSEXP, SEXP progress_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type integration_constant(integration_constantSEXP);
    Rcpp::traits::input_parameter< int >::type n_mc(n_mcSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_tf(progress_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(fullMaxDistMultiCpp(mu, sigma, k, Lambda, integration_constant, n_mc, progress_tf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_potMax_declusterCpp", (DL_FUNC) &_potMax_declusterCpp, 3},
    {"_potMax_declusterWithTimeCpp", (DL_FUNC) &_potMax_declusterWithTimeCpp, 5},
    {"_potMax_gumbelMaxDistCpp", (DL_FUNC) &_potMax_gumbelMaxDistCpp, 6},
    {"_potMax_gumbelMaxDistUncertCpp", (DL_FUNC) &_potMax_gumbelMaxDistUncertCpp, 7},
    {"_potMax_fullMaxDistCpp", (DL_FUNC) &_potMax_fullMaxDistCpp, 7},
    {"_potMax_fullMaxDistUncertCpp", (DL_FUNC) &_potMax_fullMaxDistUncertCpp, 8},
    {"_potMax_gumbelMaxDistMultiCpp", (DL_FUNC) &_potMax_gumbelMaxDistMultiCpp, 6},
    {"_potMax_fullMaxDistMultiCpp", (DL_FUNC) &_potMax_fullMaxDistMultiCpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_potMax(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
