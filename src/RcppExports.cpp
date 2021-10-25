// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logCLnormal_par
double logCLnormal_par(const NumericVector& beta, const NumericMatrix& cov, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_logCLnormal_par(SEXP betaSEXP, SEXP covSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(logCLnormal_par(beta, cov, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// eenormalmean_par
List eenormalmean_par(const NumericVector& beta, const NumericMatrix& cov, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_eenormalmean_par(SEXP betaSEXP, SEXP covSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(eenormalmean_par(beta, cov, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// GEEnormalmean_par
List GEEnormalmean_par(const NumericVector& beta, const NumericMatrix& cov_inv, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_GEEnormalmean_par(SEXP betaSEXP, SEXP cov_invSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cov_inv(cov_invSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(GEEnormalmean_par(beta, cov_inv, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// GEEnormalderivmean_par
List GEEnormalderivmean_par(const NumericMatrix& cov_inv, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_GEEnormalderivmean_par(SEXP cov_invSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cov_inv(cov_invSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(GEEnormalderivmean_par(cov_inv, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// eenormalderivmean_par
List eenormalderivmean_par(const NumericMatrix& cov, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_eenormalderivmean_par(SEXP covSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(eenormalderivmean_par(cov, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// psi_g_AR1derivvar_par
List psi_g_AR1derivvar_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_psi_g_AR1derivvar_par(SEXP betaSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_g_AR1derivvar_par(beta, sigma, rho, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// psi_g_CSderivvar_par
List psi_g_CSderivvar_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_psi_g_CSderivvar_par(SEXP betaSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_g_CSderivvar_par(beta, sigma, rho, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// psi_g_indderivvar_par
List psi_g_indderivvar_par(const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_psi_g_indderivvar_par(SEXP betaSEXP, SEXP sigmaSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_g_indderivvar_par(beta, sigma, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// eenormalCSvar_par
List eenormalCSvar_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_eenormalCSvar_par(SEXP betaSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(eenormalCSvar_par(beta, sigma, rho, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// eenormalAR1var_par
List eenormalAR1var_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_eenormalAR1var_par(SEXP betaSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(eenormalAR1var_par(beta, sigma, rho, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// eenormalindvar_par
List eenormalindvar_par(const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, const NumericMatrix& block_x, const double& m, const double& n);
RcppExport SEXP _DIQIF_eenormalindvar_par(SEXP betaSEXP, SEXP sigmaSEXP, SEXP block_ySEXP, SEXP block_xSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_y(block_ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type block_x(block_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(eenormalindvar_par(beta, sigma, block_y, block_x, m, n));
    return rcpp_result_gen;
END_RCPP
}
// increQIF_sub
List increQIF_sub(arma::mat X, arma::vec y, arma::vec nobs, String family, String corstr, arma::vec beta_old, int maxit, double tol);
RcppExport SEXP _DIQIF_increQIF_sub(SEXP XSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP familySEXP, SEXP corstrSEXP, SEXP beta_oldSEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< String >::type corstr(corstrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(increQIF_sub(X, y, nobs, family, corstr, beta_old, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}
// QIF_eval
List QIF_eval(arma::mat X, arma::vec y, arma::vec nobs, String family, String corstr, arma::vec beta_old);
RcppExport SEXP _DIQIF_QIF_eval(SEXP XSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP familySEXP, SEXP corstrSEXP, SEXP beta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< String >::type corstr(corstrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(QIF_eval(X, y, nobs, family, corstr, beta_old));
    return rcpp_result_gen;
END_RCPP
}
// getHessGrads
Rcpp::List getHessGrads(SEXP y, SEXP x, SEXP offset, SEXP doffset, SEXP w, SEXP linkwave, SEXP zsca, SEXP zcor, SEXP corp, SEXP clusz, SEXP geestr, SEXP cor, SEXP par);
RcppExport SEXP _DIQIF_getHessGrads(SEXP ySEXP, SEXP xSEXP, SEXP offsetSEXP, SEXP doffsetSEXP, SEXP wSEXP, SEXP linkwaveSEXP, SEXP zscaSEXP, SEXP zcorSEXP, SEXP corpSEXP, SEXP cluszSEXP, SEXP geestrSEXP, SEXP corSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type doffset(doffsetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type w(wSEXP);
    Rcpp::traits::input_parameter< SEXP >::type linkwave(linkwaveSEXP);
    Rcpp::traits::input_parameter< SEXP >::type zsca(zscaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type zcor(zcorSEXP);
    Rcpp::traits::input_parameter< SEXP >::type corp(corpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type clusz(cluszSEXP);
    Rcpp::traits::input_parameter< SEXP >::type geestr(geestrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cor(corSEXP);
    Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(getHessGrads(y, x, offset, doffset, w, linkwave, zsca, zcor, corp, clusz, geestr, cor, par));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DIQIF_logCLnormal_par", (DL_FUNC) &_DIQIF_logCLnormal_par, 6},
    {"_DIQIF_eenormalmean_par", (DL_FUNC) &_DIQIF_eenormalmean_par, 6},
    {"_DIQIF_GEEnormalmean_par", (DL_FUNC) &_DIQIF_GEEnormalmean_par, 6},
    {"_DIQIF_GEEnormalderivmean_par", (DL_FUNC) &_DIQIF_GEEnormalderivmean_par, 5},
    {"_DIQIF_eenormalderivmean_par", (DL_FUNC) &_DIQIF_eenormalderivmean_par, 5},
    {"_DIQIF_psi_g_AR1derivvar_par", (DL_FUNC) &_DIQIF_psi_g_AR1derivvar_par, 7},
    {"_DIQIF_psi_g_CSderivvar_par", (DL_FUNC) &_DIQIF_psi_g_CSderivvar_par, 7},
    {"_DIQIF_psi_g_indderivvar_par", (DL_FUNC) &_DIQIF_psi_g_indderivvar_par, 6},
    {"_DIQIF_eenormalCSvar_par", (DL_FUNC) &_DIQIF_eenormalCSvar_par, 7},
    {"_DIQIF_eenormalAR1var_par", (DL_FUNC) &_DIQIF_eenormalAR1var_par, 7},
    {"_DIQIF_eenormalindvar_par", (DL_FUNC) &_DIQIF_eenormalindvar_par, 6},
    {"_DIQIF_increQIF_sub", (DL_FUNC) &_DIQIF_increQIF_sub, 8},
    {"_DIQIF_QIF_eval", (DL_FUNC) &_DIQIF_QIF_eval, 6},
    {"_DIQIF_getHessGrads", (DL_FUNC) &_DIQIF_getHessGrads, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_DIQIF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
