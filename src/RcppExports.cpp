// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rou_y
arma::vec rou_y(int N, double sd_u);
RcppExport SEXP VARext_rou_y(SEXP NSEXP, SEXP sd_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type sd_u(sd_uSEXP);
    rcpp_result_gen = Rcpp::wrap(rou_y(N, sd_u));
    return rcpp_result_gen;
END_RCPP
}
// M_i_ig
arma::vec M_i_ig(arma::mat Y, int i, arma::vec a, arma::mat A, arma::mat Sigma, int print_level);
RcppExport SEXP VARext_M_i_ig(SEXP YSEXP, SEXP iSEXP, SEXP aSEXP, SEXP ASEXP, SEXP SigmaSEXP, SEXP print_levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    rcpp_result_gen = Rcpp::wrap(M_i_ig(Y, i, a, A, Sigma, print_level));
    return rcpp_result_gen;
END_RCPP
}
// var_cm
arma::mat var_cm(arma::mat Y, arma::vec a, arma::mat A);
RcppExport SEXP VARext_var_cm(SEXP YSEXP, SEXP aSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(var_cm(Y, a, A));
    return rcpp_result_gen;
END_RCPP
}
// disc_cv
arma::mat disc_cv(arma::mat Y, arma::mat M);
RcppExport SEXP VARext_disc_cv(SEXP YSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(disc_cv(Y, M));
    return rcpp_result_gen;
END_RCPP
}
// disc_skew
arma::mat disc_skew(arma::mat Y, arma::mat M);
RcppExport SEXP VARext_disc_skew(SEXP YSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(disc_skew(Y, M));
    return rcpp_result_gen;
END_RCPP
}
// disc_obj
double disc_obj(arma::vec m_i, arma::mat Y, int i, arma::vec a, arma::mat A, arma::mat Sigma, arma::vec w, int print_level);
RcppExport SEXP VARext_disc_obj(SEXP m_iSEXP, SEXP YSEXP, SEXP iSEXP, SEXP aSEXP, SEXP ASEXP, SEXP SigmaSEXP, SEXP wSEXP, SEXP print_levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type m_i(m_iSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    rcpp_result_gen = Rcpp::wrap(disc_obj(m_i, Y, i, a, A, Sigma, w, print_level));
    return rcpp_result_gen;
END_RCPP
}
// disc_grad
arma::vec disc_grad(arma::vec m_i, arma::mat Y, int i, arma::vec a, arma::mat A, arma::mat Sigma, arma::vec w);
RcppExport SEXP VARext_disc_grad(SEXP m_iSEXP, SEXP YSEXP, SEXP iSEXP, SEXP aSEXP, SEXP ASEXP, SEXP SigmaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type m_i(m_iSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(disc_grad(m_i, Y, i, a, A, Sigma, w));
    return rcpp_result_gen;
END_RCPP
}
// var_lhood_N
double var_lhood_N(arma::mat Y, arma::vec par, int lags, int print_level, bool cond);
RcppExport SEXP VARext_var_lhood_N(SEXP YSEXP, SEXP parSEXP, SEXP lagsSEXP, SEXP print_levelSEXP, SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lhood_N(Y, par, lags, print_level, cond));
    return rcpp_result_gen;
END_RCPP
}
// var_lhood
double var_lhood(arma::mat Y, arma::vec par, int lags, int print_level, bool cond);
RcppExport SEXP VARext_var_lhood(SEXP YSEXP, SEXP parSEXP, SEXP lagsSEXP, SEXP print_levelSEXP, SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lhood(Y, par, lags, print_level, cond));
    return rcpp_result_gen;
END_RCPP
}
// var_lr_mu
arma::vec var_lr_mu(arma::vec a, arma::mat A, int lags);
RcppExport SEXP VARext_var_lr_mu(SEXP aSEXP, SEXP ASEXP, SEXP lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lr_mu(a, A, lags));
    return rcpp_result_gen;
END_RCPP
}
// var_lr_mu_grad
arma::mat var_lr_mu_grad(arma::vec a, arma::mat A, int lags);
RcppExport SEXP VARext_var_lr_mu_grad(SEXP aSEXP, SEXP ASEXP, SEXP lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lr_mu_grad(a, A, lags));
    return rcpp_result_gen;
END_RCPP
}
// var_lr_variance
arma::mat var_lr_variance(arma::mat A, arma::mat Sigma, int lags, bool contemp);
RcppExport SEXP VARext_var_lr_variance(SEXP ASEXP, SEXP SigmaSEXP, SEXP lagsSEXP, SEXP contempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< bool >::type contemp(contempSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lr_variance(A, Sigma, lags, contemp));
    return rcpp_result_gen;
END_RCPP
}
// var_init_lhood_grad_numeric
arma::vec var_init_lhood_grad_numeric(arma::vec Y0, arma::vec par, int n_var, int lags);
RcppExport SEXP VARext_var_init_lhood_grad_numeric(SEXP Y0SEXP, SEXP parSEXP, SEXP n_varSEXP, SEXP lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type n_var(n_varSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(var_init_lhood_grad_numeric(Y0, par, n_var, lags));
    return rcpp_result_gen;
END_RCPP
}
// var_lhood_grad
arma::vec var_lhood_grad(arma::mat Y, arma::vec par, int lags, int print_level, bool cond);
RcppExport SEXP VARext_var_lhood_grad(SEXP YSEXP, SEXP parSEXP, SEXP lagsSEXP, SEXP print_levelSEXP, SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lhood_grad(Y, par, lags, print_level, cond));
    return rcpp_result_gen;
END_RCPP
}
// var_lhood_grad_vcv
arma::mat var_lhood_grad_vcv(arma::mat Y, arma::vec par, int lags, int print_level, bool cond);
RcppExport SEXP VARext_var_lhood_grad_vcv(SEXP YSEXP, SEXP parSEXP, SEXP lagsSEXP, SEXP print_levelSEXP, SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lhood_grad_vcv(Y, par, lags, print_level, cond));
    return rcpp_result_gen;
END_RCPP
}
// var_lhood_grad_numeric
arma::vec var_lhood_grad_numeric(arma::mat Y, arma::vec par, int lags, int print_level, bool cond);
RcppExport SEXP VARext_var_lhood_grad_numeric(SEXP YSEXP, SEXP parSEXP, SEXP lagsSEXP, SEXP print_levelSEXP, SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lhood_grad_numeric(Y, par, lags, print_level, cond));
    return rcpp_result_gen;
END_RCPP
}
// var_lhood_hessian_numeric
arma::mat var_lhood_hessian_numeric(arma::mat Y, arma::vec par, int lags, int print_level, bool cond);
RcppExport SEXP VARext_var_lhood_hessian_numeric(SEXP YSEXP, SEXP parSEXP, SEXP lagsSEXP, SEXP print_levelSEXP, SEXP condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< int >::type print_level(print_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type cond(condSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lhood_hessian_numeric(Y, par, lags, print_level, cond));
    return rcpp_result_gen;
END_RCPP
}
// var_sim_e
arma::mat var_sim_e(arma::vec a, arma::mat A, arma::mat y0, arma::mat e);
RcppExport SEXP VARext_var_sim_e(SEXP aSEXP, SEXP ASEXP, SEXP y0SEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(var_sim_e(a, A, y0, e));
    return rcpp_result_gen;
END_RCPP
}
// var_e
arma::mat var_e(arma::mat Sigma, int n_pds);
RcppExport SEXP VARext_var_e(SEXP SigmaSEXP, SEXP n_pdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pds(n_pdsSEXP);
    rcpp_result_gen = Rcpp::wrap(var_e(Sigma, n_pds));
    return rcpp_result_gen;
END_RCPP
}
// var_sim
arma::mat var_sim(arma::vec a, arma::mat A, arma::mat y0, arma::mat Sigma, int n_pds);
RcppExport SEXP VARext_var_sim(SEXP aSEXP, SEXP ASEXP, SEXP y0SEXP, SEXP SigmaSEXP, SEXP n_pdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pds(n_pdsSEXP);
    rcpp_result_gen = Rcpp::wrap(var_sim(a, A, y0, Sigma, n_pds));
    return rcpp_result_gen;
END_RCPP
}
// markov_sim
arma::vec markov_sim(const int n, const NumericMatrix M, const int s0, const int n_s);
RcppExport SEXP VARext_markov_sim(SEXP nSEXP, SEXP MSEXP, SEXP s0SEXP, SEXP n_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< const int >::type n_s(n_sSEXP);
    rcpp_result_gen = Rcpp::wrap(markov_sim(n, M, s0, n_s));
    return rcpp_result_gen;
END_RCPP
}
