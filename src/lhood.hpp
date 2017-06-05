/***********************************************************************************
 * lhood.hpp
 *
 * Interface to lhood.cpp
 *
 * 17may2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#ifndef LHOOD_HPP
#define LHOOD_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

double var_lhood( arma::mat Y, arma::vec par, int lags, int print_level, bool cond ) ;
arma::vec var_lhood_grad_numeric( arma::mat Y, arma::vec par, int lags,
                                  int print_level, bool cond ) ;
arma::mat var_lr_variance( arma::mat A, arma::mat Sigma, int lags, bool contemp ) ;
arma::vec var_lr_mu( arma::vec a, arma::mat A, int lags ) ;
arma::vec var_lhood_grad( arma::mat Y, arma::vec par, int lags, int print_level,
                          bool cond ) ;
double var_init_lhood( arma::vec Y_0, arma::vec a, arma::mat A, arma::mat Sigma, int lags ) ;
double var_init_lhood_par( arma::vec Y_0, arma::vec par, int n_var, int lags ) ;

#endif
