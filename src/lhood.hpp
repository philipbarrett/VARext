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

double var_lhood( arma::mat Y, arma::vec par, int lags, int print_level ) ;
arma::vec var_lhood_grad_numeric( arma::mat Y, arma::vec par, int lags, int print_level ) ;

#endif
