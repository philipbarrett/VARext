/***********************************************************************************
 * sim.hpp
 *
 * Interface to sim.cpp
 *
 * 17may2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#ifndef SIM_HPP
#define SIM_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat var_sim_e( arma::vec a, arma::mat A, arma::mat y0, arma::mat e ) ;
arma::mat var_e( arma::mat Sigma, int n_pds ) ;
arma::mat var_sim( arma::vec a, arma::mat A, arma::mat y0, arma::mat Sigma, int n_pds ) ;
arma::vec markov_sim( const int n, const NumericMatrix M, const int s0, const int n_s ) ;

#endif
