/***********************************************************************************
 * lhood.cpp
 *
 * Computes maximum likelihood function for a VAR
 *
 * 17may2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "lhood.hpp"

// [[Rcpp::export]]
double var_lhood( arma::mat Y, arma::vec par, int lags = 1, int print_level = 0 ){
// Computes the (negative) likelihood of the restricted VAR of the form:
//    Y_t = a + sum_{i=1}^lags A_i * Y_{t-i} + e_t,  e_t ~ N(0,Sigma)
// Where the integer matrices A_switch and Sigma_switch define the
// parameter restrictions (zeros for restricted to zero, one otherwise)

  int m = Y.n_cols ;        // Number of periods
  int n = Y.n_rows ;        // Number of variables

  mat Yt = Y.cols(lags,m-1) ;
  mat Yt_lag = zeros( lags * n, m - lags ) ;
  for( int i = 0 ; i < lags ; i ++ ){
    Yt_lag.rows( n*i, n*(i+1)-1 ) = Y.cols( lags-1-i, m-2-i ) ;
  }
      // Create the lagged variables

  vec a = zeros(n) ;
  mat A = zeros(n,n*lags) ;
  mat Sigma = zeros(n,n) ;
      // Initialize the matrices of parameters
  int counter = 0 ;
      // Counter
  for( int i = 0 ; i < n ; i ++ ){
    a(i) = par(counter) ;
    counter++ ;
  }
  for( int j = 0 ; j < n*lags ; j ++ ){
    for( int i = 0 ; i < n ; i ++ ){
      A(i,j) = par(counter) ;
      counter++ ;
    }
  }
  // Copy parameters into A

  for( int i = 0 ; i < n ; i ++ ){
    for( int j = 0 ; j <= i ; j ++ ){
      Sigma(i,j) = par(counter) ;
      Sigma(j,i) = par(counter) ;
      counter++ ;
    }
  }
  // Copy parameters into Sigma
  mat e = Yt - A * Yt_lag - a * ones<rowvec>(m-lags) ;
      // The matrix of residuals
  mat Sigma_I = Sigma.i() ;
      // Inverse of Sigma
  mat e_t = e.t() ;
      // Transpose of the error terms
  double term = 0 ;
      // The final term in the likelihood
  for( int i = 0 ; i < m - lags ; i ++ ){
    vec temp = e_t.row(i) * Sigma_I * e.col(i) ;
    term += temp(0) ;
  }
  if( print_level > 0 ){
    Rcout << "a:\n" << a <<std::endl ;
    Rcout << "A:\n" << A <<std::endl ;
    Rcout << "Sigma:\n" << Sigma <<std::endl ;
    Rcout << "Sigma_I:\n" << Sigma_I <<std::endl ;
    Rcout << "term = " << term <<std::endl ;
    Rcout << "M_PI = " << M_PI <<std::endl ;
    Rcout << "std::log( 2 * M_PI ) = " << std::log( 2 * M_PI ) <<std::endl ;
    Rcout << "det(Sigma) = " << det(Sigma) <<std::endl ;
    Rcout << "std::log(det(Sigma)) = " << std::log(det(Sigma)) <<std::endl ;
    Rcout << "Yt.cols(0,4):\n" << Yt.cols(0,4) << std::endl ;
    Rcout << "Yt_lag.cols(0,4):\n" << Yt_lag.cols(0,4) << std::endl ;
    Rcout << "e.cols(0,4):\n" << e.cols(0,4) << std::endl ;
  }

  double l = .5 * ( n * std::log( 2 * M_PI ) + std::log( det( Sigma ) ) )
    + .5 * term / ( m - lags ) ;
      // The negative log likelihood
  return l ;
}

// [[Rcpp::export]]
arma::vec var_lhood_grad_numeric( arma::mat Y, arma::vec par, int lags = 1, int print_level = 0 ){
// Computes the likelihood of the likelihood function numerically

  int n_par = par.n_elem ;                                  // The number of parameters
  vec out = zeros( n_par ) ;                                // The output vector
  vec par_2 = par ;                                         // The perturbed parameter vector
  double inc = 1e-06 ;                                      // Increment
  double l = var_lhood( Y, par, lags, print_level ) ;       // The base likelihood

  for( int i = 0 ; i < n_par ; i++ ){
    par_2 = par ;
    par_2(i) += inc ;
        // Create the perturbed vector
    double l_2 = var_lhood( Y, par_2, lags, print_level ) ;
        // The perturbed likelihood
    out(i) = ( l_2 - l ) / inc ;
  }
  return out ;
}
