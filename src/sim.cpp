/***********************************************************************************
 * sim.cpp
 *
 * Simulates a VAR
 *
 * 17may2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "sim.hpp"

// [[Rcpp::export]]
arma::mat var_sim_e( arma::vec a, arma::mat A, arma::mat y0, arma::mat e ){
// Simulates a VAR given a matrix of errors

  int lags = y0.n_cols ;                    // Number of lags
  int n_pds = e.n_cols ;                    // Number of simulation periods
  int m = a.n_elem ;                        // Number of variables
  mat sim = zeros( m, n_pds + lags ) ;      // Initialize the simulation matrix
  sim.cols( 0, lags - 1 ) = y0 ;            // Put in the simulated columns

  uvec v_rev = zeros<uvec>( lags * m ) ;
  for ( int i = 0 ; i < lags ; i++ ){
    for ( int j = 0 ; j < m ; j++ ){
      v_rev(i*m+j) = m*(lags-i-1) + j;
    }
  }
  mat A_rev = A.cols( v_rev ) ;
      // Need to reorder the columns of A (so clumsy!)

  for( int i = 0 ; i < n_pds ; i ++ ){
    sim.col(lags+i) = a + A_rev * vectorise( sim.cols(i,i+lags-1) ) + e.col(i) ;
        // Update the simulation
  }
  return sim ;
}

// [[Rcpp::export]]
arma::mat var_e( arma::mat Sigma, int n_pds ){
// Draws errors with variance-covariance matrix Sigma

  int m = Sigma.n_rows ;                    // Number of variables
  mat sd = chol(Sigma).t() ;                // Square root of sigma: sd * sd.t() = Sigma
  mat err_0 = randn( m, n_pds ) ;           // Standard normal simulations
  mat err = sd * err_0 ;                    // The matrix of errors with variance sigma
  return err ;
}

// [[Rcpp::export]]
arma::mat var_sim( arma::vec a, arma::mat A, arma::mat y0, arma::mat Sigma, int n_pds ){
// Simulates a VAR

  mat e = var_e( Sigma, n_pds ) ;           // The errors
  return var_sim_e( a, A, y0, e ) ;
}
