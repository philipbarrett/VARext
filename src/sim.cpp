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

// [[Rcpp::export]]
arma::vec markov_sim( const int n, const NumericMatrix M, const int s0,
                      const int n_s ){
  // Fast Markov simulation

  NumericVector out(n) ;
  out(0) = s0 ;
  // Initialize output
  NumericMatrix sum_M( n_s, n_s ) ;
  // The row sum of M
  for( int i=0 ; i < n_s ; i++ ){
    sum_M( i, 0 ) = M( i , 0 ) ;
    // The first column
    for( int j=1 ; j < n_s ; j++ ){
      sum_M( i, j ) = sum_M( i, j - 1 ) + M( i , j ) ;
    }
  } // Create the row sum of M

  NumericVector shks = runif( n ) ;
  // The vector of random shocks on [0,1]
  int this_s = 0 ;

  for( int i = 1 ; i < n ; i++ ){
    this_s = 0 ;
    // Initialize counters
    while( shks(i) > sum_M( out(i-1), this_s ) ){
      this_s++ ;
    }
    out( i ) = this_s ;
  }

  vec out_arma = zeros(n) ;
  for ( int i = 0 ; i < n ; i++ )
    out_arma(i) = out(i) ;
  // Convert to arma type

  return out_arma ;
}
