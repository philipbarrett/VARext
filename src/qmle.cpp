/***********************************************************************************
 * qmle.cpp
 *
 * Computes quasi-maximum likelihood function for a regime-switching model
 *
 * 14apr2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

/**
 * Do do:
 *  - Debug
 */

#include "qmle.hpp"
#include "lhood.hpp"

// // [[Rcpp::export]]
// arma::vec one_step_var_lhood( arma::mat Y, arma::mat par, int print_level = 0 ){
// // Computes the vector of conditional one-step likelihoods for the restricted
// // VAR, across different states, where the restrictions are assumed to be the
// // same across all the states.
//
//   int n_states = par.n_cols ;
//   vec out = zeros(n_states) ;
//       // Initialize the output
//   for ( int i = 0 ; i < n_states ; i++ ){
//     out(i) = var_lhood( Y, par.col(i), print_level ) ;
//   }
//   return out ;
// }
//
// // [[Rcpp::export]]
// double msw_var_lhood( arma::mat Y, arma::vec p0, arma::mat P,
//                            arma::mat par, int print_level = 0 ){
// // Computes the quasi-likelihood for the Markov switching model
//   double l_inc = 0 ;
//   double l = 0 ;
//       // Initialize the output
//   int n = p0.n_elem ;
//   int m = Y.n_cols ;
//       // Problem dimensions
//   vec l_vec = zeros(n) ;
//   vec p_vec = zeros(n) ;
//       // Vector of likelihoods and minus-log-likelihoods
//   mat P_t = P.t() ;
//       // Transpose of P
//   vec p_filter = p0 ;
//   vec p_pred = P_t * p_filter ;
//       // Initialize the prediction and filtering probabilities
//   for( int i = 1 ; i < m ; i ++ ){
//     l_vec = one_step_var_lhood( Y.cols(i-1,i), par, print_level ) ;
//     p_vec = exp( - l_vec ) ;
//         // The vector of likelihoods
//     l_inc = dot(p_pred, p_vec) ;
//     l += - std::log( l_inc ) / ( m - 1 ) ;
//         // Increment the likelihood
//     p_filter = p_pred % p_vec / l_inc ;
//     p_pred = P_t * p_filter ;
//         // Update the filtering and prediction probabilities
//   }
//   return l ;
// }
//
// // [[Rcpp::export]]
// arma::mat msw_var_lhood_p( arma::mat Y, arma::vec p0, arma::mat P,
//                                  arma::mat par, int print_level = 0 ){
// // Returns the forecast and filterin probabilities for the quasi-likelihood for
// // the Markov switching model
//   double l_inc = 0 ;
//       // Initialize the output
//   int n = p0.n_elem ;
//   int m = Y.n_cols ;
//       // Problem dimensions
//   vec l_vec = zeros(n) ;
//   vec p_vec = zeros(n) ;
//       // Vector of minus log-likelihoods and likelihoods
//   mat P_t = P.t() ;
//       // Transpose of P
//   vec p_filter = p0 ;
//   vec p_pred = P_t * p_filter ;
//       // Initialize the prediction and filtering probabilities
//   mat m_p_filter = zeros(n,m) ;
//   mat m_p_pred = zeros(n,m) ;
//       // The output matrices
//   m_p_filter.col(0) = p_filter ;
//   m_p_pred.col(0) = p_pred ;
//       // Initialization
//   for( int i = 1 ; i < m ; i ++ ){
//     l_vec = one_step_var_lhood( Y.cols(i-1,i), par, print_level ) ;
//     p_vec = exp( - l_vec ) ;
//         // The vector of likelihoods
//     l_inc = dot( p_pred, p_vec ) ;
//         // Increment the likelihood
//     p_filter = p_pred % p_vec / l_inc ;
//     p_pred = P_t * p_filter ;
//         // Update the filtering and prediction probabilities.  Transpose P to
//         // get a column vector output.
//     m_p_filter.col(i) = p_filter ;
//     m_p_pred.col(i) = p_pred ;
//   }
//   return join_cols( m_p_filter, m_p_pred ) ;
// }



