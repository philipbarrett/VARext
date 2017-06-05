/***********************************************************************************
 * disc.cpp
 *
 * Provides error functions and starting guesses for VAR discretization
 *
 * 03jun2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "disc.hpp"

// [[Rcpp::export]]
arma::vec rou_y( int N, double sd_u ){
// Calculates the Rouwenhorst (1995) grid points for a process with
// unconditional std dev sd_u
  vec out = zeros(N) ;
  for( int i = 0 ; i < N ; i ++ ){
    out(i) = - sd_u * std::sqrt( N - 1 ) + 2 * sd_u / std::sqrt( N - 1 ) *  i ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::vec M_i_ig( arma::mat Y, int i, arma::vec a, arma::mat A, arma::mat Sigma,
                  int print_level = 0 ){
// Creates an initial guess for M_i by matching the conditional mean, covariance and skew
  int n_pts = Y.n_rows ;                                    // Number of grid points
  int n_vars = Y.n_cols ;                                   // Number of variables
  vec y_i = conv_to<vec>::from( Y.row(i) ) ;                // Row i of the grid of points
  vec mu_i = a + A * y_i ;                                  // The conditional mean
  mat Y_dev = Y - ones( n_pts ) * conv_to<rowvec>::from(mu_i) ;
      // The matrix of deviations from the conditional mean
  mat Q = Y_dev % Y_dev % Y_dev ;                           // The conditional skews
  vec v_sig = zeros( .5 * n_vars * ( n_vars + 1 ) ) ;       // The vector of sigma
  mat Z = zeros( n_pts, .5 * n_vars * ( n_vars + 1 ) ) ;   // The conditional variances

  int counter = 0 ;
  for( int k = 0 ; k < n_vars ; k++ ){
    for( int j = k ; j < n_vars ; j++ ){
      v_sig(counter) = Sigma(k,j) ;
      Z.col(counter) = Y_dev.col(k) % Y_dev.col(j) ;
      counter ++ ;
    }
  }
      // Fill in the variance stuff

  mat B = zeros( n_vars + 1 + .5 * n_vars * ( n_vars + 1 ) + n_vars, n_pts ) ;
  B.rows(0,n_vars-1) = Y.t() ;
  B.row(n_vars) = ones<rowvec>(n_pts) ;
  B.rows( n_vars + 1, n_vars + 1 + .5 * n_vars * ( n_vars + 1 ) - 1 ) = Z.t() ;
  B.rows( n_vars + 1 + .5 * n_vars * ( n_vars + 1 ),
          2 * n_vars + .5 * n_vars * ( n_vars + 1 ) ) = Q.t() ;
  vec c = zeros( n_vars + 1 + .5 * n_vars * ( n_vars + 1 ) + n_vars ) ;
  c.subvec(0,n_vars-1) = mu_i ;
  c(n_vars) = 1 ;
  c.subvec( n_vars + 1, n_vars + .5 * n_vars * ( n_vars + 1 ) ) = v_sig ;
      // Create the matrices for B.m=c

  vec m_i = solve( B, c ) ;
      // Solve the equations

  if(print_level > 0 ){
    Rcout << "B:\n" << B << std::endl ;
    Rcout << "c:\n" << c << std::endl ;
    Rcout << "m_i:\n" << m_i << std::endl ;
  }

  return( m_i ) ;
}

// [[Rcpp::export]]
arma::mat var_cm( arma::mat Y, arma::vec a, arma::mat A ){
// Computes the conditional mean fo the VAR at the grid of points Y
  mat cm = a * ones<rowvec>( Y.n_rows )  + A * Y.t() ;
  return cm.t() ;
}

// [[Rcpp::export]]
arma::mat disc_cv( arma::mat Y, arma::mat M ){
// Computes the conditional variance for the Markov process with nodes Y and
// transition matrix M
  mat cm = M * Y ;                                        // The conditional mean
  int n_var = Y.n_cols ;                                  // The number of variables
  int n_pts = Y.n_rows ;                                  // The number of nodes
  mat out = zeros( .5 * n_var * ( n_var + 1 ), n_pts ) ;  // Matrix of covariances

  for ( int i_pt = 0 ; i_pt < n_pts ; i_pt++ ){
    int counter = 0 ;
    for( int i = 0 ; i < n_var ; i++ ){
      for( int j = i ; j < n_var ; j++ ){
        out(counter,i_pt) = sum( M.row(i_pt) *
                      ( ( Y.col(i) - cm(i_pt,i) ) % ( Y.col(j) - cm(i_pt,j) ) ) ) ;
        counter++ ;
      }
    }
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat disc_skew( arma::mat Y, arma::mat M ){
// Computes the conditional variance for the Markov process with nodes Y and
// transition matrix M
  mat cm = M * Y ;                                        // The conditional mean
  int n_var = Y.n_cols ;                                  // The number of variables
  int n_pts = Y.n_rows ;                                  // The number of nodes
  mat out = zeros( n_var, n_pts ) ;                       // Matrix of skews
  for ( int i_pt = 0 ; i_pt < n_pts ; i_pt++ ){
    for( int i = 0 ; i < n_var ; i++ ){
      vec Y_dev = Y.col(i) - cm(i_pt,i) ;
      out(i,i_pt) = sum( M.row(i_pt) * ( Y_dev % Y_dev % Y_dev ) ) ;
    }
  }
  return out ;
}

// [[Rcpp::export]]
double disc_obj( arma::vec m_i, arma::mat Y, int i, arma::vec a, arma::mat A,
                 arma::mat Sigma, arma::vec w, int print_level=0 ){
// The objective function for the discretized VAR
  int n_pts = Y.n_rows ;                                    // Number of grid points
  int n_vars = Y.n_cols ;                                   // Number of variables
  vec y_i = conv_to<vec>::from( Y.row(i) ) ;                // Row i of the grid of points
  vec mu_i = a + A * y_i ;                                  // The conditional mean
  mat Y_dev = Y - ones( n_pts ) * conv_to<rowvec>::from(mu_i) ;
                                  // The matrix of deviations from the conditional mean
  mat Q = Y_dev % Y_dev % Y_dev ;                           // The conditional skews
  vec v_sig = zeros( .5 * n_vars * ( n_vars + 1 ) ) ;       // The vector of sigma
  mat Z = zeros( n_pts, .5 * n_vars * ( n_vars + 1 ) ) ;    // The conditional variances
  mat W = diagmat(w) ;                                      // The vector of weights

  int counter = 0 ;
  for( int k = 0 ; k < n_vars ; k++ ){
    for( int j = k ; j < n_vars ; j++ ){
      v_sig(counter) = Sigma(k,j) ;
      Z.col(counter) = Y_dev.col(k) % Y_dev.col(j) ;
      counter ++ ;
    }
  }
  // Fill in the variance stuff

  mat B = zeros( n_vars + .5 * n_vars * ( n_vars + 1 ) + n_vars, n_pts ) ;
  vec c = zeros( n_vars + .5 * n_vars * ( n_vars + 1 ) + n_vars ) ;
      // Initialize the objective matrix
  B.rows(0,n_vars-1) = Y.t() ;
  B.rows( n_vars , n_vars + .5 * n_vars * ( n_vars + 1 ) - 1 ) = Z.t() ;
  B.rows( n_vars  + .5 * n_vars * ( n_vars + 1 ),
          2 * n_vars + .5 * n_vars * ( n_vars + 1 ) - 1 ) = Q.t() ;
  c.subvec(0,n_vars-1) = mu_i ;
  c.subvec( n_vars, n_vars + .5 * n_vars * ( n_vars + 1 ) - 1 ) = v_sig ;
      // Create the matrices for eps = B.m - c

      if(print_level > 0 ){
        Rcout << "B:\n" << B << std::endl ;
        Rcout << "c:\n" << c << std::endl ;
        Rcout << "m_i:\n" << m_i << std::endl ;
      }

  vec eps = B * m_i - c ;
  double out = .5 * sum( eps.t() * W * eps ) ;
  return out ;
}


// [[Rcpp::export]]
arma::vec disc_grad( arma::vec m_i, arma::mat Y, int i, arma::vec a, arma::mat A,
                 arma::mat Sigma, arma::vec w ){
  // The objective function for the discretized VAR
  int n_pts = Y.n_rows ;                                    // Number of grid points
  int n_vars = Y.n_cols ;                                   // Number of variables
  vec y_i = conv_to<vec>::from( Y.row(i) ) ;                // Row i of the grid of points
  vec mu_i = a + A * y_i ;                                  // The conditional mean
  mat Y_dev = Y - ones( n_pts ) * conv_to<rowvec>::from(mu_i) ;
  // The matrix of deviations from the conditional mean
  mat Q = Y_dev % Y_dev % Y_dev ;                           // The conditional skews
  vec v_sig = zeros( .5 * n_vars * ( n_vars + 1 ) ) ;       // The vector of sigma
  mat Z = zeros( n_pts, .5 * n_vars * ( n_vars + 1 ) ) ;    // The conditional variances
  mat W = diagmat(w) ;                                      // The vector of weights

  int counter = 0 ;
  for( int k = 0 ; k < n_vars ; k++ ){
    for( int j = k ; j < n_vars ; j++ ){
      v_sig(counter) = Sigma(k,j) ;
      Z.col(counter) = Y_dev.col(k) % Y_dev.col(j) ;
      counter ++ ;
    }
  }
  // Fill in the variance stuff

  mat B = zeros( n_vars + .5 * n_vars * ( n_vars + 1 ) + n_vars, n_pts ) ;
  vec c = zeros( n_vars + .5 * n_vars * ( n_vars + 1 ) + n_vars ) ;
      // Initialize the objective matrix
  B.rows(0,n_vars-1) = Y.t() ;
  B.rows( n_vars , n_vars + .5 * n_vars * ( n_vars + 1 ) - 1 ) = Z.t() ;
  B.rows( n_vars  + .5 * n_vars * ( n_vars + 1 ),
          2 * n_vars + .5 * n_vars * ( n_vars + 1 ) - 1 ) = Q.t() ;
  c.subvec(0,n_vars-1) = mu_i ;
  c.subvec( n_vars, n_vars + .5 * n_vars * ( n_vars + 1 ) - 1 ) = v_sig ;
      // Create the matrices for eps = B.m - c
  vec eps = B * m_i - c ;
  return B.t() * ( W * eps ) ;
}
