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
double var_lhood_N( arma::mat Y, arma::vec par, int lags = 1, int print_level = 0,
                    bool cond = true ){
// Computes the unscaled (negative) likelihood of the restricted VAR of the form:
//    Y_t = a + sum_{i=1}^lags A_i * Y_{t-i} + e_t,  e_t ~ N(0,Sigma)
// The boolean switch "cond" determines if the likelihood is conditional on the
// first observation

  int m = Y.n_cols ;        // Number of periods
  int n = Y.n_rows ;        // Number of variables
  return var_lhood( Y, par, lags, print_level, cond ) * ( m - lags ) ;
}

// [[Rccp::export]]
double var_init_lhood( arma::vec Y_0, arma::vec a, arma::mat A, arma::mat Sigma, int lags ){
// Compute the likelihood contribution from the initial point
  int n = A.n_rows ;                                          // Number of varirables
  mat lr_var = var_lr_variance( A, Sigma, lags, true ) ;      // Long run variance: one term only
  vec mu = var_lr_mu( a, A, lags ) ;                          // The unconditional mean
  vec temp = ( Y_0 - mu ).t() * inv( lr_var ) * ( Y_0 - mu ) ;
      // The quadratic term in the likelihood - first period only.
  return .5 * ( n * std::log( 2 * M_PI ) + std::log( det( lr_var ) )
                  + temp(0) ) ;                               // Return the likelihood
}

// [[Rccp::export]]
double var_init_lhood_par( arma::vec Y_0, arma::vec par, int n_var, int lags ){
// par form

  /** Create base likelihood **/
  vec a = zeros(n_var) ;
  mat A = zeros(n_var,n_var*lags) ;
  mat Sigma = zeros(n_var,n_var) ;
      // Initialize the matrices of parameters
  int counter = 0 ;
      // Counter
  for( int i = 0 ; i < n_var ; i ++ ){
    a(i) = par(counter) ;
    counter++ ;
  }
  for( int j = 0 ; j < n_var*lags ; j ++ ){
    for( int i = 0 ; i < n_var ; i ++ ){
      A(i,j) = par(counter) ;
      counter++ ;
    }
  }
      // Copy parameters into A

  for( int i = 0 ; i < n_var ; i ++ ){
    for( int j = 0 ; j <= i ; j ++ ){
      Sigma(i,j) = par(counter) ;
      Sigma(j,i) = par(counter) ;
      counter++ ;
    }
  }
      // Copy into Sigma
  return var_init_lhood( Y_0, a, A, Sigma, lags ) ;
}

// [[Rcpp::export]]
double var_lhood( arma::mat Y, arma::vec par, int lags = 1, int print_level = 0,
                  bool cond = true ){
// Computes the (negative) likelihood of the restricted VAR of the form:
//    Y_t = a + sum_{i=1}^lags A_i * Y_{t-i} + e_t,  e_t ~ N(0,Sigma)
// The boolean switch "cond" determines if the likelihood is conditional on the
// first observation

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
  double l = .5 * ( n * std::log( 2 * M_PI ) + std::log( det( Sigma ) ) )
    + .5 * term / ( m - lags ) ;
      // The negative log likelihood

  if( print_level > 0 ){
    Rcout << "n = " << n <<std::endl ;
    Rcout << "m = " << m <<std::endl ;
    Rcout << "a:\n" << a <<std::endl ;
    Rcout << "A:\n" << A <<std::endl ;
    Rcout << "Sigma:\n" << Sigma <<std::endl ;
    Rcout << "Sigma_I:\n" << Sigma_I <<std::endl ;
    Rcout << "term = " << term <<std::endl ;
    Rcout << "M_PI = " << M_PI <<std::endl ;
    Rcout << "std::log( 2 * M_PI ) = " << std::log( 2 * M_PI ) <<std::endl ;
    Rcout << "det(Sigma) = " << det(Sigma) <<std::endl ;
    Rcout << "std::log(det(Sigma)) = " << std::log(det(Sigma)) <<std::endl ;
    Rcout << "l = " << l <<std::endl ;
    Rcout << "Yt.cols(0,4):\n" << Yt.cols(0,4) << std::endl ;
    Rcout << "Yt_lag.cols(0,4):\n" << Yt_lag.cols(0,4) << std::endl ;
    Rcout << "e.cols(0,4):\n" << e.cols(0,4) << std::endl ;
  }

  if( !cond ){
    vec Y_0 = Y.col(lags - 1 ) ;
        // Initial term(s)
    l += var_init_lhood( Y_0, a, A, Sigma, lags ) / ( m - lags ) ;
        // Augment the likelihood
    if( print_level > 0 ){
      // Rcout << "lr_var:\n" << lr_var <<std::endl ;
      Rcout << "Y_0:\n" << Y_0 <<std::endl ;
      // Rcout << "temp(0) = " << temp(0) <<std::endl ;
      Rcout << "l = " << l <<std::endl ;
      // Rcout << "det(lr_var) = " << det(lr_var) <<std::endl ;
      // Rcout << "std::log(det(lr_var)) = " << det(lr_var) <<std::endl ;
    }
  }
  return l ;
}

// [[Rcpp::export]]
arma::vec var_lr_mu( arma::vec a, arma::mat A, int lags ){
// Computes the long run mean
  int n_var = a.n_elem ;
      // Number of variables
  mat phi = zeros( n_var, n_var ) ;
      // Matrix for inversion
  for( int i = 0 ; i < lags ; i ++ ){
    phi = phi + A.cols( i*n_var, (i+1)*n_var-1 ) ;
  }
  return( solve( eye(n_var,n_var) - phi, a ) ) ;
}

// [[Rcpp::export]]
arma::mat var_lr_mu_grad( arma::vec a, arma::mat A, int lags ){
// Computes the derivative of the long run mean
  int n_var = a.n_elem ;
      // Number of variables
  mat phi = zeros( n_var, n_var ) ;
      // Matrix for inversion
  for( int i = 0 ; i < lags ; i ++ ){
    phi = phi + A.cols( i*n_var, (i+1)*n_var-1 ) ;
  }
  mat i_phi = inv( eye(n_var,n_var) - phi ) ;
      // The inverse
  mat out = zeros( n_var * ( 1 + lags * n_var ), n_var ) ;
      // The output
  out.rows(0,n_var-1) = i_phi.t() ;
      // Derivatives wrt a
  mat d_phi = zeros( n_var, n_var * n_var ) ;
      // The derivatives wrt A
  for( int i = 0 ; i < phi.n_elem ; i++ ){
    mat temp = zeros( size( phi ) ) ;
    temp(i) = 1 ;
        // Temporary vector to extract the ith coefficient of the derivative
    d_phi.col(i) = i_phi * temp * i_phi * a ;
  }

  for( int i = 0 ; i < lags ; i ++ ){
    out.rows(n_var+i*n_var*n_var,n_var+(i+1)*n_var*n_var-1) = d_phi.t() ;
  }
  return( out ) ;
}

// [[Rcpp::export]]
arma::mat var_lr_variance( arma::mat A, arma::mat Sigma, int lags=1, bool contemp=true ){
// Computes the long-run variance of the VAR
  int n_var = A.n_rows ;                      // Number of variables
  mat phi = zeros(n_var*lags,n_var*lags) ;    // The relvant auto-regressive matrix
  mat S = zeros(n_var*lags,n_var*lags) ;      // Relevant variance (needs cahnging for lags>1)

  if( lags == 1 ){
    phi = A ;
    S = Sigma ;
  }else{
    phi.rows(0,n_var-1) = A ;
    for( int i = 1 ; i < lags ; i ++ ){
      phi.submat(i*n_var,(i-1)*n_var,(i+1)*n_var-1,i*n_var-1) = eye(n_var,n_var) ;
    }
    S.submat(0,0,n_var-1,n_var-1) = Sigma ;
        // The companion form
  }
  mat out = mat( n_var*lags,n_var*lags ) ;
      // The long-run variance
  vec v_out( out.memptr(), out.n_elem, false, false);
      // Share memory with the vectorized form
  int n_v_out = std::pow(n_var*lags,2) ;
      // The number of elements of the vectorized output
  v_out = solve( eye( n_v_out, n_v_out ) - kron( phi, phi ),
                 vectorise( S ) ) ;
      // Calculate the vectorized output

  if( contemp )
    return( out.submat(0,0,n_var-1,n_var-1) ) ;
  return out ;
}

// [[Rcpp::export]]
arma::vec var_init_lhood_grad_numeric( arma::vec Y0, arma::vec par, int n_var, int lags=1 ){
// Computes the derivative of the likelihood initial point

  vec out = 0 * par ;                 // Initialize output
  int n_par = par.n_elem ;            // Number of parameters
  vec par_2 = par ;                   // Copy to incremented version
  double inc = 1e-6 ;                 // Increment
  double l_base = var_init_lhood_par( Y0, par, n_var, lags ) ;

  /** Create new likelihood **/
  for( int k = 0 ; k < n_par ; k ++ ){
    par_2 = par ;
    par_2(k) += inc ;
        // Increment
    double l_2 = var_init_lhood_par( Y0, par_2, n_var, lags ) ;
        // The perturbed likelihood
    out(k) = ( l_2 - l_base ) / inc ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::vec var_lhood_grad( arma::mat Y, arma::vec par, int lags = 1, int print_level = 0,
                          bool cond = true ){
// Computes the derivative of the (negative) likelihood of the restricted VAR of the form:
//    Y_t = a + sum_{i=1}^lags A_i * Y_{t-i} + e_t,  e_t ~ N(0,Sigma)

  int m = Y.n_cols ;        // Number of periods
  int n = Y.n_rows ;        // Number of variables
  int n_par = par.n_elem ;  // The number of parameters

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
  vec deriv = zeros( n_par ) ;
      // The vector of derivatives the likelihood
  for( int i = 0 ; i < m - lags ; i ++ ){

    vec l_mu = Sigma_I * e.col(i) ;         // Deriv of l wrt the mean in pd i
    vec contrib_a = l_mu ;                  // Moving const terms affects means 1-for-1
    vec contrib_A = vectorise( (Yt_lag.col(i) * l_mu.t()).t() ) ;
        // Changing the A parameters affects the likelihood the mean by the
        // appropriate lagged variable
    mat m_contrib_sig = - .5 * ( Sigma_I - Sigma_I * ( e.col(i) * e.col(i).t() ) * Sigma_I ) ;
        // The covariance matrix derivative
    vec contrib_sig = zeros( .5 * n * ( 1 + n ) ) ;
    int counter = 0 ;
    for( int j = 0 ; j < n ; j ++ ){
      for( int k = 0 ; k <= j ; k ++ ){
        int mult = ( j == k ) ? 1 : 2 ;
        contrib_sig( counter ) = mult * m_contrib_sig(k,j) ;
        counter ++ ;
            // Repackage in the par form
      }
    }   // Map the lower triangle into a vector
    vec contrib = zeros(n_par) ;
    contrib.head(n) = contrib_a ;
    contrib.subvec( n, n*lags*n+n-1 ) = contrib_A ;
    contrib.subvec( n*lags*n+n, n_par - 1 ) = contrib_sig ;
    deriv += contrib ;
  }
  if( print_level > 0 ){
    Rcout << "a:\n" << a <<std::endl ;
    Rcout << "A:\n" << A <<std::endl ;
    Rcout << "Sigma:\n" << Sigma <<std::endl ;
    Rcout << "Sigma_I:\n" << Sigma_I <<std::endl ;
    Rcout << "deriv = " << deriv <<std::endl ;
    Rcout << "Yt.cols(0,4):\n" << Yt.cols(0,4) << std::endl ;
    Rcout << "Yt_lag.cols(0,4):\n" << Yt_lag.cols(0,4) << std::endl ;
    Rcout << "e.cols(0,4):\n" << e.cols(0,4) << std::endl ;
  }

  if( cond )
    return - deriv / ( m - lags ) ;
        // Conditional likelihood derivative
  return ( - deriv + var_init_lhood_grad_numeric( Y.col(lags-1), par, n, lags ) ) / ( m - lags ) ;
      // Have to use numeric derivative for the initial condition

}
// [[Rcpp::export]]
arma::mat var_lhood_grad_vcv( arma::mat Y, arma::vec par, int lags = 1,
                              int print_level = 0, bool cond = true ){
// Computes the variance of derivative of the (negative) likelihood of the
// restricted VAR of the form:
//    Y_t = a + sum_{i=1}^lags A_i * Y_{t-i} + e_t,  e_t ~ N(0,Sigma)

  if(!cond){
    Rcout << "Warning: var_lhood_grad_vcv does not yet allow for unconditional likelihood."
          << std::endl ;
  }

  int m = Y.n_cols ;        // Number of periods
  int n = Y.n_rows ;        // Number of variables
  int n_par = par.n_elem ;  // The number of parameters

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
  mat deriv2 = zeros( n_par, n_par ) ;
      // The vector of derivatives the likelihood
  for( int i = 0 ; i < m - lags ; i ++ ){

    vec l_mu = Sigma_I * e.col(i) ;         // Deriv of l wrt the mean in pd i
    vec contrib_a = l_mu ;                  // Moving const terms affects means 1-for-1
    vec contrib_A = vectorise( (Yt_lag.col(i) * l_mu.t()).t() ) ;
        // Changing the A parameters affects the likelihood the mean by the
        // appropriate lagged variable
    mat m_contrib_sig = - .5 * ( Sigma_I - Sigma_I * ( e.col(i) * e.col(i).t() ) * Sigma_I ) ;
        // The covariance matrix derivative
    vec contrib_sig = zeros( .5 * n * ( 1 + n ) ) ;
    int counter = 0 ;
    for( int j = 0 ; j < n ; j ++ ){
      for( int k = 0 ; k <= j ; k ++ ){
        int mult = ( j == k ) ? 1 : 2 ;
        contrib_sig( counter ) = mult * m_contrib_sig(k,j) ;
        counter ++ ;
            // Repackage in the par form
      }
    }   // Map the lower triangle into a vector
    vec contrib = zeros(n_par) ;
    contrib.head(n) = contrib_a ;
    contrib.subvec( n, n*lags*n+n-1 ) = contrib_A ;
    contrib.subvec( n*lags*n+n, n_par - 1 ) = contrib_sig ;
    deriv2 += (contrib * contrib.t()) ;
  }
  if( print_level > 0 ){
    Rcout << "a:\n" << a <<std::endl ;
    Rcout << "A:\n" << A <<std::endl ;
    Rcout << "Sigma:\n" << Sigma <<std::endl ;
    Rcout << "Sigma_I:\n" << Sigma_I <<std::endl ;
    Rcout << "deriv2:\n" << deriv2 <<std::endl ;
    Rcout << "Yt.cols(0,4):\n" << Yt.cols(0,4) << std::endl ;
    Rcout << "Yt_lag.cols(0,4):\n" << Yt_lag.cols(0,4) << std::endl ;
    Rcout << "e.cols(0,4):\n" << e.cols(0,4) << std::endl ;
  }

  return deriv2 / ( m - lags ) ;
}

// [[Rcpp::export]]
arma::vec var_lhood_grad_numeric( arma::mat Y, arma::vec par, int lags = 1, int print_level = 0,
                                  bool cond = true ){
// Computes the derivative of the likelihood function numerically

  int n_par = par.n_elem ;                                  // The number of parameters
  vec out = zeros( n_par ) ;                                // The output vector
  vec par_2 = par ;                                         // The perturbed parameter vector
  double inc = 1e-06 ;                                      // Increment
  double l = var_lhood( Y, par, lags, print_level, cond ) ; // The base likelihood

  for( int i = 0 ; i < n_par ; i++ ){
    par_2 = par ;
    par_2(i) += inc ;
        // Create the perturbed vector
    double l_2 = var_lhood( Y, par_2, lags, print_level, cond ) ;
        // The perturbed likelihood
    out(i) = ( l_2 - l ) / inc ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat var_lhood_hessian_numeric( arma::mat Y, arma::vec par, int lags = 1,
                                     int print_level = 0, bool cond=true ){
// Computes the second derivative of the likelihood function numerically

  int n_par = par.n_elem ;                                        // The number of parameters
  mat out = zeros( n_par, n_par ) ;                               // The output vector
  vec par_u = par ;                                               // The perturbed parameter vector
  vec par_d = par ;                                               // The perturbed parameter vector
  double inc = 1e-06 ;                                            // Increment
  // vec l_d = var_lhood_grad( Y, par, lags, print_level, cond ) ;   // The base gradient

  for( int i = 0 ; i < n_par ; i++ ){
    par_u = par ;
    par_u(i) += inc ;
    par_d = par ;
    par_d(i) -= inc ;
        // Create the perturbed vector
    vec l_u = var_lhood_grad( Y, par_u, lags, print_level, cond ) ;
    vec l_d = var_lhood_grad( Y, par_d, lags, print_level, cond ) ;
        // The perturbed likelihood
    out.col(i) = ( l_u - l_d ) / ( 2 * inc ) ;
  }
  return out ;
}

