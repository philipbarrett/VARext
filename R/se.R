####################################################################################
# se.R
#
# Computes standard errors for the VAR
# 24may2017
# Philip Barrett, Washington DC
#
####################################################################################

var.ols.se <- function( Y, lags=1 ){
# Computes the standrad errors of the autoregression coefficients
  n.var <- nrow(Y)
  n.pds <- ncol(Y)
  if(is.null(rownames(Y))) rownames(Y) <- paste0('y', 1:n.var)
      # Add some row names to revent errors
  v.est <- VAR(t(Y), lags )
      # The VAR estimates
  x.data <- as.matrix(v.est$datamat[,-(1:n.var)])
  x.data <- x.data[ , c( ncol(x.data), 1:(ncol(x.data)-1) ) ]
  x.xprime <- t(x.data) %*% x.data
  x.xprime.i <- solve( x.xprime )
      # The data matrix and associated outputs
  Sigma <- var(sapply(v.est$varresult, function(x) x$residuals ) )
  ols.var <- kronecker( Sigma, x.xprime.i )
  # ols.var <- lapply( 1:n.var, function(i) Sigma[i,i] * x.xprime.i )
      # The "standard" variance, under independent, homoskedastic errors
  n.coeff <-  1 + n.var * lags
  ols.var.nw <- matrix( 0, n.var * n.coeff, n.var * n.coeff )
  for( i in 1:n.var ){
    for( j in 1:n.var ){
      res.var.nw <- cov( x.data * v.est$varresult[[i]]$residuals,
                         x.data * v.est$varresult[[j]]$residuals ) * ( n.pds - lags * n.var )
          # The NW residual variance
      ols.var.nw[ (i-1)*n.coeff + 1:n.coeff, (j-1)*n.coeff + 1:(1+n.var*lags) ] <-
                          x.xprime.i %*% res.var.nw %*% x.xprime.i
          # The Newey-West variance estimator. Corrects for heteroskedasticity
          # (which in this case means that the variance of the innovation is not
          # constant over time).
    }
  }
  # res.var.nw <- lapply(v.est$varresult, function(x) var( x.data * x$residuals ) ) * ( n.pds - lags * n.var )
  #     # The central variance term in the Newey-West estimator: Var( X'eps | X )
  # ols.var.nw <- lapply( res.var.nw, function(v) x.xprime.i %*% v %*% x.xprime.i )
  idx <- c( sapply( 1:n.coeff, function(i) i + (1:n.var-1) * n.coeff ) )
      # Reorder in the a, A, (Sigma) form that par and MLE estimates use
  return( list( ols=ols.var[idx,idx], ols.nw=ols.var.nw[idx,idx] ) )
}

var.mle.se <- function( Y, l.var ){
# Computes the variance of the estimators for the MLE from the numerical Hessian
# of the likelihood
  par <- par.to( l.var$a, l.var$A, l.var$Sigma )              # Convert to par format
  n.free <- ncol(Y) - length(par)                             # Degrees of freedom
  lags <- ncol( l.var$A ) / nrow( l.var$A )                   # Number of lags
  I <- var_lhood_grad_vcv(sim, par, lags )                    # The information matrix
  E.H <- var_lhood_hessian_numeric(sim, par, lags )                   # The expected Hessian
  # return( E.H %*% solve(I) %*% E.H )
      # If I = E.H then this should be what we return.  But it is not.  Whai!?
  return( solve(E.H) / n.free )
}

##### NOT USED #####
# fisher.info.var <- function( a, A, Sigma ){
# # Computes the fisher information matrix for the VAR
#   n.par <- length( par.to( a, A, Sigma ) )                  # The number of parameters
#   mu.d <- mu.calc.d( a, A, Sigma )                          # The derivative of the mean
#   I.A <- mu.d %*% solve(Sigma) %*% t(mu.d)                  # The part due to the dynamic parameters
#   Sigma.i <- solve(Sigma)                                   # The inverse derivative
#   I.Sig <- matrix( sapply( Sigma.i, function(s) .5 * sum(diag( Sigma.i * s ) ) ),
#                    nrow(Sigma), ncol(Sigma) )
#       # The Sigma part: .5 * trace( Sigma.i * dSigma * Sigma.i * d.Sigma )
#   info <- as.matrix( bdiag( I.A, I.Sig ) )
# }






