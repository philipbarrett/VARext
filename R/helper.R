####################################################################################
# helper.R
#
# Various helper functions
# 20may2017
# Philip Barrett, Washington DC
#
####################################################################################

par.to.l <- function( l.var ){
# Converts to par form from list
  return(par.to(l.var$a,l.var$A,l.var$Sigma))
}

par.to <- function( a, A, Sigma ){
# Converts to par form
  par <- c( a, A, Sigma[ lower.tri(Sigma,TRUE)])
  return(par)
}

par.from <- function(par, n.var, lags ){
# Returns separate a, A, Sigma from par
  a <- par[1:n.var]
  A <- matrix( par[n.var + 1:(n.var^2*lags)], n.var, n.var * lags )
  Sigma <- matrix(NA,n.var,n.var)
  Sigma[lower.tri(Sigma,TRUE)] <- tail(par,n.var*(n.var+1)/2)
  Sigma <- pmax(Sigma, t(Sigma), na.rm=TRUE)
  return(list(a=a, A=A, Sigma=Sigma))
}

mu.calc <- function( a, A ){
# Computes the mean of a VAR with parameters a, A
  n.var <- nrow(A)
  lags <- ncol(A) / n.var
      # Problem dimensions
  A.sum <- Reduce( '+', lapply( 1:lags, function(i) A[, n.var*(i-1) + 1:n.var ] ) )
  mu <- solve( diag(n.var) - A.sum, a )
      # Create the output
  return(mu)
}

mu.calc.d <- function( a, A, Sigma ){
# Computes the numerical derivative of the mean using the ordering of the long
# form of the parameters.  Omits derivatives w.r.t Sigma terms.
  n.var <- nrow(A)                                        # Number of variables
  lags <- ncol(A) / n.var                                 # Number of lags
  inc <- 1e-06                                            # Increment
  par <- par.to( a, A, Sigma )                            # Long form of parameters
  mu <- mu.calc(a,A)                                      # The base mean
  out <- matrix( 0, (1+lags*n.var)*n.var, n.var )         # Initialize output

  for( i in 1:((1+lags*n.var)*n.var) ){
    par.inc <- par
    par.inc[i] <- par.inc[i] + inc                        # The incremented parameters
    l.par.inc <- par.from( par.inc, n.var, lags )         # As a list
    mu.inc <- mu.calc( l.par.inc$a, l.par.inc$A )         # Recompute mu
    out[i,] <- ( mu.inc - mu ) / inc                      # The derivative
  }
  return(out)
}

fitted.var <- function( Y, a, A ){
# Computes the fitted values
  n.var <- nrow(A)                                      # Number of variables
  lags <- ncol(A) / n.var                               # Number of lags
  n.pds <- ncol(Y)                                      # Number of simulation periods
  out <- matrix(0,nrow=n.var, ncol=n.pds-lags)          # Initialize output
  for( i in 1:(n.pds-lags) ){
    out[,i] <- var_sim_e( a, A, matrix(Y[,i-1+1:lags],n.var,lags),matrix(0,n.var,1))[,lags+1]
  }
  return(out)
}


