####################################################################################
# helper.R
#
# Various helper functions
# 20may2017
# Philip Barrett, Washington DC
#
####################################################################################

par.to <- function( a, A, Sigma ){
# Converts to par form
  par <- c( a, A, Sigma[ lower.tri(Sigma,TRUE)])
  return(par)
}

par.from <- function(par, n.var, lags ){
# Returns separate a, A, Sigma from par
  a <- par[1:n.var]
  A <- matrix( par[n.var + 1:(n.var^2*lags)], n.var, n.var * lags )
  Sigma <- matrix(0,n.var,n.var)
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

# se.calc <-

