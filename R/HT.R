####################################################################################
# HT.R
#
# Hypothesis test for the VAR
# 26may2017
# Philip Barrett, Washington DC
#
####################################################################################

g.deriv <- function( g, l.var, ... ){
# Computes the numerical deriviative of a restriction
  inc <- 1e-06                                    # The increment for the numerical derivative
  par <- par.to.l( l.var )                        # The long form of the parameters
  n.par <- length(par)                            # The number of parameters
  g.0 <- g(par, ...)                              # The "base" level of g
  n.r <- length(g.0)                              # The number of restirctions
  out <- matrix( 0, n.par, n.r )                  # Initialize the output

  for( i in 1:n.par ){
    par.1 <- par                                  # Initialize the incremented parameter
    par.1[i] <- par[i] + inc                      # Increment it
    g.1 <- g(par.1, ...)                          # The new constraint
    out[i,] <- ( g.1 - g.0 ) / inc                # The derivativeB
  }

  return(out)
}

g.deriv.par <- function( g, par, ... ){
# Computes the numerical deriviative of a restriction
  inc <- 1e-06                                    # The increment for the numerical derivative
  n.par <- length(par)                            # The number of parameters
  g.0 <- g(par, ...)                              # The "base" level of g
  n.r <- length(g.0)                              # The number of restirctions
  out <- matrix( 0, n.par, n.r )                  # Initialize the output

  for( i in 1:n.par ){
    par.1 <- par                                  # Initialize the incremented parameter
    par.1[i] <- par[i] + inc                      # Increment it
    g.1 <- g(par.1, ...)                          # The new constraint
    out[i,] <- ( g.1 - g.0 ) / inc                # The derivativeB
  }

  return(out)
}

wald.stat <- function( l.var, vcv, g, ... ){
# Computes the Wald statistic and critical value for the restriction
# g(a,A,Sigma) = 0 when the par-form var-covar matrix is vcv
  g.hat <- g( par.to.l( l.var ), ... )                         # The value of the restriction
  g.hat.d <- g.deriv( g, l.var, ... )            # The derivatives of the restriction
  W <- t(g.hat) %*% solve( t(g.hat.d) %*% vcv %*% g.hat.d ) %*% g.hat
      # The Wald statistic
  n.rest <- length(g.hat)                         # The number of restrictions
  p.vals <- c(.9, .95, .975, .99, .995 )
  names(p.vals) <- c(.9, .95, .975, .99, .995 )
  crit.vals <- qchisq( p.vals, n.rest )
  return( list( W=W, crit.vals=crit.vals ) )
}
