####################################################################################
# mle.R
#
# Computes the unrestricted maximum likelihood of a VAR
# 20may2017
# Philip Barrett, Washington DC
#
####################################################################################

var.ols <- function( Y, lags=1 ){
# Computes the VAR estimation via OLS
  n.var <- nrow(Y)
  if(is.null(rownames(Y))) rownames(Y) <- paste0('y', 1:n.var)
      # Add some row names to revent errors
  v.est <- VAR(t(Y), lags )
      # The VAR estimates
  a <- sapply( v.est$varresult, function(x) x$coefficients['const'] )
  A <- t( sapply( v.est$varresult, function(x) x$coefficients[1:(n.var*lags)] ) )
  Sigma <- var(sapply(v.est$varresult, function(x) x$residuals ) )
  mu <- mu.calc(a,A)
      # Extrat and format the output
  return(list(a=a, A=A, Sigma=Sigma, mu=mu))
}

var.mle.est <- function( Y, lags=1, par=NULL ){
# Compute the unresticted MLE for the VAR
  if(is.null(par)){
    l.var <- var.ols( Y, lags )
    par <- par.to( l.var$a, l.var$A, l.var$Sigma )
  }
      # Fill in initial pramaeters if required
  obj_fn <- function(x) var_lhood( Y, x, lags )               # The objective function
  grad_fn <- function(x) var_lhood_grad( Y, x, lags )         # The derivative function
  opts <- list(algorithm="NLOPT_LD_SLSQP", xtol_abs=1e-8,
               xtol_rel=1e-8, maxeval=50000 )                 # The options for nlopt
  opt <- nloptr( par, obj_fn, grad_fn, opts=opts )            # The numerical optimum
  grad <- var_lhood_grad( Y, opt$solution, lags )             # Evaluate the derivative
  out <- par.from( opt$solution, nrow(Y), lags )              # Convert to a, A, Sigma form
  out$lhood <- obj_fn(opt$solution)                           # Return the likelihood
  out$grad <- grad                                            # Return the gradient too
  out$mu <- mu.calc( out$a, out$A )                           # The long run mean
  out$status <- opt$status
  out$message <- opt$message
  return(out)
}

mu.diff.pos <- function(par, lags, n.var=2, theta=0){
# Imposes restriction that the mean difference in a 2-variable VAR is greater than theta.
  l.var <- par.from(par, n.var, lags)
  mu <- mu.calc(l.var$a, l.var$A)
  return( - ( mu[1] - mu[2] - theta ) )
}

var.mle.rest <- function( Y, g, lags=1, par=NULL, ... ){
# Compute the restricted VAR, maximizing the log likelihood subject to g(par,...) < 0
  if(is.null(par)){
    l.var <- var.ols( Y, lags )
    par <- par.to( l.var$a, l.var$A, l.var$Sigma )
  }
      # Fill in initial pramaeters if required
  obj_fn <- function(x) var_lhood( Y, x, lags )               # The objective function
  grad_fn <- function(x) var_lhood_grad( Y, x, lags )         # The derivative function
  constraint_fn <- function(x) g( x, lags, ... )
  constraint_grad_fn <- function(x) g.deriv.par( g, x, lags=lags, ... )
  opts <- list(algorithm="NLOPT_LD_SLSQP", xtol_abs=1e-8,
               xtol_rel=1e-8, maxeval=50000 )                  # The options for nlopt
  opt <- nloptr( par, obj_fn, grad_fn,
                 eval_g_ineq = constraint_fn, eval_jac_g_ineq = constraint_grad_fn,
                 opts=opts )
      # The numerical optimum
  grad <- var_lhood_grad_numeric( Y, opt$solution, lags )     # Evaluate the derivative
  out <- par.from( opt$solution, nrow(Y), lags )              # Convert to a, A, Sigma form
  out$lhood <- obj_fn(opt$solution)                           # Return the likelihood
  out$grad <- grad                                            # Return the gradient too
  out$constraint <- constraint_fn( opt$solution )              # And the constraint value
  out$mu <- mu.calc( out$a, out$A )                           # The long run mean
  out$status <- opt$status
  out$message <- opt$message
  return(out)
}




