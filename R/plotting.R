####################################################################################
# plotting.R
#
# Provides various plotting functions for a VAR
# 21may2017
# Philip Barrett, Washington DC
#
####################################################################################

cond.var.calc <- function( A, Sigma, n.pds ){
# Calculates a list of n-period-ahead conditional variances for a VAR
  n.var <- nrow(A)                                    # Number of variables
  lags <- ncol(A) / n.var                             # Number of lags
  l.A <- lapply( 1:lags, function(i) A[,(i-1)*n.var+1:n.var])
      # List of lag matrices
  l.inc <- lapply( 1:lags,
              function(i) lapply( 1:n.pds,
                function(j) (l.A[[i]] %^% j) %*% Sigma %*% (t(l.A[[i]]) %^% j )) )
      # The list of lists of increments
  out <- list()
  out[[1]] <- Sigma
  for( i in 2:n.pds ){
    out[[i]] <- out[[i-1]]
        # Cumulate
    for( j in 1:lags ){
      if( j < i ){
        out[[i]] <- out[[i]] + l.inc[[j]][[i-j]]
      }
          # Increment
    }
  }
  return(out)
}

irf.calc <- function( l.var, i, y0=NULL, n.pds=40, shk.sd=NULL, relative=FALSE ){
# Calculates the impulse responses from a VAR in dimension i
  n.var <- nrow(l.var$A)                                    # Number of variables
  lags <- ncol(l.var$A) / n.var                             # Number of lags
  if( is.null(y0) )  y0 <- matrix( mu.calc( l.var$a, l.var$A ), n.var, lags )
      # Set initial condition to the mean if required
  shk <- if(is.null(shk.sd)) .01 else shk.sd * sqrt(l.var$Sigma[i,i])
      # The initial shock, scaled by shk,sd
  e <- matrix( 0, n.var, n.pds )
  e[i,1] <- shk                                             # The shock
  fcast.mean <- var_sim_e( l.var$a, l.var$A, y0, e )        # The mean forecast
  if(relative) fcast.mean <- fcast.mean - matrix( y0, n.var, n.pds + lags )
      # Impulse response relative to y0
  fcast.sd <- cbind( NA * y0, rep(0,n.var),
        sapply( cond.var.calc( l.var$A, l.var$Sigma, n.pds - 1 ), function(x) sqrt(diag(x)) ) )
      # One standard deviation
  return( list( mean=fcast.mean, sd=fcast.sd ) )
}

irf.plot <- function( l.var, y0=NULL, n.pds=40, shk.sd=1, relative=TRUE, varnames=NULL ){
# Plots the IRF
  n.var <- nrow(l.var$A)                                            # Number of variables
  lags <- ncol(l.var$A) / n.var                                     # Number of lags
  if( is.null(varnames) ) varnames <- paste0( 'y', 1:n.var )        # Variable names
  mu <- mu.calc( l.var$a, l.var$A )                                 # The mean

  par(mfcol=c(n.var,n.var))
  for( i in 1:n.var ){
    irf <- irf.calc( l.var, i, y0, n.pds, shk.sd, relative )
        # Compute the IRFs
    for( j in 1:n.var ){
      ylim <- range( c( irf$mean[j,] + irf$sd[j,], irf$mean[j,] - irf$sd[j,] ), na.rm = TRUE )
      ylab <- if(i==1) varnames[j] else ''
      xlab <- if(j==n.var) 'Periods' else ''
      x.vals <- c( (1-lags):0, 1:n.pds )
      plot( x.vals, irf$mean[j,], lwd=2, col='blue', ylim=ylim,
            type='l', ylab=ylab, xlab=xlab )
      lines( x.vals, irf$mean[j,] + irf$sd[j,], lty=2, col='red' )
      lines( x.vals, irf$mean[j,] - irf$sd[j,], lty=2, col='red' )
      if(relative) lines( x.vals, 0 * x.vals, lty=2 ) else lines( x.vals, rep(mu[j],length(x.vals)), lty=2 )
    }
  }
  par(mfcol=c(1,1))
}

sim.fcast.plot <- function( sim, lags=1, l.var=NULL, x.vals=NULL,  n.pds=16, n.split=20,
                            start=1, varnames=NULL, plot.fit=FALSE ){
# Plots a simulation with forecasts from the VAR in l.var
  if(is.null(l.var)) l.var <- var.ols(sim, lags)        # Create VAR if none supplies
  n.var <- nrow(l.var$A)                                # Number of variables
  lags <- ncol(l.var$A) / n.var                         # Number of lags
  if( is.null(varnames) )
    varnames <- paste0( 'y', 1:n.var )                  # Variable names
  if(is.null(x.vals)) x.vals <- 1:ncol(sim)             # Create the x values if none supplied.
  mu <- mu.calc(l.var$a, l.var$A)                       # Create the mean
  fit <- fitted.var( sim, l.var$a, l.var$A )            # The fitted values

  par(mfrow=c(n.var,1))
  for( i in 1:n.var){
    plot( x.vals, sim[i,], lwd=2, type='l', col='black', xlab='', ylab=varnames[i] )
    if(plot.fit) lines( x.vals[-(1:lags)], fit[i,], col='blue', lwd=1 )
    abline( h=mu[i], lty=2, lwd=.5 )
    j <- start - lags
    while( j < ncol(sim) - n.pds - lags ){
      y0 <- matrix( sim[,j+1:lags], n.var, lags )
      irf <- irf.calc( l.var, i, y0, n.pds, 0, FALSE )
          # Create impulse response with zero shock
      lines( x.vals[j + lags - 1 + 1:n.pds], irf$mean[i,lags-1+1:n.pds], col='blue', lwd=2 )
      lines( x.vals[j + lags - 1 + 1:n.pds],
             irf$mean[i,lags-1+1:n.pds] + irf$sd[i,lags+1:n.pds], col='red', lty=2 )
      lines( x.vals[j + lags - 1 + 1:n.pds],
             irf$mean[i,lags-1+1:n.pds] - irf$sd[i,lags+1:n.pds], col='red', lty=2 )
      j <- j + n.split
    }
  }
  par(mfrow=c(1,1))
}

err.plot <- function( sim,  lags=1, l.var=NULL, x.vals=NULL, smooth=TRUE, n.pds=16,
                      varnames=NULL ){
# Plots the VAR errors with a backward-looking smoother
  if(is.null(l.var)) l.var <- var.ols(sim, lags)        # Create VAR if none supplies
  n.var <- nrow(l.var$A)                                # Number of variables
  lags <- ncol(l.var$A) / n.var                         # Number of lags
  if( is.null(varnames) )
    varnames <- paste0( 'Residuals for y', 1:n.var )                  # Variable names
  if(is.null(x.vals)) x.vals <- 1:ncol(sim)             # Create the x values if none supplied.
  fit <- fitted.var( sim, l.var$a, l.var$A )            # The fitted values
  err <- sim[,-(1:lags)] - fit                          # The errors

  par(mfrow=c(n.var,1))
  for( i in 1:n.var){
    plot( x.vals[-(1:lags)], err[i,], lwd=2, type='l', col='black', xlab='', ylab=varnames[i] )
    abline(h=0,lwd=.5)
    if(smooth){
      sm <- filter( err[i,], rep(1/n.pds,n.pds), sides=1 )
      lines( x.vals[-(1:lags)], sm, col='blue', lwd=2 )
      if(i==1)
        legend( 'topright', paste0(n.pds, ' period backward-looking average' ),
              col='blue', lwd=2, bty='n' )
    }
  }
  par(mfrow=c(1,1))
}

lr.wald.plot <- function( sim, lags, g, v.theta, offset=0, xlab=NULL,
                          p.vals=c(.9, .95, .975, .99, .995 ), drop=NULL, ... ){
# Plots the Wald and LR test as a function of a constraint defined by theta
  n.var <- nrow(sim)
      # Problem dimensions
  l.mle <- var.mle.est( sim, lags )
  l.mle.uc <- var.mle.est( sim, lags, cond=FALSE )
      # Unconstrained estimates, conditional and unconditional
  v.mle <- var.mle.se( sim, l.mle )
  v.mle.uc <- var.mle.se( sim, l.mle.uc, cond=FALSE )
      # Local variances: SHOULD INCLUDE COND SOMEHOW
  lhood.u <- var_lhood_N( sim, par.to.l(l.var.est.mle), lags )
  lhood.u.uc <- var_lhood_N( sim, par.to.l(l.mle.uc), lags, cond=FALSE )
      # The unrestricted likelihoods (conditioanl and unconditional)
  l.mle.uc.c <- l.mle.c <- l.mle
      # Initialize constrained MLE at unconstrained
  crit.vals <- qchisq( p.vals, 1 )
      # The critical values: Only allow restrictions to vary in one dimension
  v.lhood.uc <- v.lhood.t <- v.wald.uc <- v.wald.t <- NULL

  for( thet in v.theta ){
    v.wald.t <- c(v.wald.t, wald.stat( l.mle, v.mle, g, lags=lags,
                                       n.var=n.var, theta=thet-offset, ... )$W )
    v.wald.uc <- c(v.wald.uc, wald.stat( l.mle.uc, v.mle.uc, g, lags=lags,
                                       n.var=n.var, theta=thet-offset, ... )$W )
    l.mle.c <- var.mle.rest( sim, g, lags, par=par.to.l(l.mle.c), theta=thet-offset )
    lhood.c <- var_lhood_N( sim, par.to.l(l.mle.c), lags, cond=TRUE )
    v.lhood.t <- c( v.lhood.t, 2 * ( lhood.c - lhood.u ) )
        # The conditional likelihood ratio test
    l.mle.uc.c <- var.mle.rest( sim, g, lags, par=par.to.l(l.mle.uc), theta=thet-offset,
                                cond=FALSE )
    lhood.uc <- var_lhood_N( sim, par.to.l(l.mle.uc.c), lags, cond=FALSE )
    v.lhood.uc <- c( v.lhood.uc, 2 * ( lhood.uc - lhood.u.uc ) )
        # The unconditional likelihood ratio test
  }

  if(is.null(xlab)) xlab <- expression(theta)
  plot( v.theta, v.wald.t, type='l', lwd=2, col='blue', lty=2,
        xlab=xlab, ylab='Test statistic',
        ylim=c(0, max( 1.3 * max(crit.vals), min( 2 * max(crit.vals), max( v.wald.uc ) ) ) ) )
  lines( v.theta, v.wald.uc, lwd=2, col='blue' )
  v.lhood.t[v.lhood.t <= -1e-05] <- NA
  lines( v.theta, v.lhood.t, lwd=2, col='red', lty=2 )
  if( is.null(drop) ){
    lines( v.theta, v.lhood.uc, lwd=2, col='red' )
  }else{
    lines( v.theta[-drop], v.lhood.uc[-drop], lwd=2, col='red' )
  }
  abline(h=crit.vals, lty=2)
  text( min(v.theta)+.1*abs(diff(range(v.theta))), crit.vals, paste0( 'p-value = ', p.vals ),
        pos = 3, offset=0.2, adj=c(0,0)  )
  legend( 'topleft', c('Unconditional Wald', 'Conditional Wald',
                       'Unconditional likelihood ratio', 'Conditional likelihood ratio'),
          col=c('blue','blue', 'red', 'red'), lty=c(1,2), lwd=2, bty='n')

}






