library(VARext)
context("Test the standard error calculations")

test_that("Single lag errors", {
  a <- c( .2, .1 )
  A <- matrix( c( .9, -.2, -.1, .7 ), 2, 2 )
      # The VAR coefficients
  mu <- solve( diag(2) - A, a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .3, 1.5 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  y0 <- matrix( mu, nrow=2, ncol=1 )
      # Initial value of y
  set.seed(42)
  sim <- var_sim( a, A, y0, Sigma, 200 )
      # Create the simulation
  rownames(sim) <- c('y1','y2')
      # Prevents warnings
  l.var <- var.ols( sim, 1 )
  v.est <- VAR( t(sim), 1 )
      # The VAR esimate
  l <- summary( v.est )
      # Summary.  Used for standard errors.
  ols.var.est <- var.ols.se(sim)
      # The variance of the estimates
  R.se.ols <- c(sapply(1:2,function(i) sqrt(diag(l$varresult[[i]]$cov.unscaled *
                                                   l.var$Sigma[i,i]))))[c(3,6,1,4,2,5)]
      # The standard errors from the R estimates
  expect_equal( sqrt(diag(ols.var.est$ols)), R.se.ols )
})

test_that("Information matrix equality holds", {
  a <- c( .2, .1 )
  A <- matrix( c( .9, -.2, -.1, .7 ), 2, 2 )
      # The VAR coefficients
  mu <- solve( diag(2) - A, a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .3, 1.5 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  y0 <- matrix( mu, nrow=2, ncol=1 )
      # Initial value of y
  set.seed(42)
  e <- var_e( Sigma, 1e6 )
  sim <- var_sim_e( a, A, y0, e )
      # Create the simulation
  rownames(sim) <- c('y1','y2')
      # Prevents warnings
  l.var <- var.mle.est( sim, 1 )
  l.var.uc <- var.mle.est( sim, 1 )
      # The ml estimates for cond and uncond cases
  I <- var_lhood_grad_vcv( sim, par.to.l(l.var), 1 )
  E.H <- var_lhood_hessian_numeric(sim, par.to.l(l.var), 1 )
      # The information matrix and Hessian
  I.uc <- var_lhood_grad_vcv( sim, par.to.l(l.var.uc), 1 )
  E.H.uc <- var_lhood_hessian_numeric(sim, par.to.l(l.var.uc), 1 )
      # The information matrix and Hessian
  expect_true( mean(pmin( abs(I - E.H), abs(I/E.H-1) )) < 2e-03 )
  expect_true( mean(pmin( abs(I.uc - E.H.uc), abs(I.uc/E.H.uc-1) )) < 2e-03 )
      # Check that one of the absolute and % error are small
})
