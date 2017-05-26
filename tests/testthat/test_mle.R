library(VARext)
context("Test the MLE estimation")

test_that("Single lag OLS estimation", {
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
  v <- VAR(t(sim))
      # The VAR esimate
  l.var <- var.ols( sim )
      # The list of the outputs
  a.est <- sapply( v$varresult, function(x) x$coefficients['const'] )
  A.est <- t(sapply( v$varresult, function(x) x$coefficients[1:2] ))
  Sigma.est <-var( sapply( v$varresult, function(x) x$residuals ) )
  mu.est <- solve( diag(2) - A.est, a.est )
      # Replicate the parameters
  expect_equal( l.var$a, a.est )
  expect_equal( l.var$A, A.est )
  expect_equal( l.var$Sigma, Sigma.est )
  expect_equal( l.var$mu, mu.est )
})



test_that("Single lag MLE estimation is accurate", {
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
  sim <- var_sim( a, A, y0, Sigma, 500 )
      # Create the simulation
  rownames(sim) <- c('y1','y2')
  l.var <- var.ols( sim )
      # The list of the outputs
  l.mle <- var.mle.est(sim)
      # The MLE estimate
  I <- var_lhood_grad_vcv( sim, par.to.l(l.mle) )
      # The variance of the derivative
  E.H <- var_lhood_hessian_numeric( sim, par.to.l(l.mle) )
      # The expected 2nd derivative
  #### THESE SHOULD BE EQUAL BUT ARE NOT. WHY? ####
  #### If we up the number of simulations to 1e7 then they are (more-or-less) ###
  expect_true( mean(abs( l.var$a - l.mle$a )) < 2e-09 )
  expect_true( mean(abs( l.var$A - l.mle$A )) < 1e-09 )
  expect_true( mean(abs( l.var$Sigma * 199 / 200 - l.mle$Sigma )) < 1e-07 )
      # NB: The OLS vs MLE Sigmas are different due to a degrees of freedom
      # adjustment.
})

test_that("Two lag MLE estimation is accurate", {
  a <- c( .2, .1 )
  A <- matrix( c( .9, -.2, -.1, .7, -.1, .2, .05, -.1 ), 2, 4 )
      # The VAR coefficients
  mu <- solve( diag(2) - A[,1:2] - A[,3:4], a )
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
  l.var <- var.ols( sim, 2 )
      # The list of the outputs
  l.mle <- var.mle.est(sim, 2)
      # The MLE estimate
  expect_true( mean(abs( l.var$a - l.mle$a )) < 2e-04 )
  expect_true( mean(abs( l.var$A - l.mle$A )) < 1e-04 )
  expect_true( mean(abs( l.var$Sigma * 199 / 200 - l.mle$Sigma )) < 3e-03 )
      # NB: The OLS vs MLE Sigmas are different due to a degrees of freedom
      # adjustment.
})
