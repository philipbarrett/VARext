library(VARext)
context("Test the VAR likelihood function")

test_that("Single lag likelihood function matches normal error likelihood", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.1,.7), 2, 2)
      # VAR coefficients
  mu <- solve( diag(2) - A, a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  set.seed(42)
  e <- var_e( Sigma, 400 )
      # Simulate the errors
  y0 <- matrix( mu, nrow=2, ncol=1 )
      # Initial value of y
  Y <- var_sim_e( a, A, y0, e )
      # The simulation
  l.R <- - sum( dmvnorm( t(e), rep(0,2), Sigma, log=TRUE ) ) / ncol(e)
      # The likelihood of the errors.  Need subtract #lags * #vars
  par <- c( a, A, Sigma[lower.tri(Sigma,TRUE)] )
      # The parameters
  l <- var_lhood( Y, par, print_level = 0 )
      # Compute the likelihood
  l.d <- var_lhood_grad( Y, par )
  l.d.n <- var_lhood_grad_numeric( Y, par )
      # The numeric and analytic gradients
  expect_equal( l, l.R )
  expect_true( mean( abs( l.d - l.d.n ) ) < 1e-05 )
})

test_that("Two lag likelihood function matches normal error likelihood", {
  a <- c(1,1)
  A <- matrix( c( .6, -.2, -.1, .5, .1, .05, -.01, .2 ), 2, 4 )
      # The VAR coefficients
  mu <- solve( diag(2) - A[,1:2] - A[,3:4], a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  set.seed(42)
  e <- var_e( Sigma, 400 )
      # Simulate the errors
  y0 <- matrix( mu, nrow=2, ncol=2 )
      # Initial value of y
  Y <- var_sim_e( a, A, y0, e )
      # The simulation
  l.R <- - sum( dmvnorm( t(e), rep(0,2), Sigma, log=TRUE ) ) / ncol(e)
      # The likelihood of the errors.  Need subtract #lags * #vars
  par <- c( a, A, Sigma[lower.tri(Sigma,TRUE)] )
      # The parameters
  l <- var_lhood( Y, par, 2, print_level = 0 )
      # Compute the likelihood
  l.d <- var_lhood_grad( Y, par, 2 )
  l.d.n <- var_lhood_grad_numeric( Y, par, 2 )
      # The numerc and analytic gradients
  expect_equal( l, l.R )
  expect_true( mean( abs( l.d - l.d.n ) ) < 1e-05 )
})

test_that("The long-run variance calculation is accurate for one lag", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.1,.7), 2, 2)
      # VAR coefficients
  mu <- solve( diag(2) - A, a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  y0 <- matrix( mu, nrow=2, ncol=1 )
      # Initial value of y
  set.seed(42)
  sim <- var_sim( a, A, y0, Sigma, 1e6 )
      # The simulation
  lr.var <- var_lr_variance( A, Sigma, 1 )
      # The long-run variance
  expect_true( max( abs( lr.var / var(t(sim)) - 1 ) ) < .01 )
})

test_that("The long-run variance calculation is accurate for two lags", {
  a <- c(1,1)
  A <- matrix( c( .6, -.2, -.1, .5, .1, .05, -.01, .2 ), 2, 4 )
      # The VAR coefficients
  mu <- solve( diag(2) - A[,1:2] - A[,3:4], a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
        # Simulation var-covar matrix
  y0 <- matrix( mu, nrow=2, ncol=2 )
      # Initial value of y
  set.seed(42)
  sim <- var_sim( a, A, y0, Sigma, 1e7 )
      # The simulation
  lr.var <- var_lr_variance( A, Sigma, 2 )
      # The long-run variance
  expect_true( max( abs( lr.var / var(t(sim)) - 1 ) ) < .02 )
})

test_that("The unconditional likelihood for one lag is accurate",{
  a <- c(1,1)
  A <- matrix( c(.5,.2,.1,.7), 2, 2)
      # VAR coefficients
  mu <- solve( diag(2) - A, a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  lr.var <- var_lr_variance( A, Sigma, 1 )
      # The lnong-run variance
  set.seed(42)
  e <- var_e( Sigma, 400 )
      # Simulate the errors
  y0 <- matrix( rmvnorm( 1, mean=mu, sigma=lr.var), nrow=2, ncol=1 )
      # Initial value of y
  Y <- var_sim_e( a, A, y0, e )
      # The simulation
  l.R <- - sum( dmvnorm( t(e), rep(0,2), Sigma, log=TRUE ) ) / ncol(e)
  l.R <- l.R - dmvnorm( t(y0), mu, lr.var, log = TRUE ) / ncol(e)
      # The likelihood of the errors.  Need subtract #lags * #vars
  par <- par.to( a, A, Sigma )
      # The parameters
  l <- var_lhood( Y, par, print_level = 0, cond = FALSE )
      # Compute the likelihood
  l.d <- var_lhood_grad( Y, par, lags=1, cond = TRUE )
  l.d.n <- var_lhood_grad_numeric( Y, par, lags=1, cond = TRUE )
  l.d.u <- var_lhood_grad( Y, par, lags=1, cond = FALSE )
  l.d.n.u <- var_lhood_grad_numeric( Y, par, lags=1, cond = FALSE )
      # The numeric and analytic gradients
  par.2 <- par
  par.2[1] <- par[1] + 1e-06

  expect_equal( l, l.R )
  expect_true( max( abs( mu.calc.d( a, A, Sigma ) /
                           var_lr_mu_grad( a, A, 1 ) - 1 ) ) < 1e-5 )
      # The derivative of the LR mean
  expect_equal( c( l.d.u - l.d.n.u  - ( l.d - l.d.n ) ), rep( 0, 9 ) )
      # Checking that the error on the derivative is not substantially increased
      # by making it unconditional
})


test_that("The unconditional likelihood for two lags is accurate",{
  a <- c(1,1)
  A <- matrix( c( .6, -.2, -.1, .5, .1, .05, -.01, .2 ), 2, 4 )
      # The VAR coefficients
  mu <- solve( diag(2) - A[,1:2] - A[,3:4], a )
      # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  lr.var <- var_lr_variance( A, Sigma, 2 )
      # The lnong-run variance
  set.seed(42)
  e <- var_e( Sigma, 400 )
      # Simulate the errors
  y0 <- matrix( rmvnorm( 1, mean=mu, sigma=lr.var), nrow=2, ncol=2 )
      # Initial value of y
  Y <- var_sim_e( a, A, y0, e )
      # The simulation
  l.R <- - sum( dmvnorm( t(e), rep(0,2), Sigma, log=TRUE ) ) / ncol(e)
  l.R <- l.R - dmvnorm( Y[,2], mu, lr.var, log = TRUE ) / ncol(e)
      # The likelihood of the errors.  Need subtract #lags * #vars
  par <- par.to( a, A, Sigma )
      # The parameters
  l <- var_lhood( Y, par, 2, print_level = 0, cond = FALSE )
      # Compute the likelihood
  l.d <- var_lhood_grad( Y, par, lags=2, cond = TRUE )
  l.d.n <- var_lhood_grad_numeric( Y, par, lags=2, cond = TRUE )
  l.d.u <- var_lhood_grad( Y, par, lags=2, cond = FALSE )
  l.d.n.u <- var_lhood_grad_numeric( Y, par, lags=2, cond = FALSE )
      # The numeric and analytic gradients
  expect_equal( l, l.R )
  expect_true( max( abs( mu.calc.d( a, A, Sigma ) /
                           var_lr_mu_grad( a, A, 2 ) - 1 ) ) < 1e-5 )
      # The derivative of the LR mean
  expect_equal( c( l.d.u - l.d.n.u  - ( l.d - l.d.n ) ), rep( 0, 13 ) )
      # Checking that the error on the derivative is not substantially increased
      # by making it unconditional
})
