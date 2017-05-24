library(VARext)
context("Test the VAR likelihood function")

test_that("Single lag likelihood function matches normal error likelihood", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.4,.9), 2, 2)
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
      # The likelihood of the errors
  par <- c( a, A, Sigma[lower.tri(Sigma,TRUE)] )
      # The parameters
  l <- var_lhood( Y, par, print_level = 0 )
      # Compute the likelihood
  expect_equal( l, l.R )
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
      # The likelihood of the errors
  par <- c( a, A, Sigma[lower.tri(Sigma,TRUE)] )
      # The parameters
  l <- var_lhood( Y, par, 2, print_level = 0 )
      # Compute the likelihood
  expect_equal( l, l.R )
})

### To add: Test of numeric gradient
