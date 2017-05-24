library(VARext)
context("Test the helper functions")

test_that("par conversion works with one lag", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.4,.9), 2, 2)
      # VAR coefficients
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  par <- par.to(a,A,Sigma)
  l.pars <- par.from(par,2,1)
      # Convert back and forth
  expect_equal(a, l.pars$a)
  expect_equal(A, l.pars$A)
  expect_equal(Sigma, l.pars$Sigma)
})

test_that("par conversion works with two lags", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.4,.7, .2, -.1, .1, -.2), 2, 4)
      # VAR coefficients
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  par <- par.to(a,A,Sigma)
  l.pars <- par.from(par,2,2)
      # Convert back and forth
  expect_equal(a, l.pars$a)
  expect_equal(A, l.pars$A)
  expect_equal(Sigma, l.pars$Sigma)
})

test_that("mu.calc is accurate", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.4,.7, .2, -.1, .1, -.2), 2, 4)
      # VAR coefficients
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  mu <- mu.calc(a,A)
  y0 <- matrix(mu, nrow = 2, ncol = 2)
  sim <- var_sim( a, A, y0, Sigma, 1e6 )
      # Create the simulation
})

test_that("Fitted values errors are the same as those fed into the model", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,-.4,.8), 2, 2)
      # VAR coefficients
  mu <- mu.calc(a,A)
     # The unconditional average
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  e <- var_e( Sigma, 400 )
      # Simulate the errors
  y0 <- matrix( mu, nrow=2, ncol=1 )
      # Initial value of y
  Y <- var_sim_e( a, A, y0, e )
      # The simulation
  fit <- fitted.var( Y, a, A )
      # Fit the VAR
  e.0 <- Y[,-1] - fit
      # The errors of the fitted values
  expect_equal( e, e.0 )
})
