library(VARext)
context("VAR simulation")

test_that("Error simulation has correct covariance", {
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  set.seed(42)
  e <- var_e( Sigma, 1e6 )
  Sigma.num <- var(t(e))
      # Create numerical equivalent
  expect_true(max(abs(Sigma - Sigma.num)) < 1e-2)
})

test_that("1-lag VAR simulation has correct coefficients", {
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
  sim <- var_sim( a, A, y0, Sigma, 1e6 )
      # Create the simulation
  rownames(sim) <- c('y1','y2')
      # PRevents warnings
  v <- VAR(t(sim))
      # Do the VAR regression
  a.num <- sapply( v$varresult, function(x) x$coefficients['const'] )
  A.num <- t(sapply( v$varresult, function(x) x$coefficients[1:2] ))
  mu.num <- apply( sim, 1, mean )
      # Compute the numerical coefficients
  expect_true( mean( abs( a - a.num ) ) < 2e-03 )
  expect_true( mean( abs( A - A.num ) ) < 1e-03 )
  expect_true( mean( abs( mu - mu.num ) ) < 5e-03 )
})

test_that("2-lag VAR simulation has correct coefficients", {
  a <- c( 1, 1 )
  A <- matrix( c( .9, -.2, -.1, .7, -.1, .2, .4, -.2 ), 2, 4 )
      # The VAR coefficients
  mu <- solve( diag(2) - A[,1:2] - A[,3:4], a )
      # The unconditional average
      # A %*% rep( mu, 2 ) + a = mu
  sd <- matrix( c( 1, .2, .3, 1.5 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  y0 <- matrix( mu, nrow=2, ncol=2 )
      # Initial value of y
  set.seed(42)
  sim <- var_sim( a, A, y0, Sigma, 1e6 )
      # Create the simulation
  rownames(sim) <- c('y1','y2')
      # Prevents warnings
  v <- VAR(t(sim), 2)
      # Do the VAR regression
  a.num <- sapply( v$varresult, function(x) x$coefficients['const'] )
  A.num <- t(sapply( v$varresult, function(x) x$coefficients[1:4] ))
  mu.num <- apply( sim, 1, mean )
      # Compute the numerical coefficients
  expect_true( mean( abs( a - a.num ) ) < 1e-02 )
  expect_true( mean( abs( A - A.num ) ) < 2e-03 )
  expect_true( mean( abs( mu - mu.num ) ) < 1e-02 )
})

