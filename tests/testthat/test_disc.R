library(VARext)
context("Test the discretization objective function")

test_that("Single lag initial guess", {

  a <- c(1,1)
  A <- matrix( c(.5,.2,.1,.7), 2, 2)
      # VAR coefficients
  lags <- 1
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  l.var <- list( a=a, A=A, Sigma=Sigma, mu=c(var_lr_mu( a, A, 1 )) )
      # The VAR

  X <- var.disc.pts(l.var, lags, 2, 3, lb = NULL)
      # The grid of points
  M <- t(sapply( 1:nrow(X)-1, function(i) M_i_ig( X, i, a, A, Sigma ) ))
      # Create the initial transition probabilities
  f.0 <- disc_obj( M[1,], X, 0, a, A, Sigma, rep(1,7) )
      # objective function

  expect_equal( apply(M,1,sum), rep(1,nrow(X)) )
  expect_equal( M %*% X, var_cm(X, a, A ) )
  for( i in 1:nrow(X) ) expect_equal( disc_cv(X,M)[,i], Sigma[lower.tri(Sigma,TRUE)] )
      # Check initial guesses are ok
  expect_equal( f.0, 0 )
      # Objective function should be near zero

  M <- matrix( runif(13^2,0,1), 13,13)
  fn <- function(m) disc_obj( m, X, 0, a, A, Sigma, rep(1,7) )
  grad <- function(m) disc_grad( m, X, 0, a, A, Sigma, rep(1,7) )
      # Create & check derivatives
  cc <- check.derivatives( c(M[1,]), fn, grad, check_derivatives_print = 'none' )
      # Create numerical derivatives
  expect_true( max(abs(cc$relative_error)) < 5e-06 )
} )
