library(VARext)
context("Test the hypothesis testing")

test_that("Check the derivative function", {

  a <- c( .2, .1 )
  A <- matrix( c( .9, -.2, -.1, .7, -.1, .2, .05, -.1 ), 2, 4 )
      # The VAR coefficients
  sd <- matrix( c( 1, .2, .3, 1.5 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  l.var <- list( a=a, A=A, Sigma=Sigma )
      # Create list
  g <- function( par, theta ){
    l.var <- par.from( par, 2, 2 )
    theta - diff( mu.calc( l.var$a, l.var$A ) )
  }
      # Differene in long-run means
  expect_equal( - g(par.to.l(l.var),0), diff(mu.calc(a,A)) )
  expect_equal( - c( g.deriv( g, l.var, theta = 1 ) )[1:10],
                 apply( mu.calc.d( a, A, Sigma ), 1, diff ) )
} )
