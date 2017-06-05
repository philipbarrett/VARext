library(VARext)
context("Test the discretization")

test_that("Single lag discretization", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.1,.7), 2, 2)
      # VAR coefficients
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  lags <- 1
  l.var <- list( a=a, A=A, Sigma=Sigma, mu=c(var_lr_mu( a, A, 1 )) )
      # The VAR
  X <- var.disc.pts(l.var, lags, 2, 4, lb = NULL)
  disc <- var.disc.1( l.var, X )
      # Discretize the VAR
  set.seed(42)
  m.idx <- markov_sim( 1e6, disc$M, 0, disc$n.X )[-(1:1e5)]
  m.sim <- disc$X[m.idx+1, ]
      # Create the simulation
  l.var.est <- var.ols( t(m.sim), 1 )

  # Impose lower bound

  expect_equal( disc$X[1,], l.var$mu )
  expect_true( max( abs( apply(m.sim,2,mean) / l.var$mu - 1 ) ) < .005 )

} )

test_that("Two lag discretization", {
  a <- c(1,1)
  A <- matrix( c(.5,.2,.1,.7, -.1, .2, .2, -.2), 2, 4)
      # VAR coefficients
  sd <- matrix( c( 1, .2, .4, 2 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  l.var <- list( a=a, A=A, Sigma=Sigma, mu=c(var_lr_mu( a, A, 2 )) )
      # The VAR form
  mu <- var_lr_mu(a,A,2)
      # The mean
  set.seed(42)
  disc <- var.disc( l.var, 1, 2 )
      # Discretize the VAR
  m.idx <- markov_sim( 1e6, disc$trans, 0, disc$n.X )[-(1:1e5)]
  m.sim <- disc$X[m.idx+1, ]

  l.var.est <- var.ols( t(m.sim[,1:2]), 2 )
  # Impose lower bound

  expect_equal( disc$X[1,], l.var$mu )
} )
