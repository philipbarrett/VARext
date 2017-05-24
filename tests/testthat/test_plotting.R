library(VARext)
context("Test the plotting functions")

test_that("One-lag conditional variance is correct", {
  a <- c( .2, .1 )
  A <- matrix( c( .9, -.2, -.1, .7 ), 2, 2 )
      # The VAR coefficients
  sd <- matrix( c( 1, .2, .3, 1.5 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  cv <- cond.var.calc( A, Sigma, 5 )
      # Conditional variance
  expect_equal( cv[[1]], Sigma )
  expect_equal( cv[[2]], Sigma + A %*% Sigma %*% t(A) )
  expect_equal( cv[[3]], Sigma + A %*% Sigma %*% t(A) + (A %^% 2) %*% Sigma %*% (t(A)%^% 2) )
  expect_equal( cv[[4]],
                Sigma + A %*% Sigma %*% t(A) + (A %^% 2) %*% Sigma %*% (t(A)%^% 2) +
                (A %^% 3) %*% Sigma %*% (t(A)%^% 3) )
})

test_that("Two-lag conditional variance is correct", {
  a <- c( .2, .1 )
  A <- matrix( c( .9, -.2, -.1, .7, -.1, .2, .05, -.1 ), 2, 4 )
      # The VAR coefficients
  sd <- matrix( c( 1, .2, .3, 1.5 ), 2, 2 )
  Sigma <- sd %*% t(sd)
      # Simulation var-covar matrix
  cv <- cond.var.calc( A, Sigma, 5 )
      # Conditional variance
  A.1 <- A[,1:2]
  A.2 <- A[,3:4]
      # Split the A matrix
  expect_equal( cv[[1]], Sigma )
  expect_equal( cv[[2]], Sigma + A.1 %*% Sigma %*% t(A.1) )
  expect_equal( cv[[3]], Sigma + A.1 %*% Sigma %*% t(A.1) +
                  (A.1 %^% 2) %*% Sigma %*% (t(A.1)%^% 2) +
                  (A.2 %^% 1) %*% Sigma %*% (t(A.2)%^% 1) )
  expect_equal( cv[[4]], Sigma + A.1 %*% Sigma %*% t(A.1) +
                  (A.1 %^% 2) %*% Sigma %*% (t(A.1)%^% 2) +
                  (A.1 %^% 3) %*% Sigma %*% (t(A.1)%^% 3) +
                  A.2 %*% Sigma %*% (t(A.2)) + (A.2 %^% 2) %*% Sigma %*% (t(A.2)%^% 2) )
})
