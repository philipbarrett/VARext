#### Session code for a one-lag VAR ####

rm(list=ls())
library(VARext)

# Define
a <- c(1, 2)
A <- matrix( c( .6, .1, -.2, .7), 2, 2 )
m.sd <- matrix( c( 1, .2, -.3, .5), 2, 2 )
lags <- 1
Sigma <- m.sd %*% t(m.sd)
mu <- mu.calc( a, A )
l.var <- list( a=a, A=A, Sigma=Sigma, mu=mu )

# Simulate
n.sim <- 200
y0 <- matrix( mu, 2, 1 )
sim <- var_sim( a, A, y0, Sigma, n.sim )

# Estimate: OLS
l.var.est <- var.ols( sim )
irf.plot( l.var, n.pds = 20 )
irf.plot( l.var.est, n.pds = 20 )

# Estimate: MLE
l.var.est.mle <- var.mle.est( sim )
l.var.est.mle.u <- var.mle.est( sim, cond=FALSE )
irf.plot( l.var.est.mle, n.pds = 20 )
irf.plot( l.var.est.mle.u, n.pds = 20 )

# Estimate: restricted MLE
l.var.rest.mle <- var.mle.rest( sim, mu.diff.pos, theta=2-diff(l.var.est$mu) )
l.var.rest.mle.u <- var.mle.rest( sim, mu.diff.pos, cond=FALSE, theta=2-diff(l.var.est$mu) )
irf.plot( l.var.rest.mle, n.pds = 20 )
irf.plot( l.var.rest.mle.u, n.pds = 20 )

# Infer: OLS + MLE
v.ols <- var.ols.se( sim )
v.ols.ext <- lapply( v.ols, function(x) {
  out <- matrix( 0, 9, 9 )
  out[1:6,1:6] <- x
  return(out)
} )
v.mle <- var.mle.se( sim, l.var.est.mle )
v.mle.est.ineff <- var.mle.se( sim, l.var.est.mle, ineff = TRUE )
v.mle.rest <- var.mle.se( sim, l.var.rest.mle )
v.mle.rest.u <- var.mle.se( sim, l.var.rest.mle.u, cond = FALSE, ineff=TRUE )
print( cbind( sapply( v.ols, function(x) c( sqrt(diag(x)), rep(NA,3) ) ),
              v.mle=sqrt(diag(v.mle)), v.mle.est.ineff=sqrt(diag(v.mle.est.ineff)),
              v.mle.rest=sqrt(diag(v.mle.rest)), v.mle.rest.u=sqrt(diag(v.mle.rest.u)) ) )

# Plot vs. simulations
sim.fcast.plot(sim, l.var = l.var.est)
sim.fcast.plot(sim, l.var = l.var.est.mle )
sim.fcast.plot(sim, l.var = l.var.rest.mle )
sim.fcast.plot(sim, l.var = l.var.rest.mle.u )

# Create latex output
l.l.var <- list( l.var.est, l.var.est.mle, l.var.est.mle, l.var.rest.mle, l.var.rest.mle.u )
l.m.var <- list( v.ols$ols.nw, v.mle, v.mle.est.ineff, v.mle.rest, v.mle.rest.u )
v.cond <- c( T, T, T, T, F )
v.lhood <- sapply( 1:5, function(i) var_lhood_N( sim, par.to.l(l.l.var[[i]]), 1, cond=v.cond[i] ) )
var.table( l.l.var, l.m.var, file='./tests/test_onelag.tex',
           specnames = c('OLS Newey-West', 'MLE', 'MLE w/o enfficiency', 'Restricted MLE', 'Restricted UMLE' ),
           caption = 'One-lag test VAR. In restricted MLE, mean difference is one
           larger than the OLS estimate', label='tab:onelag', footer=TRUE, v.lhood=v.lhood )

# Model tests
v.theta <- seq( .01, 1, length.out=20) ^ 2
lr.wald.plot( sim, lags, mu.diff.pos, v.theta, offset=diff(l.var.est$mu) )

# Discretization
disc <- var.discretize( l.var, 1, c(-Inf,-Inf) )
    # Should move this to a testing file
