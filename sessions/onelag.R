#### Session code for a one-lag VAR ####

rm(list=ls())

# Define
a <- c(1, 2)
A <- matrix( c( .6, .1, -.2, .7), 2, 2 )
m.sd <- matrix( c( 1, .2, -.3, .5), 2, 2 )
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
irf.plot( l.var.est.mle, n.pds = 20 )

# Estimate: restricted MLE
l.var.rest.mle <- var.mle.rest( sim, mu.diff.pos, theta=.25-diff(l.var.est$mu) )
irf.plot( l.var.rest.mle, n.pds = 20 )

# Infer: OLS + MLE
v.ols <- var.ols.se( sim )
v.ols.ext <- lapply( v.ols, function(x) {
  out <- matrix( 0, 9, 9 )
  out[1:6,1:6] <- x
  return(out)
} )
v.mle <- var.mle.se( sim, l.var.est.mle )
v.mle.rest <- var.mle.se( sim, l.var.rest.mle )
print( cbind( sapply( v.ols, function(x) c( sqrt(diag(x)), rep(NA,3) ) ),
              v.mle=sqrt(diag(v.mle)), v.mle.rest=sqrt(diag(v.mle.rest)) ) )

# Plot vs. simulations
sim.fcast.plot(sim, l.var = l.var.est)
sim.fcast.plot(sim, l.var = l.var.est.mle )
sim.fcast.plot(sim, l.var = l.var.rest.mle )

# Create latex output
l.l.var <- list( l.var.est, l.var.est.mle, l.var.rest.mle )
l.m.var <- list( v.ols$ols.nw, v.mle, v.mle.rest )
v.lhood <- sapply( l.l.var, function(x) var_lhood_N( sim, par.to(x$a, x$A, x$Sigma) ) )
var.table( l.l.var, l.m.var, file='./tests/test_onelag.tex',
           specnames = c('OLS Newey-West', 'MLE', 'Restricted MLE' ),
           caption = 'ONe-lag test VAR. In restricted MLE, mean difference is one
           larger than the OLS estimate', label='tab:onelag', footer=TRUE, v.lhood=v.lhood )

# Model tests
v.theta <- seq( .01, 1, length.out=20) ^ 2
lr.wald.plot( l.var.est.mle, v.mle, mu.diff.pos, v.theta, offset=diff(l.var.est$mu) )
