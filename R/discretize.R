#####################################################################
# discretize.R
#
# Contains the functions for discretization of the VAR and related
#
# 24apr2017
# Philip Barrett, Washington DC
#####################################################################

var.disc.pts <- function(l.var, lags, n.pts, n.dirs, lb = NULL){
  # Creates the nodes for the discretization
  n.var <- length(l.var$a)
  # Number of variables
  Sig.u <- var_lr_variance( l.var$A, l.var$Sigma, lags )
  # The unconditional variance
  # l.0 <- qnorm( seq( .5 + 4 * q / 5, 1 - q / 5, length.out=n.pts ), 0, 1 )
  n.pts <- n.pts * 2 + 1
  l.0 <- seq( - sqrt(n.pts-1), sqrt(n.pts-1), length.out = n.pts )
  l.0 <- l.0[ l.0 != 0 ]
  # Vector of distances (strip out zero)
  theta.0 <- seq( 0, pi, length.out=n.dirs+1 )[-(n.dirs+1)]
  # The vector of angles
  zz.0 <- cbind( sin(theta.0), cos(theta.0) )
  pts.0 <- rbind( c(0,0),
                  do.call( 'rbind', lapply( l.0, function(l) l * zz.0 ) ) )
  # The matrix of points
  n.X.0 <- nrow(pts.0)
  X.0 <- t( t( chol( Sig.u ) ) %*% t( pts.0 ) ) + rep(1, n.X.0 ) %*% t(l.var$mu)
  # Rescale the points to the mean and unconditional variance
  if( !is.null(lb) ) X.0 <- pmax( X.0, rep(1,n.X.0) %*% t(lb) )
  # Impose the lower bound
  X <- t(apply( expand.grid( lapply( 1:lags, function(i) 1:n.X.0 ) ),
                1, function(x) c( t(X.0[x,]) ) ))
  # The grid including lags
  X <- X[order(as.integer(factor(apply(X[,-(1:n.var)], 1, toString)))),]
  # Order the points by lagged variables
  return( X )
}

var.disc.1.i <- function( l.var, X, i, m.i=NULL, w=NULL ){
# Computes the ith row of the probability matrix minimizing the errors on the
# target moments of a 1-lag VAR
  a <- l.var$a ; A <- l.var$A ; Sigma <- l.var$Sigma
      # Extract VAR coefficients
  n.var <- length(a) ; lags <- ncol(A) / n.var          # Problem dimensions
  n.m <- nrow(X)                                        # Number of points
  # if( is.null(m.i) ) m.i <- rep( 1/ n.m, n.m )
  if(class(X)!='matrix') X <- as.matrix(X)
      # SOmetimes need as X can get converted to dataframe in higher functions
  if( is.null(m.i) ) m.i <- M_i_ig( X, i-1, a, A, Sigma )
      # Initial guess.  Convert i index to c++ style counting
  if( is.null(w) ){
    n.constr <- n.var * 2 + .5 * n.var * (n.var + 1 )
    w <- rep(1,n.constr)
    w[ n.constr - n.var + 1:n.var ] <- .1
        # Down-weight the skew
  }
      # Defautl weighting vector
  obj <- function(m) disc_obj( m, X, i-1, a, A, Sigma, w )
  obj.grad <- function(m) disc_grad( m, X, i-1, a, A, Sigma, w )
      # Objective and gradient
  obj.1 <- function(m.1) disc_obj( c( m.1, 1-sum(m.1) ), X, i-1, a, A, Sigma, w ) +
                            1e7 * ( sum(m.1) > 1 ) * ( sum(m.1) - 1 ) ^ 2
    # Version with constraint imposed internally
  consrt <- function(m) sum( m ) - 1
  consrt.grad <- function(m) rep(1,n.m)
      # Constraint and gradient
  lb <- rep( 0, n.m ) ; ub <- rep( 1, n.m ) ; m.0 <- pmax( pmin( m.i, 1 ), 0 )
      # Bounds, including on initial guess
  # m.0[n.m] <- 1 - sum( m.0[-n.m] )
  #     # Guarantee constraint holds

  opts <- list(algorithm="NLOPT_LD_SLSQP", xtol_abs=1e-12,
               xtol_rel=1e-12, maxeval=500 ) #, check_derivatives=TRUE )
      # The options for nlopt
  opt <- nloptr( m.0, obj, obj.grad, lb, ub, eval_g_eq = consrt,
                 eval_jac_g_eq = consrt.grad, opts=opts )
      # The numerical optimum

  if( opt$status < 0 ){
    opts.1 <- list(algorithm="NLOPT_LN_COBYLA", xtol_abs=1e-5, xtol_rel=1e-5, maxeval=5000 )
    opt.1 <- nloptr( m.0[-n.m], obj.1, lb=lb[-n.m], ub=ub[-n.m], opts=opts.1 )
        # Approximate solution w/o derviatives
    opt <- nloptr( c(opt.1$solution, max( 1-sum(opt.1$solution), 0 ) ),
                   obj, obj.grad, lb, ub, eval_g_eq = consrt,
                   eval_jac_g_eq = consrt.grad, opts=opts )
    # if(sum(opt$solution) > 1 ) browser()
  }   # Fail-safe if cannot find an optimum to start with

  out <- list( m=opt$solution / sum(opt$solution), obj=obj(opt$solution),
               obj.grad=obj.grad(opt$solution),
               consrt=consrt(opt$solution), consrt.grad=consrt.grad(opt$solution),
               cm=t(opt$solution) %*% X, cm.targ=t(a + A %*% X[i,]),
               cv=disc_cv(X,matrix(opt$solution, n.m, n.m, byrow = TRUE ) )[,1],
               cv.targ=Sigma[lower.tri(Sigma,TRUE)],
               skew=disc_skew(X,matrix(opt$solution, n.m, n.m, byrow = TRUE ) )[,1] )
  return(out)
}

var.disc.1 <- function( l.var, X, M=NULL, w=NULL ){
# Create the discretize one-lag VAR
  if( is.null(M)){
    dta <- lapply( 1:nrow(X), function(i) var.disc.1.i( l.var, X, i, w=w ) )
  }else{
    dta <- lapply( 1:nrow(X), function(i) var.disc.1.i( l.var, X, i, M[i,], w ) )
  }
      # Create the approximations for each row
  cm.err <- sapply( dta, function(x) x$cm - x$cm.targ )
  cv.err <- sapply( dta, function(x) x$cv - x$cv.targ )
  skew.err <- sapply( dta, function(x) x$skew - x$skew )
  prob <- sapply( dta, function(x) sum( x$m ) )
      # Create the errors
  M <- t( sapply( dta, function(x) x$m ) )
      # Transition matrix
  return( list( X=X, M=M, n.X=nrow(X), cm.err=cm.err, cv.err=cv.err,
                skew.err=skew.err, prob=prob ) )
}

var.disc.lag <- function( l.var, X, M=NULL, w=NULL ){
# Computes the transition matrices for a set of points X where all but the t-1
# lag are identical

  a <- l.var$a ; A <- l.var$A ; Sigma <- l.var$Sigma      # Extract VAR coefficients
  n.var <- length(a) ; lags <- ncol(A) / n.var          # Problem dimensions
  l.var.cpy <- l.var                                      # Make  copy of the VAR
  l.var.cpy$a <- c( a + A[,-(1:n.var)] %*% t(X[1,-(1:n.var)]) )
  l.var.cpy$A <- A[,1:n.var]
      # Adjust a & A to account for the change in the prediction due to the lags
      # before period t-1
  disc <- var.disc.1( l.var.cpy, X[,1:n.var], M, w )
      # The discretization
  return(disc)
}

var.disc.trans <- function( l.var, X, M=NULL, w = NULL ){
# Puts the discretiation of the lagged process in the grand transition matrix

  n.pts <- nrow(X)                                        # Number of points (total)
  a <- l.var$a ; A <- l.var$A ; Sigma <- l.var$Sigma      # Extract VAR coefficients
  n.var <- length(a) ; lags <- ncol(A) / n.var            # Problem dimensions
  lag.idx <- as.integer(factor(apply(X[,-(1:n.var)], 1, toString)))
      # Index of lags
  # new.idx <- as.integer(factor(apply(X[,-(n.var*(lags-1)+1:n.var)], 1, toString)))
  #     # Index of realizations
  l.disc <- by( X, lag.idx, var.disc.lag, l.var=l.var, M=M, w=w )
      # The list of the transition matrices
  out <- matrix( 0, nrow(X), nrow(X))
      # The grand transition matrix
  counter <- 1
  for( i in 1:length(l.disc) ){
    Z <- X[ lag.idx==i, ]
    M <- l.disc[[i]]$M
    for( j in 1:nrow(Z) ){
      for( k in 1:ncol(M) ){
        target <- c( Z[k,1:n.var], Z[j,-((lags-1)*n.var+1:n.var)] )
            # The target to match
        match.idx <- which( apply( X, 1, function(x) all(x==target) ) )
            # The index of the match
        out[counter,match.idx] <- M[j,k] / sum(M[j,])
            # Assign to the output matrix, normalize probabilities to one where required
      }
      counter <- counter + 1                      # Increment the counter
    }
  }
  return(out)
}

var.disc <- function( l.var, n.pts, n.dirs, n.min.pts=5, lb=NULL, p.min=1e-05 ){
# Discretizes a VAR
  a <- l.var$a ; A <- l.var$A ; Sigma <- l.var$Sigma          # Extract VAR coefficients
  n.var <- length(a) ; lags <- ncol(A) / n.var                # Problem dimensions

  #### Create nodes ####
  X <- var.disc.pts(l.var, lags, n.pts, n.dirs, lb )
      # Initial discretization

  #### Iterate over solutions to strip out unecessary points ####
  n.X <- nrow(X) ; n.X.new <- n.X + 1
      # Initialize loop variables
  while( n.X != n.X.new ){
    n.X <- nrow(X)
    M <- var.disc.trans(l.var, X )                # Update the transition probabilities
    p.lr <- (M %^% 200)[1,]                       # The long-run distribution
    X <- X[ p.lr > p.min, ]                       # Eliminate rare points
    n.X.new <- nrow(X)
  }
  trans <- M / apply( M, 1, sum )
      # Regularize the probabilties just in case

  return( list( X=X, trans=trans, n.X=n.X.new, p.lr=p.lr ) )
}





var.discretize <- function(l.var, lags, lb=c(-Inf,-Inf), n.int=1e6, n.pts = 3, n.dirs = 8 ){
# Discretizes a VAR object

  X <- var.disc.pts(l.var, lags, n.pts, n.dirs, lb)
      # The grid
  n.X <- nrow(X)

  ## STATE REDUCTION ##
  # browser()
  NN <- 1e7 ; p <- .01 / 100
      # Number of simulation periods and cut-off (default 0.01%)
  sim <- var_sim( l.var$a, l.var$A, matrix(l.var$mu,nn,lags), l.var$Sigma, NN)
  if( lags > 1 ){
    for( i in 2:lags){
      sim <- rbind( sim[,-1], sim[,-ncol(sim)])
    }
  }
      # The simulation: HOW TO FIX LAGS HERE?
  nearest <- nn2( X, t(sim), 1 )
      # The nearest neighbors
  freq <- table( nearest$nn.idx ) / NN
      # Table of fequencies of the nodes
  X <- X[ as.numeric(names(freq))[freq>=p], ]
      # Select the points that show up at least twice
  n.X <- nrow(X)
      # Update number of points

  int.eps <- rmvnorm( n.int, 0 * l.var$a, l.var$Sigma )
      # The vector of integration shocks
  trans <- matrix(0, n.X, n.X)
      # Initiate the transition probability matrix
  for( i in 1:n.X ){
    X.cont <- rep(1,n.int) %*% t( l.var$a + l.var$A %*% X[i,] ) + int.eps
        # The stochastic continuation values of X
    nearest <- nn2( X, X.cont, 1 )
        # The list of nearest neighbours
    trans[i,] <- table(factor(nearest$nn.idx,levels=1:n.X) ) / n.int
  }
  m <- (trans %^% 100)[1,]
      # The unconditional distribution


  # rmg <- matrix( X[,'rfr'], nrow(X), nrow(X) ) - matrix( X[,'gth'], nrow(X), nrow(X), byrow=TRUE )
  # m.rmg <- m * trans
  # d.rmg <- cbind( c(rmg), c(m.rmg) )
  # d.rmg <- d.rmg[ order(d.rmg[,1]), ]
  #     # The unconditional distribution for r minus g
  ### Should be a separate function ###

  return( list( X=X, trans=trans, m=m, n.X=n.X ) ) #, d.rmg=d.rmg ) )
}

# var.disc.mean.var <- function( trans, pts ){
# # Computes the conditional means, conditional variances and unconditional
# # variance of a VAR discretization
#   lr.d <- (trans %^% 100)[1,]                               # The long run distribution
#   lr.mu <- c( lr.d %*% pts )                                # The unconditional mean
#   cm <- trans %*% pts                                       # The conditional mean
#   pts.dev <- ( pts - rep(1,nrow(pts)) %*% t(lr.mu) )        # Deviations from the mean
#   l.sq <- lapply( 1:nrow(trans), function(i) pts.dev[i,] %*% t(pts.dev[i,])  )
#       # The squared realizations of the deviations
#   lr.v <- Reduce(`+`,Map(`*`, l.sq, lr.d ) )                # The long-run variance
#   cv <- sapply( 1:nrow(trans), function(i)
#     c( Reduce(`+`, Map(`*`, l.sq, trans[i,] ) )[lower.tri(l.sq[[i]], TRUE)] ) )
#       # The conditional covariances
#   return( list( cm=cm, cv=cv, lr.mu=lr.mu, lr.v=lr.v ) )
# }  ## NO GOOD.  CONDITIONAL VARIANCES ARE WRONG

var.disc.mean.var.err <- function( trans, pts, a, A, Sigma ){
# Computes the error on the discretization
  lags <- ncol(A) / nrow(A)
  cm <- t( a + A %*% t( pts ) )                             # The conditional mean
  cmv.disc <- var.disc.mean.var( trans, pts )               # The conditional mean & variance
  err <- c( cmv.disc$cm - cm,
          c(cmv.disc$cv) - c( Sigma[lower.tri(Sigma, TRUE)] ) ) # The error
  return(err)
}

var.disc.mean.var.err.grad <- function( trans, pts, a, A, Sigma ){
# Computes the gradient on the error on the discretization
  lags <- ncol(A) / nrow(A)
  lr.d <- (trans %^% 100)[1,]                               # The long run distribution
  cmv.disc <- var.disc.mean.var( trans, pts )               # The conditional mean & variance
  lr.mu <- c( lr.d %*% pts )                                # The unconditional mean
  pts.dev <- ( pts - rep(1,nrow(pts)) %*% t(lr.mu) )        # Deviations from the mean
  cm.contrib <- matrix( apply(pts,1,sum), nrow(trans), ncol(trans), byrow=TRUE )
      # Because changes p(i,j) moves the conditional mean by X[1,] in each state
      # and so moves the sum by the sum of this row
  l.sq <- lapply( 1:nrow(trans), function(i) pts.dev[i,] %*% t(pts.dev[i,])  )
      # The squared realizations of the deviations
  cv.contrib <- matrix( sapply( 1:nrow(trans),
                                function(i) sum( l.sq[[i]][lower.tri(l.sq[[i]],TRUE) ] ) ),
                        nrow(trans), nrow(trans), byrow=TRUE )
      # Because the conditional variance changes with the LT f the l.sq list
  return(c(cm.contrib[,-nrow(trans)]+cv.contrib[,-nrow(trans)]))
}

min.disc.err <- function( l.var, n.pts, n.dirs, lb=NULL ){
# Computes the discretization that minimizes the square error on the conditional
# mean and variance of the Markov chain
  n.var <- nrow(l.var$A)
  lags <- ncol(l.var$A) / n.var
  X <- var.disc.pts( l.var, lags, n.pts, n.dirs, lb )
      # Create the matrix of points
  n.X <- nrow(X)
      # Number of points in discretization
  fn <- function(y){
    trans <- matrix( y, n.X, n.X-1 )
    trans <- cbind( trans, 1 - apply(trans,1,sum) )
        # Create the transition probabilities
    return( sum( .5 * var.disc.mean.var.err( trans, X, l.var$a, l.var$A, l.var$Sigma ) ^ 2 ) +
              1000 * any(trans < 0) * abs( min(trans) ) )
  }
  grad <- function(y){
    trans <- matrix( y, n.X, n.X-1 )
    trans <- cbind( trans, 1 - apply(trans,1,sum) )
        # Create the transition probabilities
    return(
              var.disc.mean.var.err.grad( trans, X, l.var$a, l.var$A, l.var$Sigma ) +
              1000 * any(trans < 0) * which.min(y) )
  }
  cm <- t( l.var$a + l.var$A %*% t( X ) )
      # The long-run variance
  M <- ls.min( t(X), t(cm) )
      # Guarantees that the mean is solved for
  while( any(M>1) || any(M<0) ){
    M <- pmin( pmax( t(M), 0 ), 1 )
    M <- M + ( 1 - apply(M,1,sum) ) / n.X
        # Make sure we have actual probabilities too
  }
  x.0 <- c( M[,-n.X] )
      # The initial guess
  x.ub <- rep( 1, n.X * (n.X-1) )
  x.lb <- rep( 0, n.X * (n.X-1) )
      # Parameter bounds
  opts <- list( algorithm='NLOPT_LD_SLSQP', maxeval=20000,
                print_level=1, check_derivatives=TRUE )
      # Optimization options
  optim <- nloptr( x.0, fn, grad, lb=x.lb, ub=x.ub, opts=opts )
      # Optomization
  trans <- cbind( matrix(optim$solution, n.X, n.X-1),
                  1 - apply( matrix(optim$solution, n.X, n.X-1), 1, sum) )
}

ls.min <- function(x, y) {
  # solves:  x %*% b = y
  d <- svd(x)
  # min-norm solution
  b.min <- d$v %*% diag(1/d$d, length(d$d)) %*% t(d$u) %*% y
  return(b.min)
}

disc.data.interp <- function( disc, dta, plot.on=TRUE, v.date=NULL, mains=NULL ){
# Interprets data as values from a VAR and plots the corresponding series
  n.vars <- ncol(dta)
  s.idx <- nn2( disc$X[,1:n.vars], dta, 1)$nn.idx
  dta.disc <- data.frame(disc$X[s.idx,])
      # Create the nearest neighbours

  if(plot.on){
    if(is.null(v.date)) v.date <- 1:nrow(dta.disc)
    if(is.null(mains)) mains <- c( paste0( 'Variable ', 1:ncol(disc$X) ),
                                   'Difference' )
    par.dft <- par('mar')
    par(mfrow=c(n.vars,1), mar=c(3,3,3,3))

    for( i in 1:n.vars ){
      plot(v.date, dta[,i], type='l', lwd=2, main=mains[i], xlab='', ylab='' )
      lines(v.date, dta.disc[,i], col='blue', lwd=2 )
      abline(h=0, lwd=.5)
    }
    par(mfrow=c(1,1), mar=par.dft)

    plot(v.date, dta[,2] - dta[,1], type='l', lwd=2, main=mains[3], xlab='', ylab='' )
    lines(v.date, dta.disc[,2] - dta.disc[,1], col='blue', lwd=2 )
    abline(h=0, lwd=.5)
  }
  return(list(s.idx=s.idx,dta.disc=dta.disc))
}

