#####################################################################
# table.R
#
# Contains the functions for creating regression tables for VARs
#
# 23may2017
# Philip Barrett, Washington DC
#####################################################################



var.table <- function( l.var, file=NULL, add.mean=FALSE, varnames=NULL,
                       label=NULL, caption=NULL, add.se=TRUE, footer=TRUE,
                       start=NULL, end=NULL ){
# Creates a nice regression table for the VAR

  ### Set up ###
  n.v <- length(l.var$varresult)      # Number of variables
  n.col <- 2 + n.v * 2 + add.mean     # Number of columns
  if( is.null(varnames) ) varnames <-
    gsub( '_', '', names(l.var$varresult[[1]]$coefficients[-(n.v+1)]) )
  # coeff <- t( sapply( l.var$varresult, function(x) x$coefficients ) )
  a <- l.var$a
  A <- l.var$A
  Sigma <- l.var$Sigma
      # Coefficients
  mu <- mu.calc( a, A )
      # The mean

  if( add.se ){
    ll <- summary(l.var)                # Need for the covariance estimates


    covres <- ll$covres
        # The covariance of the residuals
    se <- t( sapply( ll$varresult, function(x) x$coefficients[,'Std. Error'] ) )
    # The standard error
    vcv <- kronecker( covres, ll$varresult[[1]]$cov.unscaled )
    # The variance-covariance matrix of the estimates.  The sqrt of the diag
    # of this will be se.  We can do this with only the first regression in ll
    # because the data is the same for all of the (as it is a VAR).
    g.h <- grad.mu( coeff[,-(n.v+1)], coeff[,n.v+1] )
    se.mu <- sqrt( diag( t(g.h) %*% vcv %*% g.h ) )
    # The se of the means
  }else{
    covres <- cov( sapply( sim.VAR$varresult, function(x) x$residuals ) )
    # The covariance of the residuals
  }

  # Create top and tail
  head.str <- paste0( '\\begin{table}[htbp] \n\t\\centering \n\t\\begin{tabular}{',
                      paste(rep('c', n.col), collapse =''),'}' )
  header.coeff <- paste0( '\t\t\\hline\\hline\n\t\t \t\t & \\multicolumn{', n.v + 1,
                          '}{c}{Regression coefficients} &' )
  header.mean <- if(!add.mean) '' else paste0( 'Mean &' )
  header.covar <- paste0( '\\multicolumn{', n.v, '}{c}{Covariances} \\\\ \n' )
  coeff.names <- paste0( '\t\t \t\t & Const \t & ',
                         paste0( varnames, collapse=' (-1) \t & ' ), ' (-1)',
                         (if(add.mean) '\t & \t & \t' else '\t & \t' ),
                         paste0( varnames, collapse=' \t & ' ), '\\\\ \\hline \n')
  header <- paste0( header.coeff, header.mean, header.covar, coeff.names )
  # The header
  tail.str <- paste0( '\t \\hline \n \t\\end{tabular}',
                      if(is.null(caption)) '' else paste0('\t\t\\caption{',caption,'}\n'),
                      if(is.null(label)) '' else paste0('\t\t\\label{',label,'}\n'),
                      '\n\\end{table}' )
  # The closing part

  # Create the body
  line.i <- function(i){
    # Create the line(s) in the table for variable number i
    nm <- paste0( varnames[i], '\t &')
    coeff.est <- paste0( round( coeff[i,c( n.v+1, 1:n.v )], 2),
                         collapse = '\t &' )
    if( add.mean ) coeff.est <- paste0( coeff.est, ' \t & ', round( mu[i], 2) )
    sig.est <- paste0( sapply( 1:n.v, function(j) if( i >= j ) round( covres[i,j], 3) else ' \t ' ),
                       collapse = '\t &' )
    if(add.se){
      coeff.se <- paste0( round( se[i,c( n.v+1, 1:n.v )], 3) , collapse = ') \t & (' )
      sig.se <- paste0( rep( ' ', n.v ), collapse = ' \t & ' )
      if( add.mean ) coeff.se <- paste0( coeff.se, ') \t & (', round(se.mu[i],3) )
      se.line <- paste0( '\t & (', coeff.se, ') \t & ', sig.se, ' \\\\\n' )
    }else{
      se.line <- NULL
    }
    return( paste0( nm, coeff.est, '\t &', sig.est, '\\\\\n', se.line ) )
  }
  bdy <- paste0( sapply( 1:n.v, line.i  ), collapse ='\n' )

  if(footer & add.se){
    if( !is.null(start) & !is.null(end) ){
      st.dates <- paste0( ' ', year(start), 'Q', quarter(start), ':',
                          year(end), 'Q', quarter(end) )
    }else{
      st.dates <- NULL
    }
    p.val <- 1 - pnorm( 0, diff(mu), se.rmg(l.var) )
    # p.value of R-G<0 hypothesis test
    st.foot <- paste0( '\t\t\\hline\n\t\t\\multicolumn{', n.col, '}{l}{', nrow(l.var$datamat),
                       ' obs', st.dates , '\\hfill $H_0: \\mathbb E R = \\mathbb E G$ vs. $H_0: \\mathbb E R < \\mathbb E G$ p=',
                       round( p.val, 3 ), '}\\\\\n' )
  }else{
    st.foot <- NULL
  }

  out <- paste0( head.str, header, bdy, st.foot, tail.str )

  if(!is.null(file)) cat(out, file=file)

  return(out)
}
