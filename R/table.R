#####################################################################
# table.R
#
# Contains the functions for creating regression tables for VARs
#
# 23may2017
# Philip Barrett, Washington DC
#####################################################################


var.table <- function( l.l.var, l.m.var=NULL, file=NULL, varnames=NULL,
                       specnames = NULL, label=NULL, caption=NULL, footer=TRUE,
                       add.mean=TRUE, start=NULL, end=NULL, v.lhood=NULL,
                       font.size=NULL){
# Creates a nice regression table for the VARs

  ### Set up ###
  add.se <- !is.null( l.m.var )                     # Indicator for adding the mean
  n.spec <- length(l.l.var)                         # Number of specifications
  n.v <- nrow(l.l.var[[1]]$A)                       # Number of variables
  lags <- ncol(l.l.var[[1]]$A) / n.v                # Number of lags
  n.col <- 1 + n.v * n.spec                         # Number of columns
  if( is.null(varnames) ) varnames <- paste0( 'Y', 1:n.v )
  if( is.null(specnames) ) specnames <- paste0( '(', 1:n.spec, ')' )
  l.mu <- lapply( l.l.var, function(x) mu.calc( x$a, x$A ) )
      # The mean (in case not provided)
  st.font.size <- if( is.null(font.size) ) NULL else paste0('\n\t \\', font.size)

  ###  Create top and tail ###
  head.str <- paste0( '\\begin{table}[htbp] \n\t\\centering', st.font.size,
                      '\n\t\\begin{tabular}{@{\\extracolsep{4pt}}',
                      paste( c( 'l', rep('c', n.col-1) ), collapse ='' ), '@{}}' )
  header.spec <- paste0( '\t\t\\hline\\hline\n\t\t \t\t & ',
                    paste0( sapply( specnames,
                            function(x) paste0( '\\multicolumn{',
                                                n.v, '}{c}{',x,'} ' ) ), collapse='&' ),
                    '\\\\ \n' )
  header.var <- paste0( ' \t\t & ', paste0( rep( varnames, n.spec ), collapse=' \t & ' ),
                        '\\\\',
                paste0( sapply( 1:n.spec,
                    function(i) paste0( '\\cline{', 2+(i-1)*n.v,
                                        '-', 1+i*n.v, '}' ) ), collapse='' ), '\n' )
  header <- paste0( head.str, header.spec, header.var )
      # The header
  tail <- paste0( '\n \\hline \t\\end{tabular}',
                  if(is.null(caption)) '' else paste0('\t\t\\caption{',caption,'}\n'),
                  if(is.null(label)) '' else paste0('\t\t\\label{',label,'}\n'),
                  '\n\\end{table}' )
      # The closing part

  ## Create the standard devations
  n.coeff <- length( par.to( l.l.var[[1]]$a, l.l.var[[1]]$A, l.l.var[[1]]$Sigma ) )
      # Number of coefficients
  if( add.se ){
    se <- sapply( l.m.var, function(x) sqrt(diag(x))[1:n.coeff] )
        # Creates NAs where variances are not supplied (e.g. for covars)
  }
  se.par <- lapply( 1:n.spec, function(i) par.from( se[,i], n.v, lags ) )
      # Creates a list form for the standard errors

  if( add.mean & add.se ){
    se.mu <- sapply( 1:n.spec,
                     function(i){
                       n.a.A <- (n.v*(1+n.v*lags))
                       vcv <- l.m.var[[i]][1:n.a.A, 1:n.a.A]
                          # Select only the a & A parameters, as mu.calc.d
                          # doesn't return for Sigma.
                       g.h <- mu.calc.d( l.l.var[[i]]$a, l.l.var[[i]]$A, l.l.var[[i]]$Sigma )
                       return( sqrt( diag( t(g.h) %*% vcv %*% g.h ) ) )
                     }  )
    # The se of the means: Computed via numerical delta method.
  }

  ### Constant rows ###
  coeff.title <- paste0( '\\rule{0pt}{4ex} \n \\emph{Coefficients} \t ',
                         paste0( rep(' \t\t &', n.col-1 ), collapse='' ), '\\\\ \n' )
  const <- paste0( '\\quad Constant \t & ',
                   paste0( c( sapply(l.l.var, function(x) round( x$a, 2 ) ) ),
                           collapse=' \t & ' ), '\t \\\\ \n' )
  if( add.se )
    const <- paste0( const, ' \t\t & (',
                     paste0( c( sapply(se.par, function(x) round( x$a, 3) ) ),
                             collapse=') \t & (' ), ') \t \\\\ \n' )

  ### Mean rows ###
  if( add.mean ){
    mn <- paste0( '\\quad LR Mean \t & ',
                   paste0( c( sapply(l.mu, function(x) round( x, 2 ) ) ),
                           collapse=' \t & ' ), '\t \\\\ \n' )
    if( add.se )
      mn <- paste0( mn, ' \t\t & (',
                    paste0( c( round(se.mu, 3) ), collapse=') \t & (' ),
                    ') \t \\\\ \n' )
  }else{
    mn <- NULL
  }

  ### The dynamic coefficients ###
  dyn.coeff <- NULL
  for( i in 1:lags ){
    for( j in 1:n.v ){
      dyn.coeff <- paste0( dyn.coeff, '\\quad ', varnames[j], ' (-', i, ') \t &',
                           paste0( c( sapply(l.l.var, function(x) round( x$A[,(i-1)*n.v+j], 2 ) ) ),
                                   collapse=' \t & ' ), '\t \\\\ \n' )
      if( add.se )
        dyn.coeff <- paste0( dyn.coeff, ' \t\t & (',
                         paste0( c( sapply(se.par, function(x) round( x$A[,(i-1)*n.v+j], 3) ) ),
                                 collapse=') \t & (' ), ') \t \\\\ \n' )
    }
  }

  coeff <- paste0( coeff.title, const, mn, dyn.coeff )

  ### The covariances ###
  cov.title <- paste0( '\\rule{0pt}{4ex} \\emph{Innov. covar.} ',
                         paste0( rep(' \t &', n.col-1 ), collapse='' ), '\\\\ \n' )
  idx <- lower.tri(l.l.var[[1]]$Sigma, TRUE )
  l.cov.st <- lapply( l.l.var, function(x) toString(x$Sigma) )
  cov.coeff <- NULL
  for( i in 1:n.v ){
    cov.coeff <- paste0( cov.coeff, '\\quad ', varnames[i], ' \t &',
                         paste0( c( sapply(l.l.var, function(x)
                           sapply(1:n.v, function(j)
                             if(idx[i,j]) round( x$Sigma[i,j], 2 ) else '' ) ) ),
                                 collapse=' \t & ' ), '\t \\\\ \n' )
    if( add.se )
      cov.coeff <- paste0( cov.coeff, ' \t\t & ',
                           paste0( c( sapply(se.par, function(x)
                             sapply( 1:n.v, function(j)
                                if(idx[i,j] & !is.na(x$Sigma[i,j]) & x$Sigma[i,j] != 0 )
                                  paste0( '(', round( x$Sigma[i,j], 3),')' ) else '' ) ) ),
                                    collapse=' \t & ' ), ' \t \\\\ \n' )
  }
  st.cov <- paste0( cov.title, cov.coeff )

  if(footer){
    st.foot <- ' \\hline \\rule{0pt}{4ex} \n  '
    if(!is.null(v.lhood))
      st.foot <- paste0( st.foot, 'Log likelihood \t &',
                            paste0( sapply( v.lhood,
                              function(x) paste0( '\\multicolumn{', n.v, '}{c}{',
                                                  round(x,1), '}' ) ),
                                    collapse = ' \t & '), '\\\\ \n' )
    ## To add: R sq, adj R sq, start and end dates
  }

  #   if(footer & add.se){
#     if( !is.null(start) & !is.null(end) ){
#       st.dates <- paste0( ' ', year(start), 'Q', quarter(start), ':',
#                           year(end), 'Q', quarter(end) )
#     }else{
#       st.dates <- NULL
#     }
#     p.val <- 1 - pnorm( 0, diff(mu), se.rmg(l.var) )
#     # p.value of R-G<0 hypothesis test
#     st.foot <- paste0( '\t\t\\hline\n\t\t\\multicolumn{', n.col, '}{l}{', nrow(l.var$datamat),
#                        ' obs', st.dates , '\\hfill $H_0: \\mathbb E R = \\mathbb E G$ vs. $H_0: \\mathbb E R < \\mathbb E G$ p=',
#                        round( p.val, 3 ), '}\\\\\n' )
#   }else{
#     st.foot <- NULL
#   }

  # out <- paste0( head.str, header, bdy, st.foot, tail.str )
  out <- paste0( header, coeff, st.cov, st.foot, tail )

  if(!is.null(file)) cat(out, file=file)

  return(out)
}
