# TODO: add max weights constraint to EfficientFrontier()
# TODO: add computeCVaR to EfficientFrontier()
# TODO: confirm QuadProg does not have a bug (i.e. it can optimize expected returns without use dvec by adding an equality constraint)

#' Plots the efficient frontier, as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine,
#'  October 2008, p 100-106.
#'
#' @param  e  the NumPortf x 1 matrix of expected returns for each portfolio along the efficient frontier
#' @param  s  the NumPortf x 1 matrix of standard deviation of returns for each portfolio along the efficient frontier
#' @param  w  the NumPortf x N matrix of compositions (security weights) for each portfolio along the efficient frontier
#'   
#' @references 
#' A. Meucci, "Fully Flexible Views: Theory and Practice" \url{http://www.symmys.com/node/158}
#' See Meucci script for "RankingInformation/PlotFrontier.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export 


PlotFrontier = function( e, s, w )
{
  xx = dim( w )[ 1 ];
  N  = dim( w )[ 2 ];
  Data = t( apply( w, 1, cumsum ) );

  plot( c(min(s), 0), xlim = c( min(s) , max(s) ), ylim = c( 0, max(Data) ), 
    main= "frontier", xlab = " Portfolio # risk propensity", ylab = "Portfolio composition" );
  
  for( n in 1 : N )
  {
      x = rbind( min(s), s, max(s) );
      y = rbind( 0, matrix( Data[ , N-n+1 ] ), 0 );
      polygon( x, y, col = rgb( 0.9 - mod(n,3)*0.2, 0.9 - mod(n,3)*0.2, 0.9 - mod(n,3)*0.2) );
  }
}

#' Plots the results of computing the efficient frontier (Expected returns and frontier), as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine,
#' October 2008, p 100-106.
#'
#' @param  e      the NumPortf x 1 matrix of expected returns for each portfolio along the efficient frontier
#' @param  s      the NumPortf x 1 matrix of standard deviation of returns for each portfolio along the efficient frontier
#' @param  w      the NumPortf x N matrix of compositions (security weights) for each portfolio along the efficient frontier
#' @param  M      the NumPortf x 1 vector of expected returns for each asset
#' @param  Lower  constraints
#' @param  Upper  constraints 
#'   
#' @references 
#' A. Meucci, "Fully Flexible Views: Theory and Practice" \url{http://www.symmys.com/node/158}
#' See Meucci script for "RankingInformation/PlotResults.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

PlotResults = function( e, s, w, M, Lower = NULL , Upper = NULL)
{
  N = length( M );
  dev.new();
  par( mfrow = c( 1, 2 ) );
  h1 = hist( M*100, plot = F )
  barplot( h1$density, horiz = T, main = "expected returns", xlab = "", ylab = "" );
  if(length(Lower) || length(Upper))
  {
    Changed = array( 0, N );
    Changed[ union( Lower, Upper ) ] = M[ union( Lower, Upper ) ] * 100;
    h2 = hist(Changed, plot = F );
    barplot( h2$density, horiz = T, col = "red", add = T );
  }

  PlotFrontier( e*100, s*100, w );
}



#' Computes posterior probabilities to view the rankings, as it appears in A. Meucci,
#' "Fully Flexible Views: Theory and Practice", The Risk Magazine, October 2008, p 100-106.
#'
#' @param  X        a vector containing returns for all the asset classes
#' @param  p        a vector containing the prior probability values
#' @param  Lower    a vector of indexes indicating which column is lower than the corresponding column number in Upper
#' @param  Upper    a vector of indexes indicating which column is lower than the corresponding column number in Upper
#'
#' @references 
#' A. Meucci, "Fully Flexible Views: Theory and Practice" \url{http://www.symmys.com/node/158}
#' See Meucci script for "RankingInformation/ViewRanking.m"
#'
#' @author Ram Ahluwalia \email{ram@@wingedfootcapital.com}
#' @export EntropyProg

# example ViewRanking( X , p , Lower = c(3,4) , Upper = c(4,5) ) # two inequality views: asset 3 < asset 4 returns, and asset 4 < asset 5 returns

ViewRanking = function( X , p , Lower , Upper )
{
  library( matlab )
  J = nrow( X )
  N = ncol( X )
    
  K = length( Lower )
    
  # constrain probabilities to sum to one across all scenarios...
  Aeq = ones( 1 , J )
  beq = 1
    
  # ...constrain the expectations... A*x <= 0
  # X[,Lower] refers to the column of returns for Asset-lower
  # X[,Upper] refers to the column of returns for Asset-lower
  # X[ , Lower ] - X[ , Upper ] is vector returns of the "lower"" asset less the returns of the "higher" asset
  V = X[ , Lower ] - X[ , Upper ] # Jx1 vector. Expectation is assigned to each scenario
    
  A = t( V )
  b = 0 # The expectation is that (Lower - Upper)x <= 0. (i.e. The returns of upper are greater than zero for each scenario)
    
  # ...compute posterior probabilities
  p_ = EntropyProg( p , A , as.matrix(b) , Aeq , as.matrix(beq) )
    
  return( p_ )
}

#' Generates an efficient frontier based on Meucci's Ranking Information version and returns a A list with  
#' NumPortf efficient portfolios whos returns are equally spaced along the whole range of the efficient frontier,
#' as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine, October 2008, 
#' p 100-106.
#'
#' Most recent version of article and MATLAB code available at
#' http://www.symmys.com/node/158
#'
#' @param  X             a matrix with the joint-scenario probabilities by asset (rows are joint-scenarios, columns are assets)
#' @param  p             a vector of probabilities associated with each scenario in matrix X
#' @param  Options       a list of options....TBD
#'
#' @return Exps          the NumPortf x 1 vector of expected returns for each asset
#' @return Covs          the NumPortf x N vector of security volatilities along the efficient frontier
#' @return w             the NumPortf x N matrix of compositions (security weights) for each portfolio along the efficient frontier
#' @return e             the NumPortf x 1 matrix of expected returns for each portfolio along the efficient frontier
#' @return s             the NumPortf x 1 matrix of standard deviation of returns for each portfolio along the efficient frontier
#'
#' @references 
#' A. Meucci, "Fully Flexible Views: Theory and Practice" \url{http://www.symmys.com/node/158}
#' See Meucci script for "RankingInformation/EfficientFrontier.m"
#'
#' @author Ram Ahluwalia \email{ram@@wingedfootcapital.com} and Xavier Valls \email{flamejat@@gmail.com}
#' @export

RIEfficientFrontier = function( X , p , Options)
{ 

  if( !require("limSolve") ) stop("This script requieres the limSolve package installed")


  library( matlab )
    
  J = nrow( X ) # number of scenarios
  N = ncol( X ) # number of assets
    
  Exps = t(X) %*% p # probability-weighted expected return of each asset
    
  Scnd_Mom = t(X) %*% (X * ( p %*% matrix( 1, 1 , N ) ) )
  Scnd_Mom = ( Scnd_Mom + t(Scnd_Mom) ) / 2 # an N*N matrix
  Covs = Scnd_Mom - Exps %*% t( Exps )
    
  Constr = list()
    
  # constrain the sum of weights to 1
  Constr$Aeq = matrix( 1, 1 , N )
  Constr$beq = 1
    
  # constrain the weight of any security to between 0 and 1
  Constr$Aleq = rbind( diag( 1, N ) , - diag( 1, N ) ) # linear coefficients matrix A in the inequality constraint A*x <= b
  Constr$bleq = rbind( matrix( 1, N, 1 ) , matrix( 0, N, 1 ) ) # constraint vector b in the inequality constraint A*x <= b
    
  Amat = rbind( Constr$Aeq , Constr$Aleq ) # stack the equality constraints on top of the inequality constraints
  bvec = rbind( Constr$beq , Constr$bleq ) # stack the equality constraints on top of the inequality constraints
    
  ############################################################################################
  # determine return of minimum-risk portfolio
  FirstDegree  = matrix( 0, N , 1 ) # TODO: assumes that securities have zero expected returns when computing efficient frontier?
  SecondDegree = Covs
  # Why is FirstDegree "expected returns" set to 0? 
  # We capture the equality view in the equality constraints matrix
  # In other words, we have a constraint that the Expected Returns by Asset %*% Weights = Target Return
  MinVol_Weights = solve.QP( Dmat = SecondDegree , dvec = -1*FirstDegree , Amat = -1*t(Amat) , bvec = -1*bvec , meq = length( Constr$beq ) )
  MinSDev_Exp    = t( MinVol_Weights$solution ) %*% Exps
    
  ############################################################################################
  # determine return of maximum-return portfolio
  FirstDegree = -Exps
  MaxRet_Weights = linp( E = Constr$Aeq , F = Constr$beq , G = -1*Constr$Aleq , H = -1*Constr$bleq , Cost = FirstDegree , ispos = FALSE )$X
  MaxExp_Exp = t( MaxRet_Weights) %*% Exps
    
  ############################################################################################
  # slice efficient frontier in NumPortf equally thick horizontal sections
  Grid = matrix( , ncol = 0 , nrow = 0 )
  Grid = t( seq( from = Options$FrontierSpan[1] , to = Options$FrontierSpan[2] , length.out = Options$NumPortf ) )
    
  # the portfolio return varies from a minimum of MinSDev_Exp up to a maximum of MaxExp_Exp
  # We establish equally-spaced portfolio return targets and use this find efficient portfolios 
  # in the next step
  Targets = as.numeric( MinSDev_Exp ) + Grid * as.numeric( ( MaxExp_Exp - MinSDev_Exp ) ) 
    
  ############################################################################################
  # compute the NumPortf compositions and risk-return coordinates        
  FirstDegree = matrix( 0, N , 1 )
    
  w = matrix( , ncol = N , nrow = 0 )
  e = matrix( , ncol = 1 , nrow = 0 )
  s = matrix( , ncol = 1 , nrow = 0 )       
    
  for ( i in 1:Options$NumPortf )
  {
    # determine least risky portfolio for given expected return
    # Ax = b ; Exps %*% weights = Target Return
    AEq = rbind( Constr$Aeq , t( Exps ) )     # equality constraint: set expected return for each asset...
    bEq = rbind( Constr$beq , Targets[ i ] )  # ...and target portfolio return for i'th efficient portfolio
        
    Amat = rbind( AEq , Constr$Aleq ) # stack the equality constraints on top of the inequality constraints
    bvec = rbind( bEq , Constr$bleq )
        
    Weights = solve.QP( Dmat = SecondDegree , dvec = -1*FirstDegree , Amat = -1*t(Amat) , bvec = -1*bvec , meq = length( bEq ) )
        
    w = rbind( w , Weights$solution )
    s = rbind( s , sqrt( t(Weights$solution) %*% Covs %*% Weights$solution ) )
    e = rbind( e , Weights$solution %*% Exps )
  }
    
  return( list( e = e , Sdev = s , Composition = w , Exps = Exps , Covs = Covs ) )    
}
