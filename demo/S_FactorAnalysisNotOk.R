#'This script illustrates the hidden factor analysis puzzle, as described in A. Meucci,
#'"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 111 - Hidden factors: puzzle".
#'
#' See Meucci's script for "S_FactorAnalysisNotOk.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs

N = 5; # market dimension
K = 2; # factors dimension
J = 10000; # numbers of simulations

##################################################################################################################
### Define true hidden loadings B 

B = matrix( runif( N*K ), N, K ) - 0.5; 

B = B / sqrt( 1.5 * max( max( B %*% t(B) ) ) );

# define true hidden idiosyncratic variances 
D = array( 1, N) - diag( B %*% t(B) );

# define true hidden global covariance
S = B %*% t( B ) + diag( D, length(D) );

# generate normal variables with matching moments
X = MvnRnd( matrix( 0, N, 1 ), S, J );

# recover loadings FA$loadings, idiosyncratic variances FA$uniquenesness and factors FA$scores by factor analysis
#[FA$loadings, FA$uniquenesness, T_, stats, F_] 
FA = factanal(X, K, scores = "Bartlett" );

# factor analysis recovers the structure exactly however...
S_    = FA$loadings %*% t( FA$loadings ) + diag( FA$uniquenesses, length( FA$uniquenesses) );
Match = 1 - max( abs( ( S - S_) / S) );
print(Match);

# ...the systematic+idiosyncratic decomposition is NOT recovered
U_  = X - FA$scores %*% t(FA$loadings); # compute residuals
S_U = cor( U_ ); # compute correlations

# residuals are not idiosyncratic
print( S_U );