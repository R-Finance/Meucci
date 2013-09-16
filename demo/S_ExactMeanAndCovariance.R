#' Generate draws from a multivariate normal with matching mean and covariance, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_ExactMeanAndCovariance.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}


########################################################################################################
### Inputs
N = 20;  # dimension (number of risk factors)
J = 200; # number of simulations

########################################################################################################
### Generate desired population moments:

# vector of expected values M
M = matrix(runif( N ) -0.5);
# covariance matrix S
A  = matrix( runif( N * N ), c( N, N )) - 0.5;
S = A %*% t( A );

# generate sample of size J from multivariate normal N(M,S)
#X = mvnrnd(M, S, J); # no match between sample and population moments (built-in) function
X = MvnRnd( M, S, J ); # exact match between sample and population moments

########################################################################################################
### Compute sample moments and errors
M_ = matrix( apply( X, 2, mean )); #apply
S_ = ( dim( X )[1] - 1 )/ dim( X )[1] * cov( X );

########################################################################################################
### Check errors
Err_M = max( abs( M - M_ ) ) / max( abs( M ) );
Err_S = max( max( abs( S - S_) ) )/ max( max( abs( S ) ) );

print(Err_M);
print(Err_S);

