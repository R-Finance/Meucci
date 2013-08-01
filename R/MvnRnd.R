#' Generate normal simulations whose sample moments match the population moments,
#' as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   M : [vector] (N x 1) expectation
#'	@param   S : [matrix] (N x N) covariance matrix
#'	@param   J : [scalar] number of draws (even number)
#'  
#'	@return  X : [matrix] (J x N) of drawsF_U   : [vector] (J x 1) PDF values
#'
#' @references
#' \url{http://symmys.com/node/170}, \url{http://www.symmys.com/node/162}{A. Meucci - "Simulations with Exact Means and Covariances", Risk, July 2009}
#' See Meucci's script for "MvnRnd.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com} and Ram Ahluwalia \email{rahluwalia@@gmail.com}
#' @export

MvnRnd = function( M, S, J )
{
	if ( !require( "QZ" ) ) stop("QZ package installation required for this script")
	N = length(M);

	# generate antithetic variables (mean = 0)
	Y = rmvnorm( J/2, matrix( 0, N ), S );
	Y = rbind( Y, -Y );

	# compute sample covariance: NOTE defined as "cov(Y,1)", not as "cov(Y)"
	S_ = ( dim(Y)[1] - 1 )/ dim(Y)[1] * cov( Y );

	# solve Riccati equation using Schur method
	H = rbind( cbind( matrix( 0, N, N ), -S ), cbind( -S, matrix( 0, N, N ) ) );
	 
	#Schur = Schur( H );
	#U = ordschur(U_,T_,'lhp');
	
	U = ordqz( H, keyword = "lhp" )$Q;

	U_lu = U[ 1:N, 1:N ];
	U_ld = U[ (N+1):nrow(U), 1:N ];

	B = U_ld %*% solve( U_lu );

	# affine transformation to match mean and covariances
	X = Y %*% B + kronecker( matrix( 1, J, 1 ),  t( M ) );

	return( X );
}