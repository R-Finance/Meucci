#' Computes a double-decay covariance matrix.
#'
#' This function computes a double-decay covariance matrix for the risk drivers provided, as described in  
#' A. Meucci, "Personalized Risk Management: Historical Scenarios with Fully Flexible Probabilities"
#' GARP Risk Professional, Dec 2010, p 47-51
#' 
#' @param   X       matrix representing the risk drivers.
#' @param   lmd_c   numeric representing the low decay (long half-life) for the correlations.
#' @param   lmd_s   numeric representing the high decay (short half-life) for the volatilities.
#' @return  m       matrix of zeros, representing the expectation of the risk drivers.
#' @return  S       matrix representing the double-decay estimation for the correlation matrix of the risk drivers.
#' 
#' @references 
#' \url{http://www.symmys.com/node/150}
#' See Meucci script for "DoubleDecay.m"
#' 
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

DoubleDecay = function( X, lmd_c, lmd_s)
{

	N = dim( X )
	m = matrix( 0, N[2], 1);

	p_c = exp( -lmd_c * ( N[1] - t(rbind( 1:N[1] ) ) ) );

	p_c = kronecker( matrix( 1, 1, N[2] ), p_c / sum( p_c ) ); 	#  workaround on p_c=repmat( p_c/sum(p_c),1,N);

	S_1 = t( p_c * X ) %*% X;

	C = cov2cor( S_1 );

	p_s = exp( -lmd_s * ( N[1] - t(rbind( 1:N[1] ) ) ) );

	p_s = kronecker(matrix(1,1,N[2]),p_s/sum(p_s));
	S_2 = t( p_s*X ) %*% X;

	R = cov2cor(S_2) ;  
	s = c( 0.0099, 0.0538, 0.0163 );

	s = sqrt(diag(S_2));
	S = diag(s) %*% C %*% diag(s);

	return( list( m = m , S = S ) )
}