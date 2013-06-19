library(pracma);

#' Pdf of the copula of the Student t distribution at the generic point u in the unit hypercube,
#' as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   u     : [vector] (J x 1) grade
#'	@param	 nu    : [numerical] 	  degrees of freedom 
#'	@param   Mu    : [vector] (N x 1) mean
#'	@param   Sigma : [matrix] (N x N) scatter
#'	
#'  
#'	@return   F_U   : [vector] (J x 1) PDF values
#'
#' @references
#' \url{http://}
#' See Meucci's script for "StudentTCopulaPdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

StudentTCopulaPdf = function( u, nu, Mu, Sigma )
{
	N = length( u );
	s = sqrt( diag( Sigma ));

	x = Mu + s * qt( u, nu);

	z2 = t(x - Mu) %*% mldivide( Sigma, (x - Mu)); #z2 = t(x - Mu) %*% inv(Sigma) * (x-Mu);
	K  = ( nu * pi )^( -N / 2 ) * gamma( ( nu + N ) / 2 ) / gamma( nu / 2 ) * ( ( det( Sigma ) )^( -0.5 ));
	Numerator = K * (1 + z2 / nu)^(-(nu + N) / 2);
	

	fs = dt((x - Mu) / s , nu);

	Denominator = prod(fs);

	F_U = Numerator / Denominator;

	return ( F_U );
}
