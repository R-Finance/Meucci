#' @title Computes the pdf of the copula of the normal distribution at the generic point u in the unit hypercube
#'
#' @description Computes the pdf of the copula of the normal distribution at the generic point u in the unit 
#' hypercube, as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   u      [vector] (J x 1) grade
#'	@param   Mu     [vector] (N x 1) mean
#'	@param   Sigma  [matrix] (N x N) covariance
#'  
#'	@return  F_U    [vector] (J x 1) PDF values
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, 
#' "E 33 - Pdf of the normal copula".
#'
#' See Meucci's script for "NormalCopulaPdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

NormalCopulaPdf = function( u, Mu, Sigma )
{
	N = length( u );
	s = sqrt( diag( Sigma ));

	x = qnorm( u, Mu, s );

	Numerator = ( 2 * pi ) ^ ( -N / 2 ) * ( (det ( Sigma ) ) ^ ( -0.5 ) ) * exp( -0.5 * t(x - Mu) %*% mldivide( Sigma , ( x  - Mu )));

	fs = dnorm( x, Mu, s);

	Denominator = prod(fs);

	F_U = Numerator / Denominator;

	return ( F_U );
}