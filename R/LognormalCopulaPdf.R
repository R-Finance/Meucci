#' @title Computes the pdf of the copula of the lognormal distribution at the generic point u in the unit hypercube. 
#'
#' @description Computes the pdf of the copula of the lognormal distribution at the generic point u in the unit hypercube,
#' as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   u      [vector] (J x 1) grades
#'	@param   Mu     [vector] (N x 1) location parameter
#'	@param   Sigma  [matrix] (N x N) scatter parameter
#'  
#'	@return  F_U    [vector] (J x 1) PDF values
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, 
#' "E 36 - Pdf of the lognormal copula".
#'
#' See Meucci's script for "LognormalCopulaPdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

LognormalCopulaPdf = function( u, Mu, Sigma )
{	
	N = length( u );
	s = sqrt( diag( Sigma ));

	x = qlnorm( u, Mu, s );

	Numerator = ( 2 * pi ) ^ ( -N / 2 ) * ( (det ( Sigma ) ) ^ ( -0.5 ) ) / 
					prod(x) * exp( -0.5 * t(log(x) - Mu) %*% mldivide( Sigma , ( log( x ) - Mu ), pinv=FALSE ) );

	fs = dlnorm( x, Mu, s);

	Denominator = prod(fs);

	F_U = Numerator / Denominator;

	return ( F_U );
}