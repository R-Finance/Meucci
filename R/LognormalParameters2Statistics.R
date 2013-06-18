#' Compute expectation, Cov, standard deviation and Corr for a lognormal distribution, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'	@param Mu    : [vector] (N x 1) location parameter
#'	@param Sigma : [matrix] (N x N) scale parameter
#'
#'  
#'	@return Exp   : [vector] (N x 1) expectation
#'	@return Cov   : [matrix] (N x N) covariance
#'	@return Std   : [vector] (N x 1) standard deviation
#'	@return	Corr  : [matrix] (N x N) correlation
#'
#' @references
#' \url{http://}
#' See Meucci's script for "LognormalParam2Statistics.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

LognormalParam2Statistics = function(Mu, Sigma)
{

	Exp = exp( Mu + (1/2) * diag( Sigma ) );
	Cov = exp( Mu + (1/2) * diag( Sigma ) ) %*% t( exp( Mu + (1/2) * diag( Sigma ) ) ) * ( exp( Sigma ) - 1 );
	Std = sqrt( diag( Cov ) );
	Corr = diag( 1 / Std ) %*% Cov %*% diag( 1 / Std );

	return( list( Exp = Exp, Covariance = Cov, Standard_Deviation = Std, Correlation = Corr ));
}