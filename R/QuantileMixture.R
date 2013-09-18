#' @title Computes the quantile of a mixture distribution by linear interpolation/extrapolation of the cdf.
#'
#' @description Computes the quantile of a mixture distribution by linear interpolation/extrapolation of the cdf. The confidence 
#' level p can be vector. If this vector is uniformly distributed on [0,1] the sample Q is distributed as the mixture.
#' Described in A. Meucci "Risk and Asset Allocation", Springer, 2005.
#'
#'  @param	p    [scalar] in [0,1], probability
#'	@param	a    [scalar] in (0,1), mixing probability
#'  @param	m_Y  [scalar] mean of normal component
#'  @param  s_Y  [scalar] standard deviation of normal component
#'  @param  m_Z  [scalar] first parameters of the log-normal component
#'  @param  s_Z  [scalar] second parameter of the log-normal component
#'  
#'  @return	Q    [scalar] quantile
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 184 - Estimation of a quantile of a mixture I".
#'
#'See Meucci's script for "QuantileMixture.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

QuantileMixture = function( p, a, m_Y, s_Y, m_Z, s_Z )
{
	# compute first moment
	m = a * m_Y + (1 - a) * exp( m_Z + 0.5 * s_Z * s_Z); 

	# compute second moment
	Ex2 = a * (m_Y^2 + s_Y^2) + (1 - a) * exp( 2 * m_Z + 2 * s_Z * s_Z);
	s   = sqrt( Ex2 - m * m );

	# compute cdf on suitable range
	X = m + 6 * s * seq( -1, 1, 0.001 );
	F = a * pnorm( X,  m_Y, s_Y) + (1 - a) * plnorm(X, m_Z, s_Z);
	X = X[!duplicated(F)];
	F = unique(F);

	# compute quantile by interpolation
	Q = interp1( F, X, p, method = "linear");

	return( Q );

}