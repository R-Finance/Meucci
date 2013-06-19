#' Compute the mean and standard deviation of a lognormal distribution from its parameters, as described in  
#' A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'
#'	@param  e    : [scalar] expected value of the lognormal distribution
#'  @param	v    : [scalar] variance of the lognormal distribution
#'  
#'  @return	mu   : [scalar] expected value of the normal distribution
#'  @return	sig2 : [scalar] variance of the normal distribution
#'  
#'  @note	Inverts the formulas (1.98)-(1.99) in "Risk and Asset Allocation", Springer, 2005.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "LognormalMoments2Parameters.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

LognormalMoments2Parameters = function( e, v)
{
	sig2 = log( 1 + v / ( e^2 ) );
	mu = log( e ) - sig2 / 2;
	
	return( list( sigma_square = sig2 , mu = mu ) );

}