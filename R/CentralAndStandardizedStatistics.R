#' Compute central and standardized statistics, as described in A. Meucci 
#' "Risk and Asset Allocation", Springer, 2005
#'
#'	@param    X  : [vector] (J x 1) draws from the distribution
#'	@param    N  : [scalar] highest degree for the central moment
#'
#'	@return   ga : [vector] (1 x N) standardized statistics up to order N
#'	@return   mu : [vector] (1 x N) central moments up to order N
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "CentralAndStandardizedStatistics.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

CentralAndStandardizedStatistics = function( X, N )
{
	if(!require("PerformanceAnalytics")) stop("PerformanceAnalytics package required for this script");
	# compute central moments
	mu = matrix( 0, 1, N);
	mu[ 1 ] = mean(X);
	for( n in 2 : N )
	{
	    mu[ n ] = centeredmoment(X, n);
	}

	# compute standardized statistics 

	ga = mu;

	ga[ 2 ] = sqrt( mu[ 2 ]);
	for( n in 3 : N )
	{
	    ga[ n ] = mu[ n ] / (ga[ 2 ] ^ n);
	}
	
	return( list( ga = ga, mu = mu ) );

}