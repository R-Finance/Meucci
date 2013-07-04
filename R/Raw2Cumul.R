#' Map raw moments into cumulative moments, as described in A. Meucci "Risk and Asset Allocation",
#' Springer, 2005
#'
#'	@param    mu_ : [vector] (length N corresponding to order N) corresponding raw moments
#'
#'	@return   ka  : [vector] (length N corresponding to order N) cumulative moments
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "Raw2Cumul.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export
Raw2Cumul = function( mu_ )
{
	N  = length( mu_ );
	ka = mu_;

	for( n in 2 : N )
	{
	    #ka[ n ] = mu_[ n ]; Doesn't make sense

	    for( k in 1 : (n-1) )
	    {
	        ka[ n ] = ka[ n ] - choose( n-1, k-1 ) * ka[ k ] * mu_[ n-k ]; 
	    }    
	}

	return( ka );
}